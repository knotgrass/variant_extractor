import logging
from collections import defaultdict, namedtuple
from lib.pubSeqTables import threeToOneLower#, threeToOne, oneToThree, aaToDna, dnaToAa
# from lib.pm_pycbio.hgdata.Psl import Psl
from lib import pubAlg#, maxCommon
from pubMunch import maxCommon
# from pygr.seqdb import SequenceFileDB
# from pyfaidx import Fasta
# try:                import re2 as re
# except ImportError: import re
import re, sys
from rich import print


# this setting can be changed to allow protein variants
# that require a change of two base pairs. By default, it
# is off to reduce false positives
# Johannes: going for high recall, not high precision
# allowTwoBpVariants = False

# ===== DATA TYPES ========
Mention = namedtuple("Mention", "patName,start,end")

""" A mapped variant is a type-range-sequence combination from a text, 
    can be located on none, one or multiple types of sequences
    All positions are 0-based
"""
VariantFields = [
    "mutType",  # sub, del or ins
    "seqType",  # prot, dna or dbSnp
    "seqId",  # protein or nucleotide accession
    "geneId",  # entrez gene id, if a gene was found nearby in text
    "start",  # original position in text
    "end",  # end position in text
    "origSeq",  # wild type seq, used for sub and del, also used to store rsId for dbSnp variants
    "mutSeq",  # mutated seq, used for sub and ins
    "origStr",  # variant name as it occurs in the original paper
    "offset",  # offset for splicing variants, e.g., -3 for c.1184-3A>T
    "firstAa",  # first amino acid in a notation like p.A100_F102del2
    "secondAa",  # second amino acid in a notation like p.A100_F102del2
    "length",  # length of inserted or deleted string such as in p.A100_F102del2 (has length 2)
    "docId",  # document in which variant occurs (usually a pubmed ID)
    "variantType",  # missense, nonframeshift, frameshift, stopgain, stoploss, splicing
    "ivsNumber",  # IVS splicing notation number
    "vcfPos",
    "vcfRef",
    "vcfAlt",
    "correctProba"
    ]

# A completely resolved mutation
mutFields = (
    "chrom",  # chromosome
    "start",  # on chrom
    "end",  # on chrom
    "varId",  # a unique id
    "inDb",  # list of db names where was found
    "patType",  # the type of the patterns (sub, del, ins)
    "hgvsProt",  # hgvs on protein, can be multiple, separated with |
    "hgvsCoding",  # hgvs on cdna, can be multiple, separated with |
    "hgvsRna",  # hgvs on refseq, separated by "|"
    "comment",  # comment on how mapping was done
    "rsIds",  # possible rsIds, separated by "|", obtained by mapping from hgvsRna
    "protId",  # the protein ID that was used for the first mapping
    "texts",  # mutation match in text
    # "mutSupport",  # prot, dna, protDna
    # "mutCount",    # how often was this mutation mentioned?
    "rsIdsMentioned",  # mentioned dbSnp IDs that support any of the hgvsRna mutations
    "dbSnpStarts" ,  # mentioned dbSnp IDs in text, start positions
    "dbSnpEnds",  # mentioned dbSNP Ids in text, end positions

    "geneSymbol",  # symbol of gene
    "geneType",  # why was this gene selected (entrez, symNearby, symInTitle, symInAbstract)
    "entrezId",  # entrez ID of gene
    "geneStarts",  # start positions of gene mentions in document
    "geneEnds",  # end positions of gene mentions in document

    "seqType",  # the seqType of the patterns, dna or protein
    "mutPatNames",  # the names of the patterns that matched, separated by |
    "mutStarts",  # start positions of mutation pattern matches in document
    "mutEnds",  # end positions of mutation pattern matches in document
    "mutSnippets",  # the phrases around the mutation mentions, separated by "|"
    "geneSnippets",  # the phrases around the gene mentions, separated by "|"
    "dbSnpSnippets"  # mentioned dbSNP Ids in text, snippets
    )

# fields of the output file
MutRec = namedtuple("mutation_desc", mutFields)

# ======= GLOBALS ===============
# this can be used to shuffle all protein sequences before matching
# to get a random background estimate
doShuffle = False

geneData = None


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.info("Loading gene information for varFinder.py")


def parseRegex(regex_path:str = 'data/regex.txt'):
    """ parse and compile regexes to list (seqType, mutType, patName, pat) """
    # read regexes, translate placeholders to long form and compile
    replDict = {
    "sep"         : r"""(?:^|[:;\s\(\[\'\"/,])""",
    "fromPos"     : r'(?P<fromPos>[1-9][0-9]*)',
    "toPos"       : r'(?P<toPos>[1-9][0-9]*)',
    "pos"         : r'(?P<pos>[1-9][0-9]*)',
    "offset"         : r'(?P<offset>[1-9][0-9]*)',
    "fromPoss"     : r'(?P<fromPos>[1-9][0-9]+)',
    "toPoss"       : r'(?P<toPos>[1-9][0-9]+)',
    "poss"         : r'(?P<pos>[1-9][0-9]+)',
    "plusMinus"    : r'(?P<plusMinus>[+\^\-\u00FE\u2AF9\u2AFA])',  # plus often gets picked up as \xc3\xbe, whatever that is
    "arrow"        : r'(?P<arrow>[/r\>4\u2192\.!\u2B0E\u02DA\[])',  # u2192 is the right arrow in unicode; some people choose to write c.234G.A instead of c.234G>A, maybe their shift key is not working?
    "underscore"   : r'(?P<underscore>[_ ])',
    "origAaShort" : r'(?P<origAaShort>[CISQMNPKDTFAGHLRWVEYX])',
    "origAasShort" : r'(?P<origAasShort>[CISQMNPKDTFAGHLRWVEYX]+)',
    "origAaLong"  : r'(?P<origAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X))',
    "origAasLong"  : r'(?P<origAasLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X)+)',
    "mutAaShort"  : r'(?P<mutAaShort>[fCISQMNPKDTFAGHLRWVEYX*])',
    "mutAaLong"  : r'(?P<mutAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X))',
    "firstAaShort"  : r'(?P<firstAaShort>[CISQMNPKDTFAGHLRWVEY])',
    "firstAaLong"  : r'(?P<firstAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE))',
    "secondAaShort"  : r'(?P<secondAaShort>[CISQMNPKDTFAGHLRWVEY])',
    "secondAaLong"  : r'(?P<secondAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE))',
    "mutAasShort"  : r'(?P<mutAasShort>[fCISQMNPKDTFAGHLRWVEYX*]+)',
    "mutAasLong"  : r'(?P<mutAasLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X)+)',
    "dna"         : r'(?P<dna>[ACTG])',
    "dnas"         : r'(?P<dnas>[ACTG]+)',
    "noDnas"         : r'([ACTG]+)',
    "origDna"     : r'(?P<origDna>[ACTG])',
    "origDnas"     : r'(?P<origDnas>[ACTG]+)',
    "mutDna"      : r'(?P<mutDna>[ACTG])',
    "mutDnas"      : r'(?P<mutDnas>[ACTG]+)',
    "fs"          : r'(?P<fs>fs([\*\u2217\u204E\u26B9\u2731\u066D]?[0-9]*[*\u2217\u204E\u26B9\u2731\u066D]?))',
    "length"      : r'(?P<length>[1-9][0-9]*)',
    "space"       : r'([ ?]{1,3})',
    "ivsNumber"   : r'(?P<ivsNumber>([IV]+)|[0-9]+)'
    }
    

    logger.info("Parsing regexes from %s" % regex_path)
    regexList = []
    counts = defaultdict(int)
    for row in maxCommon.iterTsvRows(regex_path, commentPrefix="#"):
        logger.info("Translating %s" % row.pat)
        patName = row.patName
        if patName == "":
            patName = row.pat
        patFull = row.pat.format(**replDict)
        logger.info("full pattern is %s" % patFull)
        flags = 0
        if "Long}" in row.pat:
            flags = re.IGNORECASE
            logger.info("ignoring case for this pattern")
        # print("Full pattern: %s" % patFull)
        patComp = re.compile(patFull, flags=flags | re.UNICODE)
        regexList.append((row.seqType, row.mutType, patName, patComp))
        counts[(row.seqType, row.mutType)] += 1

    for regexType, count in counts.items():
        logger.info("regexType %s, found %d regexes" % (str(regexType), count))
    return regexList


def makeMention(match, patName):
    start = match.start()
    end = match.end()
    mention = Mention(patName, start, end)
    return mention


class VariantDescription(object):
    """ A variant description fully describes a variant
        It has at least a type (sub, del, etc), a start-end position on a
        (potentially unknown) sequence a tuple (origSeq, mutSeq) that describes
        the mutation e.g. ("R", "S"). The class can generate a descriptive name
        for the mutation like "p.R123S"

        It can optionally include a sequence ID, when the sequence ID was part of the 
        mutation description in the text, e.g. the HGVS "NM_000925:p.R123S"

    >>> VariantDescription("sub", "prot", 10, 11, "R", "S")
    VariantDescription(mutType=u'sub',seqType=u'prot',seqId=u'None',geneId=u'',start=u'10',end=u'11',origSeq=u'R',mutSeq=u'S')
    """
    __slots__ = VariantFields

    def __init__(self, mutType, seqType, start, end, origSeq, mutSeq, docId,
                 seqId=None, geneId="", offset=0, origStr="", firstAa=None, secondAa=None,
                 length=None, variantType=None, ivsNumber=None,
                 vcfPos=None, vcfRef=None, vcfAlt=None, correctProba=0.5):
        self.mutType = mutType  # sub, del, ins, dbSnp, ...
        self.seqType = seqType  # cds, rna, dna, prot,  ...
        self.seqId = seqId  # NM_*                  ...
        self.geneId = geneId  #                       ...
        self.start = int(start) if start is not None else None  #                       ...
        self.end = int(end)  if end is not None else None  #                       ...
        self.origSeq = origSeq  #                       ...
        self.mutSeq = mutSeq  #                       ...
        self.origStr = origStr  #                       ...
        self.offset = offset  #                       ...
        self.firstAa = firstAa  #                       ...
        self.secondAa = secondAa  #                       ...
        if length == None and self.origSeq:
            self.length = len(self.origSeq)
        else:
            self.length = length
        self.docId = docId
        self.variantType = variantType
        self.ivsNumber = ivsNumber
        self.vcfPos = vcfPos
        self.vcfRef = vcfRef
        self.vcfAlt = vcfAlt
        self.correctProba = correctProba

    def getName(self):
        # used to be HGVS text for variant; now just unique identifier
        return str(self)

    def asRow(self):
        row = []
        for i in self.__slots__:
            row.append(str(getattr(self, i)))
        return row

    def __repr__(self):
        # return ",".join(self.asRow())
        parts = []
        for field in self.__slots__:
            parts.append(field + "=" + repr(str(getattr(self, field))))
        return "VariantDescription(%s)" % ",".join(parts)


def isOverlapping(match, exclPos):
    posSet = set(range(match.start(), match.end()))
    if len(posSet.intersection(exclPos)) != 0:
        logger.debug("regex overlaps an excluded position (gene?)")
        return True
    return False


def parseMatchSub(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    # grab long and short versions of amino acid
    correctProba = 0.4
    if "origAaShort" in groups and "mutAaShort" in groups:
        if groups["origAaShort"].isupper() and groups["mutAaShort"].islower() or \
           groups["origAaShort"].islower() and groups["mutAaShort"].isupper():
            # usually something like p.S234f , where f is actually part of "fs" (frameshift)
            return None
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    if "origAaLong" in groups:
        origSeq = threeToOneLower[groups["origAaLong"].lower()]
        correctProba = 0.9

    if "mutAaShort" in groups:
        mutSeq = groups["mutAaShort"]
    if "mutAaLong" in groups:
        mutSeq = threeToOneLower[groups["mutAaLong"].lower()]

    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.7
    if "mutDna" in groups:
        mutSeq = groups["mutDna"]

    if "origDna" in groups and "mutDna" in groups:
        if groups["origDna"].isupper() and groups["mutDna"].islower() or \
           groups["origDna"].islower() and groups["mutDna"].isupper():
            return None

    mutSeq = mutSeq.upper()
    origSeq = origSeq.upper()

    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos

    if "toPos" in groups:
        seqEnd = int(groups["toPos"])
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1

    if "length" in groups:
        length = int(groups["length"])
    else:
        length = len(origSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif origSeq and (("*" in origSeq.lower()) or ("x" in origSeq.lower())):
        variantType = "stoploss"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif variantType is None and seqType == "prot":
        if origSeq == mutSeq:
            variantType = "synonymous"
        else:
            variantType = "missense"

    if variantType == "frameshift":
        var = VariantDescription(mutType="fs", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                                 mutSeq=mutSeq.upper(), origStr=match.group(0).strip(), length=length, docId=docId,
                                 variantType=variantType, correctProba=correctProba)
    else:
        var = VariantDescription(mutType="sub", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                                 mutSeq=mutSeq.upper(), origStr=match.group(0).strip(), length=length, docId=docId,
                                 variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDel(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    logger.debug("Parsing match del %s into: %s" % (match.group(0), str(groups)))

    correctProba = 0.4
    haveStart = False
    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos
        haveStart = True
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"].upper()
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()].upper()

    haveEnd = False
    if "toPos" in groups:
        seqEnd = int(groups["toPos"]) + 1
        haveEnd = True
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"].upper()
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()].upper()
        correctProba = 0.95

    # these are heuristics ... we might still deal with a stoploss if origSeq is not given, for example ...
    origSeq = None
    if "origAasShort" in groups:
        origSeq = groups["origAasShort"]
        correctProba = 0.9
    if "origAasLong" in groups:
        origSeq = ''.join([threeToOneLower[groups["origAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["origAasLong"]), 3)])
        correctProba = 0.9
    if "origDnas" in groups:
        origSeq = groups["origDnas"]
        correctProba = 0.7
    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.7
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    elif "origAaLong" in groups:
        origSeq = threeToOneLower[groups["origAaLong"].lower()]
        correctProba = 0.9
    if origSeq:
        origSeq = origSeq.upper()

    if "length" in groups and groups["length"]:
        length = int(groups["length"])
        if haveStart and haveEnd and (length != seqEnd - seqStart):
            logger.debug("Warning: length %d != seqStart %d - seqEnd %d; setting to seqEnd-seqStart (we might confuse reference numbers for length)" % (length, seqEnd, seqStart))
            length = seqEnd - seqStart  # toPos + 1 already up there ...
    elif origSeq:
        length = len(origSeq)
    else:
        length = seqEnd - seqStart

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif firstAa and (("*" in firstAa.lower()) or ("x" in firstAa.lower())):
        variantType = "stoploss"
    elif secondAa and (("*" in secondAa.lower()) or ("x" in secondAa.lower())):
        variantType = "stoploss"
    elif origSeq and (("*" in origSeq.lower()) or ("x" in origSeq.lower())):
        variantType = "stoploss"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if length % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="del", seqType=seqType, start=seqStart, end=seqEnd, 
                             origSeq=origSeq, mutSeq=None, origStr=match.group(0).strip(), firstAa=firstAa,
                             secondAa=secondAa, length=length, docId=docId, variantType=variantType,
                             correctProba=correctProba)
    return var


def parseMatchIns(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()

    correctProba = 0.4
    if "fromPos" in groups and groups["fromPos"]:
        pos = int(groups["fromPos"])
        seqStart = pos
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"]
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()]
        correctProba = 0.9

    if "toPos" in groups and groups["toPos"]:
        seqEnd = int(groups["toPos"]) + 1
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 2
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"]
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()]
        correctProba = 0.95

    mutSeq = None
    if "mutAasShort" in groups:
        mutSeq = groups["mutAasShort"]
    if "mutAasLong" in groups:
        mutSeq = ''.join([threeToOneLower[groups["mutAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["mutAasLong"]), 3)])
    if "dnas" in groups:
        mutSeq = groups["dnas"]
    logger.debug("match: %s" % match.group(0))
    logger.info(str([g for g in groups]))
    if mutSeq:
        mutSeq = mutSeq.upper()

    if "length" in groups:
        length = int(groups["length"])
    else:
        length = len(mutSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if length % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="ins", seqType=seqType, start=seqStart, end=seqEnd,
                             origSeq=None, mutSeq=mutSeq, origStr=match.group(0).strip(),
                             firstAa=firstAa, secondAa=secondAa, length=length, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDelIns(match, patName, seqType, docId):
    groups = match.groupdict()

    correctProba = 0.4
    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"]
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()]
        correctProba = 0.9

    if "toPos" in groups:
        seqEnd = int(groups["toPos"]) + 1
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"]
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()]
        correctProba = 0.95

    mutSeq = None
    origSeq = None
    if "mutAasShort" in groups:
        mutSeq = groups["mutAasShort"]
    if "mutAasLong" in groups:
        mutSeq = ''.join([threeToOneLower[groups["mutAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["mutAasLong"]), 3)])
    if "dnas" in groups:
        mutSeq = groups["dnas"]
    if "origDnas" in groups:
        origSeq = groups["origDnas"]
    if "origDna" in groups:
        origSeq = groups["origDna"]
    if "mutDnas" in groups:
        mutSeq = groups["mutDnas"]
    logger.debug("match: %s" % match.group(0))
    logger.info(str([g for g in groups]))
    if mutSeq:
        mutSeq = mutSeq.upper()

    if "length" in groups:
        length = int(groups["length"])
    else:
        if origSeq:
            length = abs(len(origSeq) - len(mutSeq))
        else:
            length = len(mutSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if (seqEnd - seqStart + len(mutSeq)) % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="delins", seqType=seqType, start=seqStart, end=seqEnd,
                             origSeq=origSeq, mutSeq=mutSeq, origStr=match.group(0).strip(),
                             firstAa=firstAa, secondAa=secondAa, length=length, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDup(match, patName, seqType, docId):
    groups = match.groupdict()

    correctProba = 0.4
    origSeq = None
    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.9
    elif "origDnas" in groups:
        origSeq = groups["origDnas"]
        correctProba = 0.9

    if "pos" in groups:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = seqStart + 1
        length = 1
    elif "fromPos" in groups:
        assert "toPos" in groups
        seqStart = int(groups["fromPos"])
        seqEnd = int(groups["toPos"]) + 1
        length = seqEnd - seqStart
    else:
        assert False

    if origSeq:
        origSeq = origSeq.upper()
        # duplication ... in the old Indian fashion:
        mutSeq = origSeq + origSeq
        # ! duplication ... in the Assyro-Babylonian fashion:
        # ! mutSeq = origSeq * 2
    else:
        mutSeq = None

    if origSeq:
        if len(origSeq) != length:
            # assume that we picked up a reference accidentally ...
            logger.debug(("Warning: length %d != length of origSeq %d; setting to seqEnd-seqStart " + \
                         "(we might confuse reference numbers for length)") % (length, len(origSeq)))
            length = len(origSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif seqType == "dna" and length % 3 == 0:
        variantType = "nonframeshift"
    elif seqType == "dna" and length % 3 != 0:
        variantType = "frameshift"

    var = VariantDescription(mutType="dup", seqType=seqType, start=seqStart,
                             end=seqEnd, origSeq=origSeq, mutSeq=mutSeq,
                             origStr=match.group(0).strip(), length=length,
                             docId=docId, variantType=variantType, correctProba=correctProba)
    return var


def parsePlusMinus(plusMinus):
    if plusMinus == "+" or plusMinus == u"\u00FE" or plusMinus == u"\u2AF9":
        return 1
    elif plusMinus == "-" or plusMinus == "^" or plusMinus == u"\u2AFA":
        return -1
    else:
        assert False, str(plusMinus)


def parseMatchSplicing(match, patName, seqType, docId):
    # dna splicing        {sep}c\.{pos}{plusMinus}{offset}{origDna}>{mutDna}
    groups = match.groupdict()
    seqStart = int(groups["pos"])
    seqEnd = seqEnd = seqStart + 1
    plusMinus = parsePlusMinus(groups["plusMinus"])
    offset = int(groups["offset"])
    offset *= plusMinus
    logger.info("Match %s orig offset %d, plusMinus: %s" % (match.group(0), offset, plusMinus))
    logger.info("Match %s signed offset %d" % (match.group(0), offset))
    origSeq = groups["origDna"]
    mutSeq = groups["mutDna"]
    var = VariantDescription(mutType="splicing", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                             mutSeq=mutSeq.upper(), origStr=match.group(0).strip(), offset=offset, docId=docId,
                             variantType="splicing", correctProba=0.5)
    return var


def convertIVSNumber(ivsNumber):
    if re.match(r'[IV]+', ivsNumber):
        if ivsNumber == "I":
            return 1
        elif ivsNumber == "II":
            return 2
        elif ivsNumber == "III":
            return 3
        elif ivsNumber == "IV":
            return 4
        elif ivsNumber == "V":
            return 5
        elif ivsNumber == "VI":
            return 6
        elif ivsNumber == "VII":
            return 7
        elif ivsNumber == "VIII":
            return 8
        elif ivsNumber == "IX":
            return 9
        else:
            # what is this BS ... do they think people want exercise in translating roman numerals
            logger.debug("too lazy to convert IVS roman numeral %s" % (ivsNumber))
            return None
    else:
        return int(ivsNumber)
    assert False


def parseMatchIVS(match, patName, seqType, docId):
    groups = match.groupdict()
    ivsNumber = convertIVSNumber(groups["ivsNumber"])
    if not ivsNumber:
        return None
    plusMinus = parsePlusMinus(groups["plusMinus"])
    offset = int(groups["offset"])
    offset *= plusMinus
    logger.info("Match %s orig offset %d, plusMinus: %s" % (match.group(0), offset, plusMinus))
    logger.info("Match %s signed offset %d" % (match.group(0), offset))
    origSeq = groups["origDna"]
    mutSeq = groups["mutDna"]
    var = VariantDescription(mutType="ivssub", seqType=seqType, start=None, end=None, origSeq=origSeq.upper(),
                             mutSeq=mutSeq.upper(), origStr=match.group(0).strip(), offset=offset, docId=docId,
                             variantType="splicing", ivsNumber=ivsNumber, correctProba=0.3)
    return var


def parseMatchFS(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    # grab long and short versions of amino acid
    correctProba = 0.4
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    if "origAaLong" in groups:
        correctProba = 0.9
        origSeq = threeToOneLower[groups["origAaLong"].lower()]

    origSeq = origSeq.upper()
    pos = int(groups["pos"])
    # if isBlacklisted(origSeq, pos, mutSeq):
    #     return None
    seqStart = pos
    seqEnd = pos + 1

    variantType = "frameshift"

    var = VariantDescription(mutType="fs", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                             mutSeq=None, origStr=match.group(0).strip(), length=None, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def findVariantDescriptions(text, docId, exclPos=set()):
    """ put mutation mentions from document together into dicts indexed by normal form 
        return dict of "prot"|"dna"|"dbSnp" -> list of (VariantDescription, list of Mention)
        uses global variable "regexes", see loadDb()

    >>> findVariantDescriptions("The R71G BRCA1 mutation is really a p.R71G mutation")
    {'prot': [(VariantDescription(mutType=u'sub',seqType=u'prot',seqId=u'None',geneId=u'',start=u'70',end=u'71',origSeq=u'R',mutSeq=u'G'), [Mention(patName=u'{sep}p\\\\.\\\\(?{origAaShort}{pos}{mutAaShort}{fs}', start=35, end=42), Mention(patName=u'{sep}{origAaShort}{pos}{mutAaShort}', start=3, end=8)])]}
    """
    regexes = parseRegex('data/regex.txt')

    exclPos = set(exclPos)
    varMentions = defaultdict(list)
    varDescObj = {}
    for seqType, mutType, patName, pat in regexes:
        for match in pat.finditer(text):
            logger.debug("Match: Pattern %s, text %s" % (patName, match.groups()))
            if isOverlapping(match, exclPos):
                logger.debug("Overlapping with exclPos")
                continue
            if mutType == "sub":
                variant = parseMatchSub(match, patName, seqType, docId=docId)
            elif mutType == "del":
                variant = parseMatchDel(match, patName, seqType, docId=docId)
            elif mutType == "ins":
                variant = parseMatchIns(match, patName, seqType, docId=docId)
            elif mutType == "dup":
                variant = parseMatchDup(match, patName, seqType, docId=docId)
            elif mutType == "splicing":
                variant = parseMatchSplicing(match, patName, seqType, docId=docId)
            elif mutType == "ivssub":  # a different splicing notation
                variant = parseMatchIVS(match, patName, seqType, docId=docId)
            elif mutType == "delins":
                variant = parseMatchDelIns(match, patName, seqType, docId=docId)
            elif mutType == "fs":
                variant = parseMatchFS(match, patName, seqType, docId=docId)
            else:
                logger.debug("Ignoring match %s; don't know how to handle" % str(match.group(0).encode('utf-8').strip()))
                continue
            if variant == None:
                continue

            mention = makeMention(match, patName)
            varDescObj[variant.getName()] = variant
            varMentions[variant.getName()].append(mention)
            debugSnip = pubAlg.getSnippet(text, mention.start, mention.end, maxContext=60)
            logger.debug("Found Variant: %s, snippet %s" % (str(variant), debugSnip))

    # convert to dict of "prot"|"dna"|"dbSnp" -> list (variant, mentions)
    variants = {}
    variants["prot"] = []
    variants["dna"] = []
    variants["dbSnp"] = []

    for varName, mentions in varMentions.items():
        variant = varDescObj[varName]
        variants[variant.seqType].append((variant, mentions))
    variants = dict(variants)
    return variants


if __name__ == '__main__':
    with open(sys.argv[1], 'r', encoding= 'utf-8') as f:
        text = f.read().strip()

    typeToDescription = findVariantDescriptions(text, docId="")
    with open('type_to_description.txt', 'w+') as f:
        for entry in typeToDescription:
            f.write(str(entry) + "\n")
            f.write(str(typeToDescription.get(entry, [])) + "\n")
    dbSnp = typeToDescription.get("dbSnp", [])

    print('prot', '=' * 120)
    for prot in typeToDescription['prot']: 
        print(prot)
    
    print('dna', '=' * 120)
    for dna in typeToDescription['dna']: 
        print(dna)
    
    print('dbSnp', '=' * 120)
    for dbSnp in typeToDescription['dbSnp']: 
        print(dbSnp)
