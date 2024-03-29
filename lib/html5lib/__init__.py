""" 
HTML parsing library based on the WHATWG "HTML5"
specification. The parser is designed to be compatible with existing
HTML found in the wild and implements well-defined error recovery that
is largely compatible with modern desktop web browsers.

Example usage:

import html5lib
f = open("my_document.html")
tree = html5lib.parse(f) 
"""
from __future__ import absolute_import
__version__ = "0.95-dev"
from .html5parser import HTMLParser, parse, parseFragment
from .treebuilders import getTreeBuilder
from .treewalkers import getTreeWalker
from .serializer import serialize
