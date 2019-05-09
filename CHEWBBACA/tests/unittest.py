#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:02:34 2019

@author: pcerqueira
"""


import unittest
from unittest.mock import patch
import sys

def thing_that_uses_argparse():
  return sys.argv

class TestThings(unittest.TestCase):
    def test_parse_args(self):
        testargs = ["one", "two"]
        with patch.object(sys, 'argv', testargs):
            things = thing_that_uses_argparse()
            print("Inside with: %s" % things)
            assert(things == ["one", "two"])
        print("Outside with: %s" % thing_that_uses_argparse())

if __name__ == '__main__':
    unittest.main()