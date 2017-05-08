#!/usr/bin/env python

# adapted from cogent package

import doctest, unittest
import sys, os
from glob import glob

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.

    __import__ only imports the top-level module.
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def suite():
    modules_to_test = [
            #common
            'test_util',
            'test_format',
            'test_mixture',
            'test_viterbi',
            #6 states
            'test_biot',
            'test_aiot',
            'test_detect_cnv',
            #9 states
            'test_emission',
            'test_transition',
            'test_hmm',
        ]
    alltests = unittest.TestSuite()
    for module in modules_to_test:
        if module.endswith('.rest'):
            module = os.path.join(*module.split(".")[:-1])+".rest"
            test = doctest.DocFileSuite(module, optionflags=
                doctest.REPORT_ONLY_FIRST_FAILURE |
                doctest.ELLIPSIS)
        else:
            test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        unittest.main(defaultTest='suite', argv=sys.argv)
