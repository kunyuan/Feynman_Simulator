#!/usr/bin/env python
import unittest
import weight
#from solver.test import *

if __name__ == '__main__':

    #testsuite = unittest.TestLoader().discover('.')
    #unittest.TextTestRunner(verbosity=1).run(testsuite)

    testsuite = unittest.TestLoader().discover('./solver')
    unittest.TextTestRunner(verbosity=3).run(testsuite)
    suite = unittest.TestLoader().loadTestsFromModule(weight)
    unittest.TextTestRunner(verbosity=3).run(suite)
    #unittest.main()
