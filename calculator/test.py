#!/usr/bin/env python
import unittest
import weight

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromModule(weight)
    unittest.TextTestRunner(verbosity=2).run(suite)
