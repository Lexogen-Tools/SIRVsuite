import unittest
import os

class testSIRVsuite(unittest.TestCase):
    def test_help(self):
        exit_status = os.system("SIRVsuite -h")
        assert exit_status == 0

    def test_not_help(self):
        exit_status = os.system("SIRVsuite")
        assert exit_status > 0