import enCount.encode as encode
import unittest


class TestEncode(unittest.TestCase):

    def test_get_online_list(self):
        experiments = encode.get_online_list()
        self.assertTrue(len(experiments) > 0 )

