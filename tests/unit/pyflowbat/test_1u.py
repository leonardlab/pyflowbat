import unittest
import pyflowbat as pfb


class Test1(unittest.TestCase):
    """Functional tests for sandbox CLI."""

    def test1(self):
        # ...and output value.
        self.assertEqual(
            1, 1, "first check failed"
        )
