import io
import unittest
from contextlib import redirect_stdout

from lqc.cli import main


class TestCli(unittest.TestCase):

    def test_version_flag(self):
        buf = io.StringIO()
        with self.assertRaises(SystemExit) as cm, redirect_stdout(buf):
            main(['--version'])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(buf.getvalue().strip(), 'lqc 0.0.5')


if __name__ == '__main__':
    unittest.main()
