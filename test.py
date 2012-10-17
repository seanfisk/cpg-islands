#!/usr/bin/env python

# Test runner and style checker

from __future__ import print_function
import abc
import sys
import argparse
from collections import OrderedDict

from py.io import TerminalWriter
import pytest
import pep8

CODE_DIRECTORY = 'cpg_islands'
TESTS_DIRECTORY = 'tests'
CHECK_FILES = [CODE_DIRECTORY,
               TESTS_DIRECTORY,
               'setup.py',
               'test.py']


class TestRunner(object):
    """Abstract test runner base class."""
    __metaclass__ = abc.ABCMeta

    def __init__(self, terminal_writer):
        self.terminal_writer = terminal_writer

    @abc.abstractproperty
    def name(self):
        """Return the lowercase name of the test runner.

        :return: the name
        :rtype: :class:`str`
        """
        return ''

    @abc.abstractproperty
    def title(self):
        """Return the proper title of the test runner.

        :return: the title
        :type: :class:`str`
        """
        return ''

    @abc.abstractmethod
    def run(self):
        """Run the test runner.

        :return: number of errors, or exit code (0 is success, >1 is failure)
        :rtype: :class:`int`
        """
        self.terminal_writer.sep('=', self.title())


class UnitTestRunner(TestRunner):
    """Runner for pytest unit tests."""
    def name(self):
        return 'tests'

    def title(self):
        return 'Pytest Unit Tests'

    def run(self):
        """Run all unit tests.

        :return: whether tests were successful
        :rtype: :class:`int`
        """
        super(UnitTestRunner, self).run()
        test_args = []
        # run on multiple processors if possible
        try:
            import xdist
            import multiprocessing
            num = max(multiprocessing.cpu_count() / 2, 1)
            test_args += ['-n', str(num)]
        except ImportError:
            pass  # oh well
        test_args += ['--verbose', TESTS_DIRECTORY]
        return pytest.main(test_args)


class StyleGuideRunner(TestRunner):
    def __init__(self, terminal_writer):
        self.terminal_writer = terminal_writer

    def name(self):
        return 'pep8'

    def title(self):
        return 'PEP8 Style Guide'

    def run(self):
        """Run PEP8 style guide checker on code and test files.

        :return: the number of errors
        :rtype: :class:`int`
        """
        super(StyleGuideRunner, self).run()
        pep8_style = pep8.StyleGuide()
        report = pep8_style.check_files(CHECK_FILES)
        return report.total_errors


def main(argv=None):
    # for printing separators; idea/code stolen from pytest
    tw = TerminalWriter()
    test_runners = [runner(tw) for runner in TestRunner.__subclasses__()]
    success = True

    if argv is None:
        argv = sys.argv

    # parse arguments
    if len(argv) == 1:
        for runner in test_runners:
            success &= runner.run() == 0
    else:
        test_runners_dict = OrderedDict([(runner.name, runner) for
                                         runner in test_runners])
        check_help = 'What check to run ({0}).'.format(
            ' | '.join(test_runners_dict.iterkeys()))
        arg_parser = argparse.ArgumentParser(
            description='Test runner and style checker')
        arg_parser.add_argument('check', metavar='CHECK',
                                choices=test_runners_dict.iterkeys(),
                                help=check_help)
        args = arg_parser.parse_args(args=argv[1:])
        success = test_runners_dict[args.check].run() == 0

    return int(not success)

if __name__ == '__main__':
    sys.exit(main())
