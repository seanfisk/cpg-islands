#!/usr/bin/env python

# Test runner and style checker

from __future__ import print_function
import abc
import sys
from collections import OrderedDict

from shovel import task
from py.io import TerminalWriter
# have to import these under different names since we have functions
# with these names
import pytest
import pep8

sys.path.append('.')

CODE_DIRECTORY = 'cpg_islands'
TESTS_DIRECTORY = 'tests'
CODE_FILES = [CODE_DIRECTORY,
              TESTS_DIRECTORY,
              'setup.py',
              'shovel.py']


class MetaTestRunner(object):
    """Abstract test runner base class."""
    __metaclass__ = abc.ABCMeta

    def __init__(self, terminal_writer=TerminalWriter()):
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


class UnitTestRunner(MetaTestRunner):
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
        # DISABLED FOR NOW
        # try:
        #     import xdist
        #     import multiprocessing
        #     num = max(multiprocessing.cpu_count() / 2, 1)
        #     test_args += ['-n', str(num)]
        # except ImportError:
        #     pass  # oh well
        test_args += ['--verbose', TESTS_DIRECTORY]
        return pytest.main(test_args)


class StyleGuideRunner(MetaTestRunner):
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
        report = pep8_style.check_files(CODE_FILES)
        if report.total_errors == 0:
            print('No style errors')
        return report.total_errors


@task
def test():
    """Run all pytest unit tests."""
    sys.exit(UnitTestRunner().run())


@task
def style():
    """Perform a PEP8 style check on the code."""
    sys.exit(StyleGuideRunner().run())


@task
def test_all():
    """Perform a style check and run all unit tests."""
    success = True
    for runner in MetaTestRunner.__subclasses__():
        success &= runner().run() == 0
    sys.exit(int(not success))
