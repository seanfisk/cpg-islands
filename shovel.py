# Test runner and style checker

from __future__ import print_function
from abc import ABCMeta, abstractmethod, abstractproperty
import sys
import subprocess

from shovel import task
from py.io import TerminalWriter
import pytest
from pyfiglet import Figlet

sys.path.append('.')

from cpg_islands import metadata

CODE_DIRECTORY = 'cpg_islands'
TESTS_DIRECTORY = 'tests'
CODE_FILES = [CODE_DIRECTORY,
              TESTS_DIRECTORY,
              'setup.py',
              'shovel.py']


class MetaTestRunner(object):
    """Abstract test runner base class."""
    __metaclass__ = ABCMeta

    def __init__(self, terminal_writer=TerminalWriter()):
        self.terminal_writer = terminal_writer

    @abstractproperty
    def name(self):
        """Return the proper name of the test runner.

        :return: the name
        :type: :class:`str`
        """
        return ''

    @abstractmethod
    def run(self):
        """Run the test runner.

        :return: number of errors, or exit code (0 is success, >1 is failure)
        :rtype: :class:`int`
        """
        self.terminal_writer.sep('=', self.name)


class LintRunner(MetaTestRunner):
    @property
    def name(self):
        return 'Flake8: PyFlakes, PEP8, and McCabe complexity'

    def run(self):
        """Run flake8 on code and test files.

        :return: the number of errors
        :rtype: :class:`int`
        """
        super(LintRunner, self).run()
        # Flake8 doesn't have an easy way to run checks using a Python
        # function, so just fork off another process to do it.
        retcode = subprocess.call(['flake8', '--max-complexity=10'] +
                                  CODE_FILES)
        if retcode == 0:
            print('No style errors')
        return retcode


class UnitTestRunner(MetaTestRunner):
    """Runner for pytest unit tests."""
    @property
    def name(self):
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
        test_args.append(TESTS_DIRECTORY)
        return pytest.main(test_args)


def _test_all():
    """Abstraction function which returns a code instead of exiting."""
    success = True
    for runner in MetaTestRunner.__subclasses__():
        success &= runner().run() == 0
    return int(not success)


@task
def test():
    """Run all pytest unit tests."""
    sys.exit(UnitTestRunner().run())


@task
def lint():
    """Perform PEP8 style check, run PyFlakes, and run McCabe
    complexity metrics on the code.
    """
    sys.exit(LintRunner().run())


@task
def test_all():
    """Perform a style check and run all unit tests."""
    retcode = _test_all()
    if retcode == 0:
        text = 'PASSED'
    else:
        text = 'FAILED'
    print(Figlet(font='starwars').renderText(text))
    sys.exit(retcode)


@task
def travis():
    """Perform setup for Travis-CI and run the tests."""
    subprocess.check_call(['python', 'setup.py', 'build_ext', '--inplace'])
    test_all()


@task
def commit():
    """Commit only if all the tests pass."""
    if _test_all() == 0:
        subprocess.check_call(['git', 'commit'])
    else:
        print('\nTests failed, not committing.')


@task
def coverage():
    """Run tests and show test coverage report."""
    pytest.main(['--cov', CODE_DIRECTORY,
                 '--cov-report', 'term-missing',
                 TESTS_DIRECTORY])


@task
def qt():
    """Run the Qt-based version of the program."""
    from cpg_islands.qt.main import main
    main([])


@task
def search(query):
    """Search all project files for a certain string."""
    # Stolen directly from
    # http://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
    git_proc = subprocess.Popen(['git', 'ls-files'], stdout=subprocess.PIPE)
    ack_proc = subprocess.Popen(['xargs', 'ack', '--pager=less -R', query],
                                stdin=git_proc.stdout)
    # Allow git_proc to receive a SIGPIPE if ack_proc exits.
    git_proc.stdout.close()
    ack_proc.communicate()


@task
def mac_app():
    """Build an application bundle for Mac OS X."""
    subprocess.check_call(['python', 'setup.py', 'py2app'])


@task
def mac_dmg():
    """Build a disk image (installer) for Mac OS X."""
    mac_app()
    subprocess.check_call(['scripts/yoursway-create-dmg/create-dmg',
                           '--window-size', '500', '300',
                           '--volname', metadata.nice_title,
                           '--app-drop-link', '380', '205',
                           'dist/{0}.dmg'.format(metadata.nice_title),
                           'dist/{0}.app'.format(metadata.nice_title)])
