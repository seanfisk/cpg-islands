import mock

from cpg_islands.models import ApplicationModel


def pytest_funcarg__model(request):
    return ApplicationModel()


def pytest_generate_tests(metafunc):
    if 'helparg' in metafunc.funcargnames:
        metafunc.parametrize('helparg', ['-h', '--help'])


class TestModels:
    class TestRun:
        # capfd argument allows capture of stdout/stderr based on file
        # descriptors
        # see <http://pytest.org/latest/capture.html> for more information
        def test_noargs(self, model, capfd):
            model.run(['progname'])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert  'CpG Island Finder' in out
            assert 'Author:' in out
            assert 'URL:' in out

        def test_help(self, model, helparg, capfd):
            with mock.patch('sys.exit') as mock_exit:
                model.run(['progname', helparg])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert 'usage' in out
            assert  'CpG Island Finder' in out
            assert 'Author:' in out
            assert 'URL:' in out
            mock_exit.assert_called_once_with(0)
