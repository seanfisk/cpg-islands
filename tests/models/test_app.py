import pytest
from mock import create_autospec, call, sentinel, MagicMock, patch

from cpg_islands import metadata
from cpg_islands.models import AppModel, MetaSeqInputModel, MetaEntrezModel


@pytest.fixture
def model():
    mock_seq_input_model = create_autospec(MetaSeqInputModel, spec_set=True)
    mock_entrez_model = create_autospec(MetaEntrezModel, spec_set=True)
    return AppModel(mock_seq_input_model, mock_entrez_model)


@pytest.fixture(params=['-h', '--help'])
def helparg(request):
    return request.param


@pytest.fixture(params=['-V', '--version'])
def versionarg(request):
    return request.param


class TestAppModel:
    def test_register_for_events(self, model):
        model.register_for_events()
        assert model.seq_input_model.mock_calls == [
            call.islands_computed.append(model.islands_computed)]
        assert model.entrez_model.mock_calls == [
            call.seq_loaded.append(model.seq_loaded)]

    class TestRun:
        # `capfd' argument allows capture of stdout/stderr based on
        # file descriptors.
        #
        # See <http://pytest.org/latest/capture.html> for more
        # information
        def test_noargs(self, model, capfd):
            model.run(['progname'])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert out == ''

        def test_help(self, model, helparg, capfd):
            with patch('sys.exit', autospec=True, spec_set=True) as mock_exit:
                mock_exit.side_effect = Exception(
                    'fake exception to stop execution')
                with pytest.raises(Exception):
                    model.run(['progname', helparg])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert 'usage' in out
            assert 'CpG Island Finder' in out
            assert 'Author:' in out
            assert 'URL:' in out
            assert mock_exit.mock_calls == [call(0)]

        def test_version(self, model, versionarg, capfd):
            with patch('sys.exit', autospec=True, spec_set=True) as mock_exit:
                mock_exit.side_effect = Exception(
                    'fake exception to stop execution')
                with pytest.raises(Exception):
                    model.run(['progname', versionarg])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert err == '{0} {1}\n'.format(metadata.nice_title,
                                             metadata.version)
            assert mock_exit.mock_calls == [call(0)]

        def test_started_and_defaults_set_not_called_after_exit(
                self, model, helparg):
            started_callback = MagicMock()
            model.started.append(started_callback)
            with patch('sys.exit', autospec=True, spec_set=True) as mock_exit:
                mock_exit.side_effect = Exception(
                    'fake exception to stop execution')
                with pytest.raises(Exception):
                    model.run(['progname', helparg])
            assert model.seq_input_model.mock_calls == []
            assert started_callback.mock_calls == []

        def test_started_defaults_set_and_algos_loaded_called_normally(
                self, model):
            started_callback = MagicMock()
            model.started.append(started_callback)
            model.run(['progname'])
            assert started_callback.mock_calls == [call()]
            assert (model.seq_input_model.mock_calls ==
                    [call.set_island_definition_defaults(),
                     call.load_algorithms()])

    def test_load_file(self, model):
        model.load_file(sentinel.file_path)
        assert (model.seq_input_model.mock_calls ==
                [call.load_file(sentinel.file_path)])
