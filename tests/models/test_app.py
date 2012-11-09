import pytest
from mock import create_autospec, call, sentinel, MagicMock, patch

from cpg_islands import metadata
from cpg_islands.models import AppModel, MetaSeqInputModel, MetaResultsModel


@pytest.fixture
def model():
    mock_seq_input_model = create_autospec(MetaSeqInputModel, spec_st=True)
    mock_results_model = create_autospec(MetaResultsModel, spec_set=True)
    return AppModel(mock_seq_input_model, mock_results_model)


@pytest.fixture(params=['-h', '--help'])
def helparg(request):
    return request.param


@pytest.fixture(params=['-V', '--version'])
def versionarg(request):
    return request.param


class TestAppModel:
    def test_register_for_events(self, model):
        model.register_for_events()
        assert (model.results_model.mock_calls ==
                [call.locations_computed.append(model._locations_computed)])
        assert model.seq_input_model.mock_calls == []

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
            assert (model.seq_input_model.mock_calls ==
                    [call.set_island_definition_defaults()])
            assert model.results_model.mock_calls == []

        def test_help(self, model, helparg, capfd):
            with patch('sys.exit') as mock_exit:
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
            assert model.seq_input_model.mock_calls == []
            assert model.results_model.mock_calls == []

        def test_version(self, model, versionarg, capfd):
            with patch('sys.exit') as mock_exit:
                mock_exit.side_effect = Exception(
                    'fake exception to stop execution')
                with pytest.raises(Exception):
                    model.run(['progname', versionarg])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert err == '{0} {1}\n'.format(metadata.nice_title,
                                             metadata.version)
            assert mock_exit.mock_calls == [call(0)]
            assert model.seq_input_model.mock_calls == []
            assert model.results_model.mock_calls == []

        def test_started_and_defaults_set_called_after_exit(
                self, model, helparg):
            started_callback = MagicMock()
            model.started.append(started_callback)
            with patch('sys.exit') as mock_exit:
                mock_exit.side_effect = Exception(
                    'fake exception to stop execution')
                with pytest.raises(Exception):
                    model.run(['progname', helparg])
            assert model.seq_input_model.mock_calls == []
            assert started_callback.mock_calls == []
            assert model.results_model.mock_calls == []

        def test_started_and_defaults_set_called(self, model):
            started_callback = MagicMock()
            model.started.append(started_callback)
            model.run(['progname'])
            assert started_callback.mock_calls == [call()]
            assert (model.seq_input_model.mock_calls ==
                    [call.set_island_definition_defaults()])
            assert model.results_model.mock_calls == []

    def test_load_file(self, model):
        model.load_file(sentinel.file_path)
        assert (model.seq_input_model.mock_calls ==
                [call.load_file(sentinel.file_path)])
        assert model.results_model.mock_calls == []

    def test_locations_computed(self, model):
        callback = MagicMock()
        model.locations_computed.append(callback)
        model._locations_computed(['unused', 'list', 'of', 'features'])
        assert callback.mock_calls == [call()]
