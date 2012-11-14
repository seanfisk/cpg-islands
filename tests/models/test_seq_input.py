import pytest
from mock import patch, create_autospec, MagicMock, call, sentinel

from cpg_islands.models import SeqInputModel, MetaResultsModel
from tests.helpers import fixture_file, read_fixture_file


@pytest.fixture
def model():
    mock_results_model = create_autospec(MetaResultsModel, spec_set=True)
    return SeqInputModel(mock_results_model)


class TestSeqInputModel:
    def test_set_island_definition_defaults(self, model):
        callback = MagicMock()
        model.island_definition_defaults_set.append(callback)
        model.set_island_definition_defaults()
        assert callback.mock_calls == [call(200, 0.6)]
        assert model.results_model.mock_calls == []

    class TestLoadFile:
        def test_load_jx500709_1(self, model):
            file_loaded_callback = MagicMock()
            model.file_loaded.append(file_loaded_callback)
            model.load_file(fixture_file('JX500709.1.gb'))
            assert (file_loaded_callback.mock_calls ==
                    [call(read_fixture_file('JX500709.1.flattened'))])
            assert model.results_model.mock_calls == []

        def test_load_genbank_no_dna(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('JX500709.1.no-dna.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('Premature end of line during sequence data')])
            assert model.results_model.mock_calls == []

        def test_load_genbank_empty(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('empty.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('No records found in handle')])
            assert model.results_model.mock_calls == []

        def test_load_genbank_two_records(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('U49845.1-and-JX500709.1.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('More than one record found in handle')])
            assert model.results_model.mock_calls == []

    @patch('cpg_islands.models.algorithms')
    class TestComputeIslands:
        def test_results_set(self, mock_algorithms, model):
            first_algo = mock_algorithms.registry[0].algorithm
            first_algo.return_value = sentinel.seq_record
            model.compute_islands(
                sentinel.seq, sentinel.island_size, sentinel.min_gc_ratio)
            assert (model.results_model.mock_calls ==
                    [call.set_results(sentinel.seq_record)])

        def test_islands_computed_called(self, mock_algorithms, model):
            callback = MagicMock()
            model.islands_computed.append(callback)
            model.compute_islands(
                sentinel.fake, sentinel.fake, sentinel.fake)
            assert callback.mock_calls == [call()]
