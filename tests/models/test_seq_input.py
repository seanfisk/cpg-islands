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
        assert callback.mock_calls == [call(200, 0.5, 0.6)]
        assert model.results_model.mock_calls == []

    @patch('cpg_islands.models.algorithms', autospec=True, spec_set=True)
    def test_load_algorithms(self, mock_algorithms, model):
        registry = []
        for i in xrange(5):
            algo = MagicMock()
            algo.name = 'fakealgo{0}'.format(i)
            registry.append(algo)
        mock_algorithms.registry = registry
        callback = MagicMock()
        model.algorithms_loaded.append(callback)
        model.load_algorithms()
        assert (callback.mock_calls ==
                [call(['fakealgo{0}'.format(i) for i in xrange(5)])])

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

    @patch('cpg_islands.models.algorithms', autospec=True, spec_set=True)
    class TestComputeIslands:
        @pytest.mark.parametrize('algo_index', range(5))
        def test_correct_algorithm_called(self, mock_algorithms,
                                          model, algo_index):
            unselected_algo_indices = set(range(5)) - set([algo_index])
            registry = [MagicMock() for _ in xrange(5)]
            for i in unselected_algo_indices:
                registry[i].algorithm.return_value = sentinel.wrong_results
            registry[algo_index].algorithm.return_value = \
                sentinel.correct_results
            mock_algorithms.registry = registry

            model.compute_islands(sentinel.seq_record,
                                  sentinel.island_size,
                                  sentinel.min_gc_ratio,
                                  sentinel.min_obs_exp_cpg_ratio,
                                  algo_index)
            for i in unselected_algo_indices:
                assert mock_algorithms.registry[i].mock_calls == []
            assert (mock_algorithms.registry[algo_index].mock_calls ==
                    [call.algorithm(sentinel.seq_record,
                                    sentinel.island_size,
                                    sentinel.min_gc_ratio,
                                    sentinel.min_obs_exp_cpg_ratio)])

        def test_results_set(self, mock_algorithms, model):
            # Mock out algorithm return value.
            first_algo = MagicMock()
            first_algo.algorithm.return_value = sentinel.results
            first_algo.name = sentinel.algo_name
            mock_algorithms.registry = [first_algo]

            # Mock out timing.
            class DefaultTimer(object):
                RETVALS = [24, 35]

                def __init__(self):
                    self.times_called = 0

                def __call__(self):
                    retval = self.RETVALS[self.times_called]
                    self.times_called += 1
                    return retval

            mock_default_timer = MagicMock()
            mock_default_timer.side_effect = DefaultTimer()

            with patch('timeit.default_timer', mock_default_timer):
                model.compute_islands(sentinel.seq_record,
                                      sentinel.island_size,
                                      sentinel.min_gc_ratio,
                                      sentinel.min_obs_exp_cpg_ratio,
                                      0)
            assert mock_default_timer.mock_calls == [call() for _ in xrange(2)]
            assert (model.results_model.mock_calls ==
                    [call.set_results(
                        sentinel.results, sentinel.algo_name, 11)])

        def test_islands_computed_called(self, mock_algorithms, model):
            callback = MagicMock()
            model.islands_computed.append(callback)
            model.compute_islands(sentinel.fake, sentinel.fake,
                                  sentinel.fake, sentinel.fake, sentinel.fake)
            assert callback.mock_calls == [call()]
