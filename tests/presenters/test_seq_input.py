from mock import create_autospec, sentinel, call
import pytest

from cpg_islands.models import MetaSeqInputModel
from cpg_islands.views import BaseSeqInputView
from cpg_islands.presenters import SeqInputPresenter


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaSeqInputModel, spec_set=True)
    mock_view = create_autospec(BaseSeqInputView, spec_set=True)
    return SeqInputPresenter(mock_model, mock_view)


class TestSeqInputPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert (presenter.model.mock_calls ==
                [call.island_definition_defaults_set.append(
                    presenter._island_definition_defaults_set),
                 call.file_loaded.append(presenter.view.set_seq),
                 call.error_raised.append(presenter.view.show_error),
                 call.algorithms_loaded.append(presenter.view.set_algorithms)])
        assert (presenter.view.mock_calls ==
                [call.submitted.append(presenter._user_submits)])

    def test_island_defintion_defaults_set(self, presenter):
        presenter._island_definition_defaults_set(343, 0.5, 0.65)
        assert (presenter.view.mock_calls ==
                [call.set_island_size('343'),
                 call.set_min_gc_ratio('0.5'),
                 call.set_min_obs_exp_cpg_ratio('0.65')])

    class TestUserSubmits:
        def test_valid_values(self, presenter):
            """When the user clicks submit with valid values, the island
            locations are shown."""
            seq_str = 'ATATGCGCATAT'
            presenter._user_submits(
                seq_str, '4', '0.5', '0.65', sentinel.algo_index)
            # we cannot use assert_called_once_with or mock_cals because
            # these two Seq's use object comparison, and therefore are not
            # "equal"
            assert presenter.model.compute_islands.call_count == 1
            # call_args[0] is ordered arguments, call_args[1] is
            # keyword arguments
            args = presenter.model.compute_islands.call_args[0]
            assert len(args) == 5
            assert str(args[0].seq) == seq_str
            assert args[1:] == (4, 0.5, 0.65, sentinel.algo_index)
            assert presenter.view.mock_calls == []

        def test_invalid_sequence(self, presenter):
            """When the user submits a sequence that does not contain
            valid bases, they are shown an error."""
            presenter._user_submits('ABCD', '3', '0.5', '0.65', 0)
            assert (presenter.view.mock_calls == [call.show_error(
                'Sequence letters not within alphabet:\n'
                '  Alphabet: GATC\n'
                '  Sequence: ABCD')])
            assert presenter.model.mock_calls == []

        def test_invalid_island_size_type(self, presenter):
            """When the user submits an invalid type of island size, they
            are shown an error."""
            presenter._user_submits(
                'ATATGCGC', 'invalid size', '0.5', '0.65', 0)
            assert (presenter.view.mock_calls ==
                    [call.show_error(
                        'Invalid integer for island size: invalid size')])
            assert presenter.model.mock_calls == []

        def test_invalid_gc_type(self, presenter):
            """When the user submits an invalid type for GC ratio, they
            are shown an error.
            """
            presenter._user_submits('ATATGCGC', '3', 'invalid gc', '0.65', 0)
            assert (presenter.view.mock_calls ==
                    [call.show_error('Invalid ratio for GC: invalid gc')])
            assert presenter.model.mock_calls == []

        def test_invalid_min_obs_exp_cpg_type(self, presenter):
            """When the user submits an invalid type for minimum
            observed/expected CpG ratio, they are shown an error.
            """
            presenter._user_submits(
                'ATATGCGC', '3', '0.5', 'invalid obsexpcpg', 0)
            assert (presenter.view.mock_calls ==
                    [call.show_error(
                        'Invalid ratio for minimum observed/expected '
                        'CpG ratio: invalid obsexpcpg')])
            assert presenter.model.mock_calls == []

        def test_lowercase_sequence(self, presenter):
            """When the user submits a sequence with lowercase letters,
            these should be gracefully handled."""
            seq_str = 'ATatgcGCAtaT'
            presenter._user_submits(
                seq_str, '4', '0.5', '0.65', sentinel.algo_index)
            assert presenter.model.compute_islands.call_count == 1
            args = presenter.model.compute_islands.call_args[0]
            assert len(args) == 5
            assert str(args[0].seq) == 'ATATGCGCATAT'
            assert args[1:] == (4, 0.5, 0.65, sentinel.algo_index)
            assert presenter.view.mock_calls == []

    class TestLoadFile:
        def test_valid_sequence(self, presenter):
            presenter.model.load_file.return_value = sentinel.file_contents
            presenter._file_loaded(sentinel.file_path)
            assert (presenter.model.mock_calls ==
                    [call.load_file(sentinel.file_path)])
            assert (presenter.view.mock_calls ==
                    [call.set_seq(sentinel.file_contents)])

        def test_invalid_sequence(self, presenter):
            presenter.model.load_file.side_effect = ValueError(
                'this is a fake message')
            presenter._file_loaded(sentinel.file_path)
            assert (presenter.model.mock_calls ==
                    [call.load_file(sentinel.file_path)])
            assert (presenter.view.mock_calls ==
                    [call.show_error(
                        'Sequence parsing error: this is a fake message')])
