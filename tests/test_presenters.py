from mock import create_autospec, sentinel, call
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from pytest import raises

from cpg_islands.models import MetaApplicationModel
from cpg_islands.views import BaseApplicationView
from cpg_islands.presenters import ApplicationPresenter
from helpers import make_features


def pytest_funcarg__presenter(request):
    mock_model = create_autospec(MetaApplicationModel, spec_set=True)
    mock_view = create_autospec(BaseApplicationView, spec_set=True)
    presenter = ApplicationPresenter(mock_model, mock_view)
    return presenter


class TestPresenters:
    def test_register_for_events(self, presenter):
        """When the application starts, the presenter registers for
        events fired by the view and the model."""
        presenter.register_for_events()
        presenter.model.started.append.assert_called_once_with(
            presenter.view.start)
        presenter.view.submitted.append.assert_called_once_with(
            presenter._user_submits)
        presenter.view.file_loaded.append.assert_called_once_with(
            presenter._file_loaded)

    def test_user_submits_valid_values(self, presenter):
        """When the user clicks submit with valid values, the island
        locations are shown."""
        feature_tuples = [(0, 5), (1, 6), (3, 8)]
        presenter.model.annotate_cpg_islands.return_value = \
            make_features(feature_tuples)
        seq_str = 'ATATGCGCATAT'
        presenter._user_submits(seq_str, '4', '0.5')
        seq = Seq(seq_str, IUPAC.unambiguous_dna)
        # we cannot use assert_called_once_with because these two
        # Seq's use object comparison, and therefore are not "equal"
        #presenter.model.annotate_cpg_islands.\
        #    assert_called_once_with(seq, 4, 0.5)
        assert presenter.model.annotate_cpg_islands.call_count == 1
        annotate_args = presenter.model.annotate_cpg_islands.call_args[0]
        assert len(annotate_args) == 3
        assert str(annotate_args[0]) == str(seq)
        assert annotate_args[1] == 4
        assert annotate_args[2] == 0.5
        presenter.view.set_locations.assert_called_once_with(feature_tuples)

    def test_user_submits_invalid_sequence(self, presenter):
        """When the user submits a sequence that does not contain
        valid bases, they are shown an error."""
        presenter._user_submits('ABCD', '3', '0.5')
        presenter.view.show_error.assert_called_once_with(
            '''Sequence letters not within alphabet:
  Alphabet: GATC
  Sequence: ABCD''')

    def test_user_submits_invalid_island_size_type(self, presenter):
        """When the user submits an invalid type of island size, they
        are shown an error."""
        presenter._user_submits('ATATGCGC', 'invalid size', '0.5')
        presenter.view.show_error.assert_called_once_with(
            'Invalid integer for island size: invalid size')

    def test_user_submits_invalid_gc_type(self, presenter):
        """When the user submits an invalid type for GC ratio, they
        are shown an error.
        """
        presenter._user_submits('ATATGCGC', '3', 'invalid gc')
        presenter.view.show_error('Invalid ratio for GC: invalid gc')

    def test_user_submits_lowercase_seqence(self, presenter):
        """When the user submits a sequence with lowercase letters,
        these should be gracefully handled."""
        feature_tuples = [(0, 5), (1, 6), (3, 8)]
        presenter.model.annotate_cpg_islands.return_value = \
            make_features(feature_tuples)
        seq_str = 'ATatgcGCAtaT'
        presenter._user_submits(seq_str, '4', '0.5')
        seq = Seq('ATATGCGCATAT', IUPAC.unambiguous_dna)
        assert presenter.model.annotate_cpg_islands.call_count == 1
        annotate_args = presenter.model.annotate_cpg_islands.call_args[0]
        assert len(annotate_args) == 3
        assert str(annotate_args[0]) == str(seq)
        assert annotate_args[1] == 4
        assert annotate_args[2] == 0.5
        presenter.view.set_locations.assert_called_once_with(feature_tuples)

    def test_user_loads_file(self, presenter):
        presenter.model.load_file.return_value = sentinel.file_contents
        presenter._file_loaded(sentinel.file_path)
        presenter.model.load_file.assert_called_once_with(sentinel.file_path)
        presenter.view.set_sequence.assert_called_once_with(
            sentinel.file_contents)

    def test_user_loads_file_invalid_sequence(self, presenter):
        presenter.model.load_file.side_effect = ValueError(
            'this is a fake message')
        presenter._file_loaded(sentinel.file_path)
        assert (presenter.model.mock_calls ==
                [call.load_file(sentinel.file_path)])
        assert presenter.view.mock_calls == [call.show_error(
            'Sequence parsing error: this is a fake message')]
