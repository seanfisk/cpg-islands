from mock import create_autospec, sentinel
from Bio.Seq import Seq

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
        presenter.register_for_events()
        presenter.model.started.append.assert_called_once_with(
            presenter.view.start)
        presenter.view.submitted.append.assert_called_once_with(
            presenter._user_submits)

    def test_user_submits_valid_values(self, presenter):
        """When the user clicks submit with valid values, the island
        locations are shown."""
        feature_tuples = [(0, 5), (1, 6), (3, 8)]
        presenter.model.annotate_cpg_islands.return_value = \
            make_features(feature_tuples)
        seq_str = 'atatgcgcatat'
        presenter._user_submits(seq_str, '4', '0.5')
        seq = Seq(seq_str)
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
