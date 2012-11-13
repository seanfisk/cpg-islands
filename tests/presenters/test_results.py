from mock import create_autospec, call
import pytest

from cpg_islands.models import MetaResultsModel
from cpg_islands.views import BaseResultsView
from cpg_islands.presenters import ResultsPresenter
from tests.helpers import make_features


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaResultsModel, spec_set=True)
    mock_view = create_autospec(BaseResultsView, spec_set=True)
    presenter = ResultsPresenter(mock_model, mock_view)
    return presenter


class TestResultPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert (presenter.view.mock_calls ==
                [call.global_highlight.append(
                    presenter._get_seq), call.feature_selected.append(
                        presenter._get_feature)])

    def test_locations_computed(self, presenter):
        feature_tuples = [(0, 5), (1, 6), (3, 8)]
        feature_locations = make_features(feature_tuples)
        presenter._locations_computed(feature_locations)
        assert (presenter.view.mock_calls ==
                [call.set_locations(feature_tuples)])

        assert presenter.model.mock_calls == []
