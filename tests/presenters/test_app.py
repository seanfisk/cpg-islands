import pytest
from mock import create_autospec, call

from cpg_islands.models import MetaAppModel
from cpg_islands.views import BaseAppView
from cpg_islands.presenters import AppPresenter


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaAppModel, spec_set=True)
    mock_view = create_autospec(BaseAppView, spec_set=True)
    return AppPresenter(mock_model, mock_view)


class TestApplicationPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert (presenter.model.mock_calls ==
                [call.started.append(presenter.view.start)])
        assert (presenter.view.mock_calls ==
                [call.file_load_requested.append(
                    presenter.model.load_file)])
