import mock

from cpg_islands.models import MetaApplicationModel
from cpg_islands.views import BaseApplicationView
from cpg_islands.presenters import ApplicationPresenter


def pytest_funcarg__presenter(request):
    mock_model = mock.create_autospec(MetaApplicationModel, spec_set=True)
    mock_view = mock.create_autospec(BaseApplicationView, spec_set=True)
    presenter = ApplicationPresenter(mock_model, mock_view)
    return presenter


class TestPresenters:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        presenter.model.started.append.assert_called_once_with(
            presenter.view.start)
