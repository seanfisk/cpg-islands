from mock import patch, call, sentinel, MagicMock

from cpg_islands.models import (MetaAppModel,
                                MetaSeqInputModel,
                                MetaResultsModel)
from cpg_islands.views import (BaseAppView,
                               BaseSeqInputView,
                               BaseResultsView)

# We don't want to require the `PySide' module for testing, since we
# are not testing our user interface. Just patch the entire module.
@patch.dict('sys.modules', {'PySide': MagicMock()})
class TestComposers:
    # Keep in mind that the order of mock passed as arguments starts
    # from the bottom up.
    @patch('cpg_islands.qt.composers.ResultsPresenter',
           autospec=True, spec_set=True)
    @patch('cpg_islands.qt.composers.ResultsView',
           autospec=BaseResultsView, spec_set=True)
    @patch('cpg_islands.qt.composers.ResultsModel',
           autospec=MetaResultsModel, spec_set=True)
    @patch('cpg_islands.qt.composers.SeqInputPresenter',
           autospec=True, spec_set=True)
    @patch('cpg_islands.qt.composers.SeqInputView',
           autospec=BaseSeqInputView, spec_set=True)
    @patch('cpg_islands.qt.composers.SeqInputModel',
           autospec=MetaSeqInputModel, spec_set=True)
    @patch('cpg_islands.qt.composers.AppPresenter',
           autospec=True, spec_set=True)
    @patch('cpg_islands.qt.composers.AppView',
           autospec=BaseAppView, spec_set=True)
    @patch('cpg_islands.qt.composers.AppModel',
           autospec=MetaAppModel, spec_set=True)
    def test_create_qt_presenter(
            self,
            mock_app_model, mock_app_view, mock_app_pres,
            mock_seq_input_model, mock_seq_input_view, mock_seq_input_pres,
            mock_results_model, mock_results_view, mock_results_pres):
        mock_results_model.return_value = sentinel.results_model
        mock_results_view.return_value = sentinel.results_view
        mock_seq_input_model.return_value = sentinel.seq_input_model
        mock_seq_input_view.return_value = sentinel.seq_input_view
        app_model = mock_app_model.return_value
        mock_app_view.return_value = sentinel.app_view

        app_pres = mock_app_pres.return_value
        # seq_input_pres = mock_seq_input_pres.return_value
        # results_pres = mock_results_pres.return_value

        from cpg_islands.qt.composers import create_app_presenter
        retval = create_app_presenter(sentinel.argv)
        assert retval == app_pres

        assert (mock_app_model.mock_calls ==
                call(sentinel.seq_input_model).run(sentinel.argv).call_list())
        assert (mock_app_view.mock_calls ==
                [call(sentinel.seq_input_view, sentinel.results_view)])
        assert (mock_app_pres.mock_calls ==
                call(app_model,
                     sentinel.app_view).register_for_events().call_list())
        assert (mock_seq_input_model.mock_calls ==
                [call(sentinel.results_model)])
        assert (mock_seq_input_view.mock_calls == [call()])
        assert (mock_seq_input_pres.mock_calls ==
                call(sentinel.seq_input_model,
                     sentinel.seq_input_view
                 ).register_for_events().call_list())
        assert (mock_results_model.mock_calls == [call()])
        assert (mock_results_view.mock_calls == [call()])
        assert (mock_results_pres.mock_calls ==
                call(sentinel.results_model,
                     sentinel.results_view).register_for_events().call_list())
