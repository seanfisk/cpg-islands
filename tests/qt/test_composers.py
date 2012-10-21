import mock

from cpg_islands.qt.composers import create_presenter


class TestComposers:
    @mock.patch('cpg_islands.qt.composers.ApplicationPresenter',
                autospec=True, spec_set=True)
    @mock.patch('cpg_islands.qt.composers.ApplicationView',
                autospec=True, spec_set=True)
    @mock.patch('cpg_islands.qt.composers.ApplicationModel',
                autospec=True, spec_set=True)
    def test_create_cli_presenter(self, mock_model, mock_view, mock_presenter):
        model = mock_model.return_value
        mock_view.return_value = mock.sentinel.view
        presenter = mock_presenter.return_value

        retval = create_presenter()
        assert retval == presenter

        mock_model.assert_called_once_with()
        mock_view.assert_called_once_with()
        expected_model_calls = mock.call().run().call_list()
        expected_presenter_calls = mock.call(model, mock.sentinel.view).\
            register_for_events().call_list()
        assert mock_model.mock_calls == expected_model_calls
        assert mock_presenter.mock_calls == expected_presenter_calls
