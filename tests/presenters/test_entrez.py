from mock import create_autospec, call, sentinel, Mock
import pytest

from cpg_islands.models import MetaEntrezModel
from cpg_islands.views import BaseEntrezView
from cpg_islands.presenters import EntrezPresenter


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaEntrezModel, spec_set=False)
    mock_view = create_autospec(BaseEntrezView, spec_set=True)
    return EntrezPresenter(mock_model, mock_view)


class TestEntrezPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert presenter.view.mock_calls == [
            call.text_changed.append(presenter._text_changed),
            call.searched.append(presenter._user_submits),
            call.result_selected.append(presenter._user_selected),
            call.load.append(presenter._load_selected)]

    class TestUserSubmits:
        def test_valid_values(self, presenter):
            """When the user clicks search with a valid string,
                the search results are shown."""
            record = {'IdList': sentinel.id_list,
                      'QueryTranslation': sentinel.query_translation}
            presenter.model.search.return_value = record
            presenter._user_submits(sentinel.search_str)
            assert (presenter.model.mock_calls ==
                    [call.search(sentinel.search_str)])
            assert (presenter.view.mock_calls ==
                    [call.set_result(sentinel.id_list),
                     call.set_query(sentinel.query_translation)])

    class TestUserSelected:
        def test_selected(self, presenter):
            sentinel.seq.seq = 'GCGC'
            presenter.model.results = {'IdList': [sentinel.seq]}
            presenter.model.get_seq.return_value = sentinel.seq
            presenter._user_selected(0)
            assert (presenter.model.mock_calls ==
                    [call.get_seq(sentinel.seq)])
            assert (presenter.view.mock_calls ==
                    [call.set_seq(sentinel.seq.seq)])

    class TestTextChanged:
        def test_changed(self, presenter):
            sentinel.results = {'CorrectedQuery': sentinel.corrected}
            presenter.model.suggest.return_value = sentinel.results
            presenter._text_changed(sentinel.text)
            assert (presenter.model.mock_calls ==
                    [call.suggest(sentinel.text)])
            assert (presenter.view.mock_calls ==
                    [call.set_suggestion(sentinel.corrected)])

    class TestLoadSelected:
        def test_changed(self, presenter):
            presenter.model.seq_input_model = Mock()
            presenter.model.seq_input_model.file_loaded = Mock(
                return_value=None)
            presenter._load_selected(sentinel.seq)
            assert (presenter.model.mock_calls ==
                    [call.load_seq(sentinel.seq)])
