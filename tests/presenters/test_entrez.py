from mock import create_autospec, call, sentinel
import pytest

from cpg_islands.models import MetaEntrezModel
from cpg_islands.views import BaseEntrezView
from cpg_islands.presenters import EntrezPresenter
from tests.helpers import make_seq_record


@pytest.fixture
def presenter():
    mock_model = create_autospec(MetaEntrezModel, spec_set=True)
    mock_view = create_autospec(BaseEntrezView, spec_set=True)
    return EntrezPresenter(mock_model, mock_view)


class TestEntrezPresenter:
    def test_register_for_events(self, presenter):
        presenter.register_for_events()
        assert presenter.view.mock_calls == [
            call.query_changed.append(presenter._query_changed),
            call.search_requested.append(presenter._user_submits),
            call.result_selected.append(presenter._user_selected),
            call.load_requested.append(presenter._load_selected)]

    class TestUserSubmits:
        def test_valid_values(self, presenter):
            """When the user clicks search with a valid string,
                the search results are shown."""
            presenter.model.search.return_value = (sentinel.id_list,
                                                   sentinel.query_translation)
            presenter._user_submits(sentinel.search_str)
            assert (presenter.model.mock_calls ==
                    [call.search(sentinel.search_str)])
            assert (presenter.view.mock_calls ==
                    [call.set_result(sentinel.id_list),
                     call.set_query_translation(sentinel.query_translation)])

    def test_user_selected(self, presenter):
        seq_str = 'ATATACGCGCATATA'
        seq_id = 'NG_032827.2'
        seq_desc = "It's a pretty cool sequence"
        seq_record = make_seq_record(seq_str)
        seq_record.id = seq_id
        seq_record.description = seq_desc
        presenter.model.get_seq_record.return_value = seq_record
        presenter._user_selected(sentinel.index)
        assert presenter.model.mock_calls == [
            call.get_seq_record(sentinel.index)]
        assert presenter.view.mock_calls == [
            call.set_seq_locus(
                seq_id, 'http://www.ncbi.nlm.nih.gov/nuccore/NG_032827.2'),
            call.set_seq_desc(seq_desc),
            call.set_seq_len('15 bases'),
            call.set_selected_seq(seq_str)]

    def test_user_changed(self, presenter):
        presenter.model.suggest.return_value = sentinel.suggestion
        presenter._query_changed(sentinel.text)
        assert presenter.model.mock_calls == [call.suggest(sentinel.text)]
        assert presenter.view.mock_calls == [
            call.set_suggestion(sentinel.suggestion)]

    def test_load_selected(self, presenter):
        presenter._load_selected()
        assert presenter.model.mock_calls == [call.load_seq()]
