from __future__ import division

import pytest
from mock import sentinel, call, create_autospec

from cpg_islands.models import EntrezModel, MetaSeqInputModel


@pytest.fixture
def model():
    mock_seq_input_model = create_autospec(MetaSeqInputModel, spec_set=True)
    return EntrezModel(mock_seq_input_model)


class TestEntrezModel:
    def test_search(self, model):
        test_model = create_autospec(model, spec_set=True)
        test_model.search(sentinel.search)
        assert test_model.mock_calls == [call.search(sentinel.search)]

    def test_suggest(self, model):
        test_model = create_autospec(model, spec_set=True)
        test_model.suggest(sentinel.suggest)
        assert test_model.mock_calls == [call.suggest(sentinel.suggest)]

    def test_get_seq(self, model):
        test_model = create_autospec(model, spec_set=True)
        test_model.get_seq(sentinel.seq_record)
        assert test_model.mock_calls == [call.get_seq(sentinel.seq_record)]
