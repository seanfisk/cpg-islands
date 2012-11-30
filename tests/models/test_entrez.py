from __future__ import division

import pytest
from mock import sentinel, call, create_autospec, patch, MagicMock

from cpg_islands.models import EntrezModel, MetaSeqInputModel


@pytest.fixture
def model():
    mock_seq_input_model = create_autospec(MetaSeqInputModel, spec_set=True)
    return EntrezModel(mock_seq_input_model)


class TestEntrezModel:
    def test_search(self, model):
        with patch('cpg_islands.models.Entrez') as mock_entrez:
            mock_entrez.esearch.return_value = sentinel.handle
            mock_entrez.read.return_value = {
                'IdList': sentinel.id_list,
                'QueryTranslation': sentinel.query_translation}
            results = model.search(sentinel.search)
        assert results == {
            'IdList': sentinel.id_list,
            'QueryTranslation': sentinel.query_translation}
        assert mock_entrez.mock_calls == [
            call.esearch(db='nucleotide', term=sentinel.search),
            call.read(sentinel.handle)]

    def test_suggest(self, model):
        with patch('cpg_islands.models.Entrez') as mock_entrez:
            mock_entrez.espell.return_value = sentinel.handle
            mock_entrez.read.return_value = sentinel.result
            results = model.suggest(sentinel.text)
        assert results == sentinel.result
        assert mock_entrez.mock_calls == [
            call.espell(db='pubmed', term=sentinel.text),
            call.read(sentinel.handle)]

    def test_get_seq(self, model):
        with patch('cpg_islands.models.Entrez') as mock_entrez:
            with patch('cpg_islands.models.SeqIO') as mock_seqio:
                handle = MagicMock()
                mock_entrez.efetch.return_value = handle
                mock_seqio.read.return_value = sentinel.record
                record = model.get_seq(sentinel.id)
        assert record == sentinel.record
        assert mock_entrez.mock_calls == call.efetch(
            db='nucleotide',
            id=sentinel.id,
            rettype='gb',
            retmode='text').close().call_list()
        assert mock_seqio.mock_calls == [call.read(handle, 'genbank')]

    class TestLoadSeq:
        def test_file_loaded_called(self, model):
            model.load_seq(sentinel.seq)
            assert model.seq_input_model.mock_calls == [
                call.file_loaded(sentinel.seq)]

        def test_seq_loaded_event_called(self, model):
            callback = MagicMock()
            model.seq_loaded.append(callback)
            model.load_seq(sentinel.seq)
            assert callback.mock_calls == [call()]
