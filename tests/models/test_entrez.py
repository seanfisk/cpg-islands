from __future__ import division

import pytest
from mock import sentinel, call, create_autospec, patch, Mock

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
                handle = Mock()
                handle.close = Mock(return_value=True)
                mock_entrez.efetch.return_value = handle
                mock_record = Mock()
                mock_record.seq = sentinel.seq
                mock_seqio.read.return_value = mock_record
                record = model.get_seq(sentinel.id)
        assert record == mock_record
        assert mock_entrez.mock_calls == [
            call.efetch(db='nucleotide', id=sentinel.id, rettype='gb',
                        retmode='text'),
            call.efetch().close()]
        assert mock_seqio.mock_calls == [call.read(handle, 'genbank')]
