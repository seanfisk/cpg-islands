import pytest
from mock import create_autospec, MagicMock, call

from cpg_islands.models import SeqInputModel, MetaResultsModel
from tests.helpers import fixture_file, read_fixture_file

@pytest.fixture
def model():
    mock_results_model = create_autospec(MetaResultsModel, spec_set=True)
    return SeqInputModel(mock_results_model)

class TestSeqInputModel:
    # fix all of these, they don't raise errors anymore
    class TestLoadFile:
        def test_load_jx500709_1(self, model):
            file_loaded_callback = MagicMock()
            model.file_loaded.append(file_loaded_callback)
            computed = model.load_file(fixture_file('JX500709.1.gb'))
            assert (file_loaded_callback.mock_calls ==
                    [call(read_fixture_file('JX500709.1.flattened'))])

        def test_load_genbank_no_dna(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('JX500709.1.no-dna.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('Premature end of line during sequence data')])

        def test_load_genbank_empty(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('empty.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('No records found in handle')])

        def test_load_genbank_two_records(self, model):
            error_raised_callback = MagicMock()
            model.error_raised.append(error_raised_callback)
            model.load_file(fixture_file('U49845.1-and-JX500709.1.gb'))
            assert (error_raised_callback.mock_calls ==
                    [call('More than one record found in handle')])
