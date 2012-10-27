from __future__ import division

import mock

from pytest import raises
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from cpg_islands.models import ApplicationModel
from helpers import make_features


def pytest_funcarg__model(request):
    return ApplicationModel()


def pytest_generate_tests(metafunc):
    if 'helparg' in metafunc.funcargnames:
        metafunc.parametrize('helparg', ['-h', '--help'])


def extract_features(feature_list, sequence):
    return [str(feature.extract(sequence)) for feature in feature_list]


def assert_features_equal(computed_features, expected_features, sequence):
    computed_features_extracted = extract_features(computed_features, sequence)
    expected_features_extracted = extract_features(expected_features, sequence)
    assert computed_features_extracted == expected_features_extracted


class TestModels:
    class TestRun:
        # capfd argument allows capture of stdout/stderr based on file
        # descriptors
        # see <http://pytest.org/latest/capture.html> for more information
        def test_noargs(self, model, capfd):
            model.run(['progname'])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert out == ''

        def test_help(self, model, helparg, capfd):
            with mock.patch('sys.exit') as mock_exit:
                model.run(['progname', helparg])
            out, err = capfd.readouterr()
            # some basic tests to check output
            assert 'usage' in out
            assert  'CpG Island Finder' in out
            assert 'Author:' in out
            assert 'URL:' in out
            mock_exit.assert_called_once_with(0)

        def test_started_called(self, model):
            started_callback = mock.MagicMock()
            model.started.append(started_callback)
            model.run(['progname'])
            started_callback.assert_called_once_with()

    class TestAnnotate:
        def test_zero_island_size(self, model):
            with raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           0, 0)
            # exc_info.value returns the actual exception
            assert str(exc_info.value) == 'Invalid island size: 0'

        def test_empty_sequence(self, model):
            with raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           1, 0)
                assert (str(exc_info.value) ==
                        'Island size (1) must be less than or '
                        'equal to sequence length (0)')

        def test_negative_island_size(self, model):
            with raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('', IUPAC.unambiguous_dna),
                                           -1, 0)
            assert str(exc_info.value) == 'Invalid island size: -1'

        def test_single_base(self, model):
            seq = Seq('C', IUPAC.unambiguous_dna)
            computed = model.annotate_cpg_islands(seq, 1, 1)
            expected = make_features([(0, 1)])
            assert_features_equal(computed, expected, seq)

        def test_ratio_one(self, model):
            seq = Seq('ATGCCGATTTTA', IUPAC.unambiguous_dna)
            computed = model.annotate_cpg_islands(seq, 4, 1)
            expected = make_features([(2, 6)])
            assert_features_equal(computed, expected, seq)

        def test_ratio_half(self, model):
            seq = Seq('ATATGCTAAT', IUPAC.unambiguous_dna)
            computed = model.annotate_cpg_islands(seq, 4, 0.5)
            expected = make_features([(2, 6),
                                      (3, 7),
                                      (4, 8)])
            assert_features_equal(computed, expected, seq)

        def test_ratio_third(self, model):
            seq = Seq('TTATATGCTAATAT', IUPAC.unambiguous_dna)
            size = 6
            computed = model.annotate_cpg_islands(seq, size, 1 / 3)
            expected = make_features([(i, i + size) for i in xrange(2, 7)])
            assert_features_equal(computed, expected, seq)

        def test_island_size_less_than_sequence_size(self, model):
            """When the user submits an island size greater than the
            sequence size, they are shown an error.
            """
            with raises(ValueError) as exc_info:
                model.annotate_cpg_islands(Seq('ATATGCGC',
                                               IUPAC.unambiguous_dna),
                                           9, 0.5)
            assert ((str(exc_info.value)) ==
                    'Island size (9) must be less than or '
                    'equal to sequence length (8)')
