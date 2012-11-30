import os.path

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from cpg_islands.algorithms import IslandMetadata, AlgoResults


def _make_feature(start, end):
    """Create a feature.

    :param start: start index
    :type start: :class:`int`
    :param end: end index
    :type end: :class:`int`
    :return: the created feature
    :rtype: :class:`SeqFeature`
    """
    return SeqFeature(FeatureLocation(start, end))


def make_seq_record(seq_str='', feature_tuples=[]):
    """Create a SeqRecord from a string.

    :param seq_str: the sequence as a string
    :type seq_str: :class:`str`
    :return: the built SeqRecord
    :rtype: :class:`SeqRecord`
    """
    return SeqRecord(
        seq=Seq(seq_str, IUPAC.unambiguous_dna),
        features=[_make_feature(start, end) for start, end in feature_tuples])


def make_algo_results(seq_str='', island_metadata_tuples=[]):
    island_features = []
    island_metadata_list = []
    for start, end, gc_ratio, obs_exp_cpg_ratio in island_metadata_tuples:
        island_features.append(_make_feature(start, end))
        island_metadata_list.append(
            IslandMetadata(gc_ratio, obs_exp_cpg_ratio))
    seq_record = SeqRecord(seq=Seq(seq_str, IUPAC.unambiguous_dna),
                           features=island_features)
    return AlgoResults(seq_record, island_metadata_list)


def fixture_file(basename):
    return os.path.join('tests', 'fixtures', basename)


def read_fixture_file(basename):
    return open(fixture_file(basename)).read()
