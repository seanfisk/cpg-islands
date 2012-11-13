import os.path

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


def make_features(tuples):
    """Build a list of features from a list of tuples.

    :param tuples: list of tuples
    :type tuples: :class:`list` of :class:`tuple`
    :return: list of features
    :rtype: :class:`list` of :class:`SeqFeature`
    """
    return [SeqFeature(FeatureLocation(a, b)) for a, b in tuples]


def make_seq_record(seq_str='', feature_tuples=[]):
    """Create a SeqRecord from a string.

    :param seq_str: the sequence as a string
    :type seq_str: :class:`str`
    :return: the built SeqRecord
    :rtype: :class:`SeqRecord`
    """
    return SeqRecord(seq=Seq(seq_str, IUPAC.unambiguous_dna),
                     features=make_features(feature_tuples))


def fixture_file(basename):
    return os.path.join('tests', 'fixtures', basename)


def read_fixture_file(basename):
    return open(fixture_file(basename)).read()
