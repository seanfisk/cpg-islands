from Bio.SeqFeature import SeqFeature, FeatureLocation


def make_features(tuples):
    """Build a list of features from a list of tuples.

    :param tuples: list of tuples
    :type tuples: :class:`list` of :class:`tuple`
    :return: list of features
    :rtype: :class:`list` of :class:`SeqFeature`
    """
    return [SeqFeature(FeatureLocation(a, b)) for a, b in tuples]
