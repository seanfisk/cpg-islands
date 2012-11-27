from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet
from Bio.SeqRecord import SeqRecord


class AppPresenter(object):
    def __init__(self, model, view):
        """Constructor.

        :param model: application model
        :type model: :class:`MetaAppModel`
        :param view: application view
        :type view: :class:`MetaAppView`
        """
        self.model = model
        self.view = view

    def register_for_events(self):
        self.model.started.append(self.view.start)
        self.model.islands_computed.append(self.view.show_results)
        self.view.file_load_requested.append(self.model.load_file)


class SeqInputPresenter(object):
    def __init__(self, model, view):
        """Constructor.

        :param model: sequence input model
        :type model: :class:`MetaSeqInputModel`
        :param view: sequence input view
        :type view: :class:`MetaSeqInputView`
        """
        self.model = model
        self.view = view

    def register_for_events(self):
        self.model.island_definition_defaults_set.append(
            self._island_definition_defaults_set)
        self.model.file_loaded.append(self.view.set_seq)
        self.model.error_raised.append(self.view.show_error)
        self.model.algorithms_loaded.append(self.view.set_algorithms)
        self.view.submitted.append(self._user_submits)

    def _island_definition_defaults_set(
            self, island_size, min_gc_ratio, min_obs_exp_cpg_ratio):
        """Called when island definition defaults are set.

        :param island_size: number of bases in the island
        :type island_size: :class:`int`
        :param min_gc_ratio: minimum ratio of Guanine/Cytosine
        :type min_gc_ratio: :class:`float`
        """
        self.view.set_island_size(str(island_size))
        self.view.set_min_gc_ratio(str(min_gc_ratio))
        self.view.set_min_obs_exp_cpg_ratio(str(min_obs_exp_cpg_ratio))

    def _user_submits(
            self, seq_str, island_size_str, min_gc_ratio_str,
            min_obs_exp_cpg_ratio_str, algo_index):
        """Called when the user submits the form.

        :param seq_str: the sequence as a string
        :type seq_str: :class:`str`
        :param island_size_str: number of bases which an island may contain
        :type island_size_str: :class:`str`
        :param min_gc_ratio_str: the ratio of GC to other bases
        :type min_gc_ratio_str: :class:`str`
        :param algo_index: the algorithm chosen
        :type algo_index: :class:`int`
        """
        seq_mixed_case = Seq(seq_str, IUPAC.unambiguous_dna)
        seq = seq_mixed_case.upper()
        # Using `_verify_alphabet' is somewhat questionable, since it
        # is marked as private. However, there are no other documented
        # ways to verify the sequence.
        if not _verify_alphabet(seq):
            self.view.show_error(
                '''Sequence letters not within alphabet:
  Alphabet: {0}
  Sequence: {1}'''.format(seq.alphabet.letters, str(seq)))
            return
        try:
            island_size = int(island_size_str)
        except ValueError:
            self.view.show_error(
                'Invalid integer for island size: {0}'.format(island_size_str))
            return
        try:
            min_gc_ratio = float(min_gc_ratio_str)
        except ValueError:
            self.view.show_error(
                'Invalid ratio for GC: {0}'.format(min_gc_ratio_str))
            return
        try:
            min_obs_exp_cpg_ratio = float(min_obs_exp_cpg_ratio_str)
        except ValueError:
            self.view.show_error(
                'Invalid ratio for minimum observed/expected '
                'CpG ratio: {0}'.format(min_obs_exp_cpg_ratio_str))
            return
        self.model.compute_islands(
            SeqRecord(seq), island_size, min_gc_ratio,
            min_obs_exp_cpg_ratio, algo_index)

    def _file_loaded(self, file_path):
        """Called when the user loads a file.

        :param file_path: the path to the file
        :type file_path: :class:`str`
        """
        try:
            self.view.set_seq(self.model.load_file(file_path))
        except ValueError as error:
            self.view.show_error('Sequence parsing error: {0}'.format(error))


class ResultsPresenter(object):
    def __init__(self, model, view):
        """Constructor.

        :param model: results model
        :type model: :class:`MetaResultsModel`
        :param view: results view
        :type view: :class:`MetaResultsView`
        """
        self.model = model
        self.view = view

    def register_for_events(self):
        self.model.islands_computed.append(self._islands_computed)
        self.view.island_selected.append(self._island_selected)

    def _seq_feature_to_tuple(self, seq_feature):
        """Convert a SeqFeature to a tuple.

        :param seq_feature: the feature
        :type seq_feature: :class:`SeqFeature`
        :return: the tuple
        :rtype: :class:`tuple` of :class:`int` of length 2
        """
        return (seq_feature.location.start.position,
                seq_feature.location.end.position)

    def _islands_computed(self, seq_record, algo_name, exec_time):
        """Called after island locations have been computed.

        :param seq_record: seq record with features
        :type seq_record: :class:`Bio.SeqRecord.SeqRecord`
        :param algo_name: the name of the algorithm used
        :type algo_name: :class:`str`
        :param exec_time: the algorithm's execution time
        :type exec_time: :class:`float`
        """
        island_tuples = [self._seq_feature_to_tuple(f) for f in
                         seq_record.features]
        self.view.clear_local_seq()
        self.view.clear_global_seq()
        self.view.set_islands(island_tuples)
        self.view.set_algo_name(algo_name)
        self.view.set_exec_time('{0} seconds'.format(exec_time))

    def _island_selected(self, island_index):
        """Called when the selected island is changed.

        :param island_index: index of the requested island
        :type island_index: :class:`int`
        """
        seq_record = self.model.get_seq_record()
        island = seq_record.features[island_index]
        island_tuple = self._seq_feature_to_tuple(island)
        self.view.set_global_seq(str(seq_record.seq), island_tuple)
        self.view.set_local_seq(str(island.extract(seq_record.seq)))
