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
        self.model.seq_loaded.append(self.view.show_seq_input)
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

    def _format_zero_indexed(self, index):
        """Format an index with a note that the index is zero-indexed.

        :param index: the index to format
        :type index: :class:`int`
        """
        return '{0} (zero-indexed)'.format(index)

    def _islands_computed(self, global_seq, feature_tuples, algo_name,
                          exec_time):
        """Called after island locations have been computed.

        :param global_seq: the full sequence
        :type global_seq: :class:`str`
        :param feature_tuples: tuples containing locations of the features
        :type feature_tuples: :class:`list` of :class:`tuple`
        :param algo_name: the name of the algorithm used
        :type algo_name: :class:`str`
        :param exec_time: algorithm's execution duration
        :type exec_time: :class:`float`
        """
        self.view.clear_subseq()
        self.view.set_global_seq(global_seq)
        self.view.set_islands(feature_tuples)
        self.view.set_algo_name(algo_name)
        self.view.set_exec_time('{0} seconds'.format(exec_time))

    def _island_selected(self, island_index):
        """Called when the selected island is changed.

        :param island_index: index of the requested island
        :type island_index: :class:`int`
        """
        island_info = self.model.get_island_info(island_index)
        self.view.highlight_global_seq(island_info.start, island_info.end)
        self.view.set_start(self._format_zero_indexed(island_info.start))
        self.view.set_end(self._format_zero_indexed(island_info.end))
        self.view.set_length('{0} bases'.format(island_info.length))
        self.view.set_subseq(island_info.subseq)
        self.view.set_gc_ratio('{0}%'.format(island_info.gc_ratio * 100))
        self.view.set_obs_exp_cpg_ratio(str(island_info.obs_exp_cpg_ratio))


class EntrezPresenter(object):
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
        """Connect view methods to presenter methods."""
        self.view.query_changed.append(self._query_changed)
        self.view.search_requested.append(self._user_submits)
        self.view.result_selected.append(self._user_selected)
        self.view.load_requested.append(self._load_selected)

    def _user_submits(self, text):
        """Handle user submission.

        :param text: text to search
        :type text: :class:`str`
        """
        id_list, query_translation = self.model.search(text)
        self.view.set_result(id_list)
        self.view.set_query_translation(query_translation)

    def _user_selected(self, index):
        """Handle user submission.

        :param index: list index of selected item on view
        :type index: :class:`int`
        """
        result = self.model.get_seq(index)
        self.view.set_selected_seq(str(result.seq))

    def _query_changed(self, query):
        """Handle user suggestions.

        :param text: text to query for suggestion
        :type text: :class:`str`
        """
        suggestion = self.model.suggest(query)
        self.view.set_suggestion(suggestion)

    def _load_selected(self):
        """Handle loading sequences."""
        self.model.load_seq()
