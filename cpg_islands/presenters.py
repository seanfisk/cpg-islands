from Bio.Seq import Seq

class ApplicationPresenter(object):
    def __init__(self, model, view):
        self.model = model
        self.view = view

    def register_for_events(self):
        self.model.started.append(self.view.start)

    def _user_submits(self, seq_str, island_size_str, minimum_gc_ratio_str):
        """Called when the user submits the form.

        :param seq_str: the sequence as a string
        :type seq_str: :class:`str`
        :param island_size_str: number of bases which an island may contain
        :type island_size_str: :class:`str`
        :param minimum_gc_ratio_str: the ratio of GC to other bases
        :type minimum_gc_ratio_str: :class:`str`
        """
        seq = Seq(seq_str)
        island_size = int(island_size_str)
        minimum_gc_ratio = float(minimum_gc_ratio_str)
        locations = self.model.annotate_cpg_islands(seq,
                                                    island_size,
                                                    minimum_gc_ratio)
        self.view.set_locations(locations)
