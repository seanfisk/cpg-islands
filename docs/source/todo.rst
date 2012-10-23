Todo
====

.. todolist::

.. todo::

   Add CpG island annotator to the application model.

   .. function:: annotate_cpg_islands(seq, island_size, minimum_gc_ratio)

      Island finder function.

      :param seq: the sequence to annotate
      :type seq: :class:`Bio.SeqRecord`
      :param island_size: the number of bases which an island may contain
      :type island_size: :class:`int`
      :param minimum_gc_ratio: the ratio of GC to other bases
      :type minimum_gc_ratio: :class:`float`
      :return: not sure yet

.. todo::

   Basic tested command-line application.