Todo
====

.. todolist::

.. todo::

   BUG: Analyze one sequence with at least one island. Then analyze a
   sequence without any islands. In the results view, the global and
   local sequence are still set to the old sequences.

.. todo::

   BUG: When Cancel is clicked in the file open dialog, an IOError (no
   such file or directory) is printed to the console. This should go
   to the user or just plain be ignored. Either way, it should be
   handled and should *not* be sent to the console.

.. todo::

   "Load File" should change to "Open" and should be given the
   shortcut Ctrl/Cmd-O.

.. todo::

   Fix sequence input text area in SeqInputView to expand to its
   maximum size.
   
.. todo::

   Sean - Remove some redundancy in valid algorithm tests in ``test_results.py``.

.. todo::

   Sean - Have mutiple algorithm options available.

.. todo::

   Refine sequence highlight to look nicer, use fixed with text, etc.

.. todo::

   Gray - ResultsView:

   * list of CpG islands on the left
   * big picture view of sequence on the top right
   * detailed view of sequence on the bottom right

.. todo::

   Measure speed of algorithm.

.. todo::

   C-based extension model implementing sliding window.

.. todo::

   Entrez database search (new MVP triad)

.. todo::

   Sean - Distribution with pyinstaller.

.. todo::

   Make documentation organization nicer.

.. todo::

   Tabbed interface for sequence files.

.. todo::

   Analysis run history.

.. todo::

   Load fasta sequence from file

.. todo::

   Perform "client-side" validation for the Sequence box:

   * Prevent from typing unwanted characters
   * Prevent invalid pastes - **rejected on the basis that some
     cleanup of the sequence might be needed**
