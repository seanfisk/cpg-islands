Todo
====

.. todolist::

.. todo::

   Add lots more tests for algorithms.

.. todo::

   Re-factor the Entrez model to something more elegant. It's sort of
   a mess.

.. todo::

   Cython sliding window and accumulating sliding window.

.. todo::

   Rename ``file_loaded`` event to something that makes more sense.

.. todo::

   Loading a file or a sequence from Entrez should prompt if
   overwriting a sequence currently typed in.

.. todo::

   Handle no Internet on the Entrez view.

.. todo::

   Switch assert mock calls to Gray's new style.

.. todo::

   Throttle calls to ESpell.
   
.. todo::

   Refactor error shower code as shown in the Presenter First paper.

.. todo::

   Measure speed of algorithm. Decided to add a "Timing" tab which
   gets populated in the same way as Results when ``compute_islands``
   is called. We will compose the SeqInputModel with a
   TimingModel. The SeqInputModel calls a setter on the TimingModel
   which then sends the information to the TimingView. The timing view
   will list the algorithms run. The algorithms which are run can be
   chosen from a multiple select box within the SeqInputModel. The
   timings will initially be graphed using a bar graph for single
   runs, and then be changed to a box plot for multiple runs of the
   same algorithm.

.. todo::

   Decrease the cyclomatic complexity of the Python accumulator algorithm.

.. todo::

   Refine global sequence text edit to have an auto-zoom.

.. todo::

   C-based extension model implementing sliding window.

.. todo::

   Refine global sequence text edit with numbered lines, better
   format, etc.

.. todo::

   In the open file dialog, the dialog should remember the directory
   last visited. To do this between opens, use a variable. To do this
   between runs (better), store it in a QSettings instance.

.. todo::

   Figure out why flake8>1.5 is not working with NOQA comments. For
   now we are just staying at flake8==1.5.

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
