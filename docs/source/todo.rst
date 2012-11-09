Todo
====

.. todolist::

.. todo::

   Remove some redundancy in valid algorithm tests in ``test_results.py``.
   
.. todo::

   Gray - ResultsView:

   * list of CpG islands on the left
   * big picture view of sequence on the top right
   * detailed view of sequence on the bottom right

.. todo::

   Dock the SequenceInputView and ResultsView into the ApplicationView.
     
.. todo::

   Entrez database search (new MVP triad)

.. todo::

   New algorithm.

.. todo::

   Measure speed of algorithm.

.. todo::

   Have mutiple algorithm options available.

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

.. todo::

   Currently the graphics portion of the application is initialized
   before command-line arguments are parsed. This causes focus to
   "flash" momentarily (at least on OS X) while a windowing
   environment is created then destroyed. This is obvious when passing
   a --help or --version argument that shouldn't use any
   graphics. This is annoying, but low priority.
