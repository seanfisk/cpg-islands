Todo
====

.. todolist::

.. todo::

   Sean - Improve/make docstrings consistent.      

.. todo::

   Remove **Show this message again** from the error boxes.

.. todo::

   Sean - Create new MVP triads:

   * Application{Model,View,Presenter}
   * SequenceInput{Model,View,Presenter}
   * Results{Model,View,Presenter}

    Initially, these will be implemented as different
    windows. Specifically, the ResultsView will have to be a modal dialog.

.. todo::

   In mock tests, switch from using ``assert_called_once_with()`` to
   testing *all* mock calls with ``mock_calls``. This will ensure that
   expected functions are called and unexpected functions are *not* called.
    
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

   Sean - Distribution with pyinstaller.

.. todo::

   Make documentation organization nicer.

.. todo::

   Load fasta sequence from file

.. todo::

   Perform "client-side" validation for the Sequence box:

   * Prevent from typing unwanted characters
   * Prevent invalid pastes - **rejected on the basis that some
     cleanup of the sequence might be needed**
