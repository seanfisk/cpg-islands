==================
CpG Island Locator
==================

.. image:: https://secure.travis-ci.org/seanfisk/cpg-islands.png
   :target: https://secure.travis-ci.org/seanfisk/cpg-islands

-------
Authors
-------
* Sean Fisk
* Gray Gwizdz

------------
Description
------------
A CpG island locator using written in the Python_ programming
language.

This is graduate research project associated with CIS 661 with
Professor Jonathan Leidig at `Grand Valley State University`_.

.. _Python: http://python.org/
.. _Grand Valley State University: http://www.gvsu.edu/

-------
Install
-------

The program's graphical interface depends upon PySide_, Python
bindings to the `Qt`_ libraries. Install this first.

.. code::

    python setup.py install
    cpg_islands

.. _PySide: http://www.pyside.org
.. _Qt: http://www.qt-project.org/

-----------
Development
-----------

Run
===

The program can be run locally by following these steps:

- Install PySide_.
- ``pip install -r requirements/dev.txt``
- ``shovel qt``

Documentation
=============

Documentation is generated on `Read the Docs`_ at
https://cpg-islands.readthedocs.org/. Build it locally by running the following::

    cd docs
    make html

Then open ``build/html/index.html``.

.. _Read the Docs: https://readthedocs.org/

Test
====

Tests are written using pytest_ and mock_. To run unit tests and PEP8
enforcement, run::

    shovel test_all

Continuous integration is provided by Travis-CI_. Find the build
status at https://secure.travis-ci.org/#!/seanfisk/cpg-islands.

.. _pytest: http://pytest.org/
.. _mock: http://www.voidspace.org.uk/python/mock/
.. _Travis-CI: https://travis-ci.org/

Test Coverage
=============

To view the test coverage report, run::

    shovel coverage

-------
License
-------

.. image:: http://www.gnu.org/graphics/gplv3-127x51.png
   :target: `GNU General Public License version 3`_

CpG Islands is free software licensed under the `GNU General Public
License version 3`_.

.. _GNU General Public License version 3: http://www.gnu.org/licenses/gpl.html#content

-------
Credits
-------

CpG Islands makes use of the following libraries/tools/services:

- Python_ programming language
- Biopython_ for sequence fetching and parsing
- Qt_ for graphical interface
- PySide_ for Python bindings to Qt
- git_ version control system
- GitHub_ for git hosting
- pytest_ test framework
- mock_ for creating mock objects
- coverage.py_ and pytest-cov_ for test coverage statistics
- flake8_ as a lint tool: enforced PEP8_ compliance, running PyFlakes_ and `McCabe
  complexity check`_
- Travis-CI_ for continuous integration
- Sphinx_ and docutils_ for documentation generation
- `Read the Docs`_ for documentation hosting
- shovel_ for running miscellaneous tasks

.. _Biopython: http://biopython.org/
.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _Sphinx: http://sphinx.pocoo.org/
.. _docutils: http://docutils.sourceforge.net/
.. _coverage.py: http://nedbatchelder.com/code/coverage/
.. _pytest-cov: http://pypi.python.org/pypi/pytest-cov
.. _flake8: http://pypi.python.org/pypi/flake8
.. _PEP8: https://github.com/jcrocholl/pep8/
.. _PyFlakes: http://pypi.python.org/pypi/pyflakes
.. _McCabe complexity check: http://nedbatchelder.com/blog/200803/python_code_complexity_microtool.html
.. _shovel: https://github.com/seomoz/shovel
