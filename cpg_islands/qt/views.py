# -*- coding: utf-8 -*-

# The encoding declaration is necessary to tell Python that we are
# using Unicode in this file. See
# <http://www.python.org/dev/peps/pep-0263/>.
""":mod:`cpg_islands.views.qt` --- Views based on Q toolkit
"""

from xml.sax import saxutils

from PySide import QtGui, QtCore

from cpg_islands import metadata
from cpg_islands.views import (BaseAppView,
                               BaseSeqInputView,
                               BaseResultsView,
                               BaseEntrezView)


class LeftFormLayout(QtGui.QFormLayout):
    """A form whose contents are left-aligned. This overrides the
    default Mac OS X behavior for QFormLayout, however, it looks more
    natural with our application."""
    def __init__(self, parent=None):
        super(LeftFormLayout, self).__init__(parent)
        self.setFormAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)


class SeqTextEdit(QtGui.QPlainTextEdit):
    """Text editor specifically for DNA sequences. Sets up special
    formatting.
    """
    def __init__(self, parent=None):
        super(SeqTextEdit, self).__init__(parent)
        font = QtGui.QFont('Consolas', 16)
        font.setStyleHint(QtGui.QFont.TypeWriter)
        self.setFont(font)


class StatsLabel(QtGui.QLabel):
    """Label for displaying statistics."""
    def __init__(self, parent=None):
        super(StatsLabel, self).__init__(parent)
        self.setTextInteractionFlags(
            QtCore.Qt.TextSelectableByMouse |
            QtCore.Qt.TextSelectableByKeyboard)


class AppView(QtGui.QMainWindow, BaseAppView):
    def __init__(self, entrez_view, seq_input_view, results_view, parent=None):
        """Initialize the main application view with docked
        SeqenceInputView and ResultsView.

        :param entrez_view: the Entrez search view
        :type entrez_view: :class:`EntrezView`
        :param seq_input_view: the input view
        :type seq_input_view: :class:`SeqInputView`
        :param results_view: the results view
        :type results_view: :class:`ResultsView`
        """
        super(AppView, self).__init__(parent)

        # Window
        self.setWindowTitle(metadata.nice_title)

        # Menu
        self.menu_bar = QtGui.QMenuBar()
        self.file_menu = self.menu_bar.addMenu('&File')
        self.file_action = self.file_menu.addAction('&Open...')
        self.file_action.setShortcut(QtGui.QKeySequence.Open)
        self.file_action.triggered.connect(self._load_file)
        self.quit_action = self.file_menu.addAction('&Quit')
        self.quit_action.setShortcut(QtGui.QKeySequence.Quit)
        self.quit_action.triggered.connect(self.close)
        self.help_menu = self.menu_bar.addMenu('&Help')
        self.about_action = self.help_menu.addAction('&About')
        self.about_action.triggered.connect(self._about)
        self.setMenuBar(self.menu_bar)

        # Main
        self.entrez_view = entrez_view
        self.seq_input_view = seq_input_view
        self.results_view = results_view
        self.tab_widget = QtGui.QTabWidget(self)
        self.tab_widget.addTab(self.entrez_view, 'Entre&z Search')
        self.tab_widget.addTab(self.seq_input_view, '&Sequence Input')
        self.tab_widget.addTab(self.results_view, 'Analysis &Results')
        self.tab_widget.setCurrentWidget(self.seq_input_view)
        self.setCentralWidget(self.tab_widget)

    def start(self):
        self.showMaximized()
        self.raise_()

    def show_results(self):
        self.tab_widget.setCurrentWidget(self.results_view)

    def show_seq_input(self):
        self.tab_widget.setCurrentWidget(self.seq_input_view)

    def _about(self):
        """Create and show the about dialog."""
        AboutDialog(self).exec_()

    def _load_file(self):
        """Create and show the file dialog."""
        # `getOpenFileName' returns a tuple of (file_name, filter).
        file_name_and_filter_tuple = QtGui.QFileDialog.getOpenFileName(
            self,
            caption='Load GenBank File...',
            filter='GenBank Sequence File (*.gb)')

        # According to the Qt documentation on `getOpenFileName'', "If
        # the user presses Cancel, it returns a null string." That's
        # not quite true for PySide, where it will return a tuple of
        # two empty strings. Handle it here and don't let it reach the
        # file loading code.
        file_name = file_name_and_filter_tuple[0]
        if len(file_name) > 0:
            self.file_load_requested(file_name)


class SeqInputView(QtGui.QWidget, BaseSeqInputView):
    def __init__(self, parent=None):
        super(SeqInputView, self).__init__(parent)

        self.top_layout = QtGui.QVBoxLayout(self)
        self.form_layout = QtGui.QFormLayout()

        self.island_size_input = QtGui.QLineEdit(self)
        self.island_size_validator = QtGui.QIntValidator()
        self.island_size_validator.setBottom(0)
        self.island_size_input.setValidator(self.island_size_validator)
        self.form_layout.addRow(u'&Island Size ≥', self.island_size_input)

        self.min_gc_ratio_input = QtGui.QLineEdit(self)
        self.min_gc_ratio_validator = QtGui.QDoubleValidator()
        self.min_gc_ratio_validator.setBottom(0)
        self.min_gc_ratio_validator.setTop(1)
        self.min_gc_ratio_input.setValidator(self.min_gc_ratio_validator)
        self.form_layout.addRow(u'&GC Ratio ≥', self.min_gc_ratio_input)

        self.min_obs_exp_cpg_ratio_input = QtGui.QLineEdit(self)
        self.min_obs_exp_ratio_validator = QtGui.QDoubleValidator()
        self.min_obs_exp_ratio_validator.setBottom(0)
        self.min_obs_exp_cpg_ratio_input.setValidator(
            self.min_obs_exp_ratio_validator)
        self.form_layout.addRow(u'&Observed/Expected CpG Ratio ≥',
                                self.min_obs_exp_cpg_ratio_input)

        self.algorithms_combo_box = QtGui.QComboBox(self)
        self.form_layout.addRow('&Algorithm', self.algorithms_combo_box)

        self.top_layout.addLayout(self.form_layout)

        self.seq_input_label = QtGui.QLabel('S&equence')
        self.seq_input_label.setAlignment(QtCore.Qt.AlignHCenter)
        self.top_layout.addWidget(self.seq_input_label)

        self.seq_input = SeqTextEdit(self)
        self.seq_input.setTabChangesFocus(True)
        self.seq_input_label.setBuddy(self.seq_input)
        self.top_layout.addWidget(self.seq_input)

        self.submit_button = QtGui.QPushButton('&Compute Islands', self)
        self.submit_button.clicked.connect(self._submit_clicked)
        self.top_layout.addWidget(self.submit_button)

    def _get_seq(self):
        """Return the widget's entered text.

        :return: the text
        :rtype: :class:`str`
        """
        return self.seq_input.toPlainText()

    def set_seq(self, seq_str):
        """Set the sequence text.

        :param seq_str: the sequence in string form
        :type seq_str: :class:`str`
        """
        self.seq_input.setPlainText(seq_str)

    def _get_min_gc_ratio(self):
        """Return the widget's entered GC ratio.

        :return: the ratio
        :rtype: :class:`str`
        """
        return self.min_gc_ratio_input.text()

    def set_min_gc_ratio(self, min_gc_ratio_str):
        self.min_gc_ratio_input.setText(min_gc_ratio_str)

    def _get_min_obs_exp_cpg_ratio(self):
        """Return the widget's entered observed/expected CpG ratio.

        :return: the ratio
        :rtype: :class:`str`
        """
        return self.min_obs_exp_cpg_ratio_input.text()

    def set_min_obs_exp_cpg_ratio(self, min_obs_exp_cpg_ratio_str):
        self.min_obs_exp_cpg_ratio_input.setText(min_obs_exp_cpg_ratio_str)

    def _get_island_size(self):
        """Return the widget's entered island size.

        :return: the key
        :rtype: :class:`str`
        """
        return self.island_size_input.text()

    def _get_algorithm_index(self):
        """Return the currently selected algorithm's index.

        :return: the algorithm index
        :rtype: :class:`int`
        """
        return self.algorithms_combo_box.currentIndex()

    def set_island_size(self, island_size):
        self.island_size_input.setText(island_size)

    def set_algorithms(self, algorithm_names):
        self.algorithms_combo_box.clear()
        self.algorithms_combo_box.addItems(algorithm_names)

    def show_error(self, message):
        """Show the user an error dialog.

        :param message: error message
        :type message: :class:`str`
        """
        QtGui.QMessageBox.critical(self, metadata.nice_title, message)

    def _submit_clicked(self):
        """Submit the input to the model."""
        try:
            self.submitted(self._get_seq(),
                           self._get_island_size(),
                           self._get_min_gc_ratio(),
                           self._get_min_obs_exp_cpg_ratio(),
                           self._get_algorithm_index())
        except ValueError as error:
            self.show_error(str(error))


class ResultsView(QtGui.QWidget, BaseResultsView):
    def __init__(self, parent=None):
        super(ResultsView, self).__init__(parent)

        self.top_layout = QtGui.QVBoxLayout(self)

        # Global statistics
        self.global_stats_layout = LeftFormLayout()
        self.algo_name_label = StatsLabel(self)
        self.exec_time_label = StatsLabel(self)
        self.global_stats_layout.addRow(
            '&Algorithm Used:', self.algo_name_label)
        self.global_stats_layout.addRow(
            '&Execution Time:', self.exec_time_label)
        self.top_layout.addLayout(self.global_stats_layout)

        # Islands list
        self.islands_list_container = QtGui.QWidget(self)
        self.islands_list_label = QtGui.QLabel('Islands &List', self)
        self.islands_list = QtGui.QListWidget(self)
        self.islands_list_label.setBuddy(self.islands_list)
        self.islands_list.currentRowChanged.connect(self._island_selected)
        self.islands_list.setFrameShape(QtGui.QFrame.StyledPanel)
        self.islands_list_layout = QtGui.QVBoxLayout(
            self.islands_list_container)
        self.islands_list_layout.addWidget(self.islands_list_label)
        self.islands_list_layout.addWidget(self.islands_list)

        # Global sequence
        self.global_seq_container = QtGui.QWidget(self)
        self.global_seq_label = QtGui.QLabel(
            'Island &Highlighted Within Full Sequence', self)
        self.global_seq = SeqTextEdit(self)
        self.global_seq_label.setBuddy(self.global_seq_label)
        self.global_seq.setReadOnly(True)
        self.global_seq_layout = QtGui.QVBoxLayout(self.global_seq_container)
        self.global_seq_layout.addWidget(self.global_seq_label)
        self.global_seq_layout.addWidget(self.global_seq)

        # Subsequence
        self.subseq_container = QtGui.QWidget(self)
        self.subseq_layout = QtGui.QVBoxLayout(self.subseq_container)
        ## Subseq stats
        self.subseq_stats_layout = LeftFormLayout()
        self.subseq_start = StatsLabel(self)
        self.subseq_stats_layout.addRow('S&tart Index:', self.subseq_start)
        self.subseq_end = StatsLabel(self)
        self.subseq_stats_layout.addRow('E&nd Index:', self.subseq_end)
        self.subseq_length = StatsLabel(self)
        self.subseq_stats_layout.addRow('Leng&th:', self.subseq_length)
        self.subseq_gc_ratio = StatsLabel(self)
        self.subseq_stats_layout.addRow('&GC Ratio:', self.subseq_gc_ratio)

        self.subseq_obs_exp_cpg_ratio = StatsLabel(self)
        self.subseq_stats_layout.addRow('&Observed/Expected CpG Ratio:',
                                        self.subseq_obs_exp_cpg_ratio)
        self.subseq_layout.addLayout(self.subseq_stats_layout)

        ## Subseq bases
        self.subseq_label = QtGui.QLabel(
            'Island &Bases', self)
        self.subseq = SeqTextEdit(self)
        self.subseq_label.setBuddy(self.subseq_label)
        self.subseq.setReadOnly(True)
        self.subseq_layout.addWidget(self.subseq_label)
        self.subseq_layout.addWidget(self.subseq)

        self.seq_splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.seq_splitter.addWidget(self.global_seq_container)
        self.seq_splitter.addWidget(self.subseq_container)
        self.seq_splitter.setSizes([100, 100])

        self.list_seq_splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.list_seq_splitter.addWidget(self.islands_list_container)
        self.list_seq_splitter.addWidget(self.seq_splitter)
        self.list_seq_splitter.setSizes([100, 490])

        self.top_layout.addWidget(self.list_seq_splitter, 1)

    def set_islands(self, islands):
        self.islands_list.clear()
        for i, (start, end) in enumerate(islands):
            self.islands_list.addItem(
                'Island {0} ({1}-{2})'.format(i, start, end))

    def set_algo_name(self, algo_name):
        self.algo_name_label.setText(algo_name)

    def set_exec_time(self, exec_time_str):
        self.exec_time_label.setText(exec_time_str)

    def set_global_seq(self, seq_str):
        self.global_seq.setPlainText(seq_str)

    def highlight_global_seq(self, start, end):
        # Create a text cursor from the existing text cursor.
        cursor = self.global_seq.textCursor()
        # Set the starting position. Don't do any selection yet.
        cursor.setPosition(start)
        # Now set this cursor as the main cursor for the widget.
        self.global_seq.setTextCursor(cursor)
        # Scroll such that the beginning of the sequence is now in the
        # center of the widget.
        self.global_seq.centerCursor()
        # Now do the selection, and use this cursor for the extra
        # selection, which creates a permanent selection.
        cursor.setPosition(end, QtGui.QTextCursor.KeepAnchor)
        extra_sel = QtGui.QTextEdit.ExtraSelection()
        extra_sel.format.setBackground(QtCore.Qt.green)
        extra_sel.format.setFontOverline(True)
        extra_sel.format.setFontUnderline(True)
        extra_sel.cursor = cursor
        self.global_seq.setExtraSelections([extra_sel])

    def set_start(self, start_str):
        self.subseq_start.setText(start_str)

    def set_end(self, end_str):
        self.subseq_end.setText(end_str)

    def set_length(self, length_str):
        self.subseq_length.setText(length_str)

    def set_subseq(self, seq_str):
        self.subseq.setPlainText(seq_str)

    def clear_subseq(self):
        self.subseq.clear()

    def set_gc_ratio(self, gc_ratio_str):
        self.subseq_gc_ratio.setText(gc_ratio_str)

    def set_obs_exp_cpg_ratio(self, obs_exp_cpg_ratio_str):
        self.subseq_obs_exp_cpg_ratio.setText(obs_exp_cpg_ratio_str)

    def _island_selected(self, current_row):
        # According to the Qt docs, "If there is no current item, the
        # currentRow is -1." Handle this case.
        if current_row >= 0:
            self.island_selected(current_row)


class EntrezView(QtGui.QWidget, BaseEntrezView):
    def __init__(self, parent=None):
        """Construct a entrez widget.

        :param parent: widget parent
        :type parent: :class:`QtGui.QWidget`
        """
        super(EntrezView, self).__init__(parent)

        self.top_layout = QtGui.QVBoxLayout(self)

        # Top Search Form
        self.search_layout = QtGui.QFormLayout()
        self.query_input = QtGui.QLineEdit(self)
        self.query_input.textChanged.connect(self._query_input_changed)
        self.search_layout.addRow('&Query', self.query_input)

        self.suggestion_display = QtGui.QLabel(self)
        self.search_layout.addRow('&Suggested Query', self.suggestion_display)

        self.query_translation_display = QtGui.QLabel(self)
        self.search_layout.addRow(
            'Query &Translation', self.query_translation_display)

        self.submit_button = QtGui.QPushButton('Se&arch', self)
        self.submit_button.clicked.connect(self._submit_clicked)
        self.search_layout.addRow(self.submit_button)
        self.top_layout.addLayout(self.search_layout)

        # Bottom Results Form
        self.results_splitter = QtGui.QSplitter(self)
        ## Left Side Results List
        self.results_list_container = QtGui.QWidget(self)
        self.results_list_layout = QtGui.QVBoxLayout(
            self.results_list_container)
        self.results_list_label = QtGui.QLabel('Res&ults', self)
        self.results_list_layout.addWidget(self.results_list_label)
        self.results_list = QtGui.QListWidget(self)
        self.results_list_label.setBuddy(self.results_list)
        self.results_list.currentRowChanged.connect(self._result_selected)
        self.results_list_layout.addWidget(self.results_list)
        self.results_splitter.addWidget(self.results_list_container)

        self.seq_display_container = QtGui.QWidget(self)
        self.seq_display_layout = QtGui.QVBoxLayout(self.seq_display_container)
        self.seq_info_layout = LeftFormLayout()
        self.seq_locus_display = QtGui.QLabel(self)
        self.seq_locus_display.setOpenExternalLinks(True)
        self.seq_info_layout.addRow(
            QtGui.QLabel(
                "Click the locus to visit the sequence on NCBI's website.",
                self))
        self.seq_info_layout.addRow(
            '&Locus (NCBI identifier):', self.seq_locus_display)
        self.seq_desc_display = QtGui.QLabel(self)
        self.seq_info_layout.addRow('&Description:', self.seq_desc_display)
        self.seq_len_display = QtGui.QLabel(self)
        self.seq_info_layout.addRow('Len&gth:', self.seq_len_display)
        self.seq_display_layout.addLayout(self.seq_info_layout)
        self.seq_display_label = QtGui.QLabel('Seque&nce', self)
        self.seq_display_layout.addWidget(self.seq_display_label)
        self.seq_display = SeqTextEdit(self)
        self.seq_display_label.setBuddy(self.seq_display)
        self.seq_display.setReadOnly(True)
        self.seq_display_layout.addWidget(self.seq_display)
        self.results_splitter.addWidget(self.seq_display_container)
        # Set self.seq_display_layout to have a stretch factor of
        # 1. The idea is that the sequence text area should be bigger
        # than the results list text area.
        self.results_splitter.setStretchFactor(1, 1)

        self.top_layout.addWidget(self.results_splitter, 1)

        self.load_button = QtGui.QPushButton('Load', self)
        self.load_button.clicked.connect(self._load_clicked)
        self.top_layout.addWidget(self.load_button)

    def set_suggestion(self, suggestion):
        self.suggestion_display.setText(suggestion)

    def set_query_translation(self, query_translation):
        self.query_translation_display.setText(query_translation)

    def set_seq_locus(self, locus, ncbi_url):
        # Do the right thing and escape these
        # fields. `saxutils.quoteattr()' automatically adds quotes
        # around the text for you.
        self.seq_locus_display.setText('<a href={0}>{1}</a>'.format(
            saxutils.quoteattr(ncbi_url),
            saxutils.escape(locus)))

    def set_seq_desc(self, desc):
        self.seq_desc_display.setText(desc)

    def set_seq_len(self, len_str):
        self.seq_len_display.setText(len_str)

    def set_selected_seq(self, seq_str):
        self.seq_display.setPlainText(seq_str)

    def set_result(self, results):
        self.results_list.clear()
        self.results_list.addItems(results)

    def _get_query(self):
        """Return the query widget's entered text.

        :return: the text
        :rtype: :class:`str`
        """
        return self.query_input.text()

    def _submit_clicked(self):
        """Submit the entered term."""
        self.search_requested(self._get_query())

    def _load_clicked(self):
        """Submit the entered term."""
        self.load_requested()

    def _query_input_changed(self):
        """Submit the entered term."""
        self.query_changed(self._get_query())

    def _result_selected(self, current_row):
        """Pulls the selected index."""
        if current_row >= 0:
            self.result_selected(current_row)


class AboutDialog(QtGui.QDialog):
    """Shows information about the program."""
    def __init__(self, parent=None):
        """Construct the dialog.

        :param parent: the widget's parent
        :type parent: :class:`QtGui.QWidget`
        """
        super(AboutDialog, self).__init__(parent)
        self.setWindowTitle('About ' + metadata.nice_title)
        self.layout = QtGui.QVBoxLayout(self)
        self.title_label = QtGui.QLabel(metadata.nice_title, self)
        self.layout.addWidget(self.title_label)
        self.version_label = QtGui.QLabel('Version ' + metadata.version, self)
        self.layout.addWidget(self.version_label)
        self.copyright_label = QtGui.QLabel('Copyright (C) ' +
                                            metadata.copyright, self)
        self.layout.addWidget(self.copyright_label)
        self.url_label = QtGui.QLabel(
            'Source: <a href="{0}">{0}</a>'.format(metadata.url), self)
        self.url_label.setOpenExternalLinks(True)
        self.layout.addWidget(self.url_label)
        self.documentation_label = QtGui.QLabel(
            'Documentation: <a href="{0}">{0}</a>'.format(
            metadata.documentation_url), self)
        self.documentation_label.setOpenExternalLinks(True)
        self.layout.addWidget(self.documentation_label)
