# -*- coding: utf-8 -*-

# The encoding declaration is necessary to tell Python that we are
# using Unicode in this file. See
# <http://www.python.org/dev/peps/pep-0263/>.
""":mod:`cpg_islands.views.qt` --- Views based on Q toolkit
"""

from PySide import QtGui, QtCore

from cpg_islands import metadata
from cpg_islands.views import (BaseAppView,
                               BaseSeqInputView,
                               BaseResultsView)


class SeqTextEdit(QtGui.QPlainTextEdit):
    """Text editor specifically for DNA sequences. Sets up special
    formatting.
    """
    def __init__(self, parent=None):
        super(SeqTextEdit, self).__init__(parent)
        font = QtGui.QFont('Consolas', 16)
        font.setStyleHint(QtGui.QFont.TypeWriter)
        self.setFont(font)


class AppView(QtGui.QMainWindow, BaseAppView):
    def __init__(self, seq_input_view, results_view, parent=None):
        """Initialize the main application view with docked
        SeqenceInputView and ResultsView.

        :param seqe_input_view: the input view
        :type seq_input_view: :class:`SeqInputView`
        :param results_view: the results view
        :type results_view: :class:`ResultsView`
        """
        super(AppView, self).__init__(parent)

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
        self.input_results_tab_widget = QtGui.QTabWidget(self)
        self.input_results_tab_widget.addTab(seq_input_view, '&Sequence Input')
        self.input_results_tab_widget.addTab(results_view, '&Results')
        self.setCentralWidget(self.input_results_tab_widget)

    def start(self):
        self.showMaximized()
        self.raise_()

    def show_results(self):
        self.input_results_tab_widget.setCurrentIndex(1)

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
        self.form_layout.addRow(u'Island Size ≥', self.island_size_input)

        self.min_gc_ratio_input = QtGui.QLineEdit(self)
        self.min_gc_ratio_validator = QtGui.QDoubleValidator()
        self.min_gc_ratio_validator.setBottom(0)
        self.min_gc_ratio_validator.setTop(1)
        self.min_gc_ratio_input.setValidator(self.min_gc_ratio_validator)
        self.form_layout.addRow(u'GC Ratio ≥', self.min_gc_ratio_input)

        self.min_obs_exp_cpg_ratio_input = QtGui.QLineEdit(self)
        self.min_obs_exp_ratio_validator = QtGui.QDoubleValidator()
        self.min_obs_exp_ratio_validator.setBottom(0)
        self.min_obs_exp_cpg_ratio_input.setValidator(
            self.min_obs_exp_ratio_validator)
        self.form_layout.addRow(
            u'Observed/Expected CpG Ratio ≥', self.min_obs_exp_cpg_ratio_input)

        self.algorithms_combo_box = QtGui.QComboBox(self)
        self.form_layout.addRow('Algorithm', self.algorithms_combo_box)

        self.top_layout.addLayout(self.form_layout)

        self.seq_input_label = QtGui.QLabel('Sequence')
        self.seq_input_label.setAlignment(QtCore.Qt.AlignHCenter)
        self.top_layout.addWidget(self.seq_input_label)

        self.seq_input = SeqTextEdit(self)
        self.seq_input.setTabChangesFocus(True)
        self.top_layout.addWidget(self.seq_input)

        self.submit_button = QtGui.QPushButton('Compute Islands', self)
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

        self.layout = QtGui.QVBoxLayout(self)

        self.form_layout = QtGui.QFormLayout()
        self.exec_time_label = QtGui.QLabel()
        self.exec_time_label.setTextInteractionFlags(
            QtCore.Qt.TextSelectableByMouse |
            QtCore.Qt.TextSelectableByKeyboard)
        self.form_layout.addRow('Execution Time:', self.exec_time_label)
        self.layout.addLayout(self.form_layout)

        self.islands_list = QtGui.QListWidget(self)
        self.islands_list.currentRowChanged.connect(self._island_selected)
        self.islands_list.setFrameShape(QtGui.QFrame.StyledPanel)
        self.global_seq = SeqTextEdit(self)
        self.global_seq.setReadOnly(True)

        self.local_seq = SeqTextEdit(self)
        self.local_seq.setReadOnly(True)

        self.sequences = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.sequences.addWidget(self.global_seq)
        self.sequences.addWidget(self.local_seq)
        self.sequences.setSizes([100, 100])

        self.holder = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.holder.addWidget(self.islands_list)
        self.holder.addWidget(self.sequences)
        self.holder.setSizes([100, 490])

        self.layout.addWidget(self.holder, 1)

    def set_islands(self, islands):
        self.islands_list.clear()
        for start, end in islands:
            self.islands_list.addItem('{0}, {1}'.format(start, end))

    def set_exec_time(self, exec_time):
        self.exec_time_label.setText(exec_time)

    def _island_selected(self, current_row):
        # According to the Qt docs, "If there is no current item, the
        # currentRow is -1." Handle this case.
        if current_row >= 0:
            self.island_selected(current_row)

    def set_local_seq(self, seq_str):
        self.local_seq.setPlainText(seq_str)

    def clear_local_seq(self):
        self.local_seq.clear()

    def set_global_seq(self, seq_str, island_location):
        self.global_seq.setPlainText(seq_str)
        # Extra selections are already set in the constructor, so we
        # know the list has at least a length of one.
        cursor = self.global_seq.textCursor()
        cursor.setPosition(island_location[0])
        cursor.setPosition(island_location[1], QtGui.QTextCursor.KeepAnchor)
        extra_sel = QtGui.QTextEdit.ExtraSelection()
        extra_sel.format.setBackground(QtCore.Qt.green)
        extra_sel.format.setFontOverline(True)
        extra_sel.format.setFontUnderline(True)
        extra_sel.cursor = cursor
        self.global_seq.setExtraSelections([extra_sel])

    def clear_global_seq(self):
        self.global_seq.clear()


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
