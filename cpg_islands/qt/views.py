""":mod:`cpg_islands.views.qt` --- Views based on Q toolkit
"""

from PySide import QtGui

from cpg_islands.views import BaseApplicationView


class ApplicationView(QtGui.QMainWindow, BaseApplicationView):
    def start(self):
        """Show and raise the window."""
        self.show()
        self.raise_()
