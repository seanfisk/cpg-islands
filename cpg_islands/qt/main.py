#!/usr/bin/env python
""":mod:`cpg_islands.qt.main` -- CpG Islands Locator Qt interface
"""

import sys

from PySide import QtGui

from cpg_islands import metadata
from cpg_islands.qt.composers import create_application_presenter


def main(argv=None):
    if argv is None:
        argv = sys.argv

    app = QtGui.QApplication(argv)
    app.setOrganizationName(metadata.organization)
    app.setOrganizationDomain(metadata.organization_domain)
    app.setApplicationName(metadata.nice_title)
    app.setApplicationVersion(metadata.version)

    presenter = create_application_presenter(argv)

    return app.exec_()

if __name__ == '__main__':
    sys.exit(main())
