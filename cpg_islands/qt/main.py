#!/usr/bin/env python
""":mod:`cpg_islands.qt.main` -- CpG Islands Locator Qt interface
"""

from cpg_islands.qt.composers import create_presenter


def main(argv=None):
    presenter = create_presenter()

if __name__ == '__main__':
    main()
