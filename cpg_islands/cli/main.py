#!/usr/bin/env python
""":mod:`cpg_islands.cli.main` -- CpG Islands Locator command-line interface
"""

from cpg_islands.cli.composers import create_presenter


def main(argv=None):
    presenter = create_presenter()

if __name__ == '__main__':
    main()
