""":mod:`cpg_islands.models` --- Application models
"""

from __future__ import print_function
import sys
import argparse

from cpg_islands import metadata


class ApplicationModel(object):
    def run(self, argv=None):
        if argv is None:
            argv = sys.argv

        author_strings = []
        for name, email in zip(metadata.authors, metadata.emails):
            author_strings.append('Author: {0} <{1}>'.format(name, email))

        epilog = '''{title} {version}

{authors}
URL: <{url}>
'''.format(
            title=metadata.nice_title,
            version=metadata.version,
            authors='\n'.join(author_strings),
            url=metadata.url)

        arg_parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=metadata.description,
            epilog=epilog)

        args = arg_parser.parse_args(args=argv[1:])

        print(epilog)
