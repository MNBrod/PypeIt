"""
Script to install user extinction file into the PypeIt cache.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase
from pypeit import cache

class InstallExtinctfile(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Script to install user-created extinction file',
                                    width=width)
        parser.add_argument('files', type=str, nargs='+',
                            help='Filename(s) of the extinction file(s) to be installed '
                                 'in the PypeIt cache')
        return parser

    @staticmethod
    def main(args):
        import os

        # Loop through the files passed
        for file in args.files:

            # Copy the user-created file to the cache
            cache.write_file_to_cache(file, os.path.basename(file), 'extinction')
