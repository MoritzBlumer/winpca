"""
WinPCA. A package for windowed PC analysis.
"""

## IMPORT PACKAGES
import sys


## CLASSES

class Log():
    '''
    Logger.
    '''

    def __init__(self):
        pass

    def newline(self):
        '''
        Print ERROR message to STDERR and exit.
        '''

        print('',
              file=sys.stderr,
              flush=True,
              )


    def info(self, message):
        '''
        Print INFO message to STDERR.
        '''

        print(f'[INFO] {message}',
              file=sys.stderr,
              flush=True,
              )


    def error(self, message):
        '''
        Print ERROR message to STDERR and exit.
        '''

        print(f'[ERROR] {message}\n',
              file=sys.stderr,
              flush=True,
              )
        sys.exit(1)

    def error_nl(self, message):
        '''
        Print ERROR message to STDERR and exit.
        '''

        print(f'\n[ERROR] {message}\n',
              file=sys.stderr,
              flush=True,
              )
        sys.exit(1)
