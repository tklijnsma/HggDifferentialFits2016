#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os


########################################
# Main
########################################

def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( 'txtFile', type=str, default='default', help='default string' )
    # parser.add_argument( '--boolean', action='store_true', help='boolean')
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()

    with open( args.txtFile, 'r' ) as txtFp:
        lines = txtFp.readlines()


    output = []

    for line in lines:
        line = line.strip()

        if len(line) == 0:
            continue
        elif line.startswith('#'):
            continue

        output.append( line )


    outputStr = ','.join( output )


    print outputStr


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()