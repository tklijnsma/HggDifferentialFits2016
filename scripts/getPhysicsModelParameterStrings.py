#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import argparse

########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'nBins', type=int, help='default string' )
    parser.add_argument( 'nCats', type=int, default=0, help='default string', nargs='?' )
    args = parser.parse_args()

    ret = []

    if args.nCats == 0:
        for iBin in xrange(args.nBins):
            ret.append( 'r{0}=1'.format(iBin) )
    else:
        for iCat in xrange(args.nCats):
            for iBin in xrange(args.nBins):
                ret.append( 'r{0}_cat{1}=1'.format(iBin,iCat) )

    print ','.join(ret)


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()