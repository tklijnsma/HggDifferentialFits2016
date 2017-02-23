#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, shutil
from os.path import *

from time import strftime
datestr = strftime( '%b%d' )

import argparse

########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'resultdir', type=str, default='default', help='default string' )
    parser.add_argument( '--suffix', type=str, help='default string' )
    parser.add_argument( '--copy', action='store_true', help='boolean')
    args = parser.parse_args()


    outdir = 'results_{0}_{1}'.format( datestr, dirname(args.resultdir).replace('/','') )
    if args.suffix:
        outdir += '_' + args.suffix

    if isdir( outdir ):
        outdir += '_{0}'
        i = 1
        while isdir( outdir.format(i) ):
            i += 1
        outdir = outdir.format(i)

    os.makedirs( outdir )



    allFiles = os.listdir( args.resultdir )

    for f in allFiles:

        if basename(f).startswith('CMS-HGG'): continue

        if args.copy or f.endswith('.txt'):
            shutil.copyfile(
                join( args.resultdir, basename(f) ),
                join( outdir, basename(f) )
                )
        else:
            os.rename(
                join( args.resultdir, basename(f) ),
                join( outdir, basename(f) )
                )











########################################
# End of Main
########################################
if __name__ == "__main__":
    main()