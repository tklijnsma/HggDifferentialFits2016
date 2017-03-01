#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, glob, re
from os.path import *


########################################
# Main
########################################

def main():

    baseDir = abspath( join( os.environ['CMSSW_BASE'], 'src' ) )
    os.chdir( baseDir )
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--asjob', action='store_true', help='boolean')
    parser.add_argument( 'jobdirs', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()



    for jobdir in args.jobdirs:

        try:

            jobdir = abspath( jobdir )
            os.chdir( jobdir )

            seed = re.search( r'seed(\d+)', jobdir ).group(1)

            cmd = ''
            cmd += ' combine  higgsCombinebestfitToy.MultiDimFit.mH125.root  '
            cmd += ' -M MultiDimFit  '
            cmd += ' -n REDOcovmat  '
            cmd += ' -m 125  '
            cmd += ' --saveWorkspace   '
            cmd += ' --dataset higgsCombinegeneratedToy.GenerateOnly.mH125.{0}.root:toys/toy_1  '.format(seed)
            cmd += ' --algo none  '
            cmd += ' --snapshotName MultiDimFit  '
            # cmd += ' --setPhysicsModelParameters r0_cat0=1,r1_cat0=1,r2_cat0=1,r3_cat0=1,r4_cat0=1,r5_cat0=1,r6_cat0=1,r0_cat1=1,r1_cat1=1,r2_cat1=1,r3_cat1=1,r4_cat1=1,r5_cat1=1,r6_cat1=1,r0_cat2=1,r1_cat2=1,r2_cat2=1,r3_cat2=1,r4_cat2=1,r5_cat2=1,r6_cat2=1  '
            cmd += ' --freezeNuisances "$(python  /afs/cern.ch/work/t/tklijnsm/Regularization/CMSSW_7_4_7/src/getCommaSeparatedListOfVariablesToFreeze.py  higgsCombinebestfitToy.MultiDimFit.mH125.root )"  '
            cmd += ' -v 2 '


            print '\n\nIn {0}, doing:'.format(jobdir)
            print cmd

            if not args.asjob:
                os.system(cmd)
            else:

                shFile = 'covmatRedo.sh'
                with open( shFile, 'w' ) as shFp:
                    w = lambda text: shFp.write( text + '\n' )
                    w( 'cd {0}'.format( baseDir ) )
                    w( 'eval `scramv1 runtime -sh`' )
                    w( 'cd {0}'.format( jobdir ) )
                    w( cmd )

                os.system( 'chmod 777 {0}'.format(shFile) )
                os.system( 'bsub {0}'.format(shFile) )

        finally:

            os.chdir( baseDir )








########################################
# End of Main
########################################
if __name__ == "__main__":
    main()