#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, shutil, sys
from os.path import *
import argparse
import tempfile
from time import strftime

import ROOT

from getCommaSeparatedListOfVariablesToFreeze import getCommaSeparatedListOfVariablesToFreeze


########################################
# Parser
########################################

parser = argparse.ArgumentParser()

parser.add_argument( '--clean', action='store_true', help='boolean')
parser.add_argument( '--test', action='store_true', help='boolean')
parser.add_argument( '--printCommands', action='store_true', help='boolean')
parser.add_argument( '--local', action='store_true', help='boolean')

parser.add_argument( '--moreMemory', action='store_true', help='boolean')
parser.add_argument( '--freezeAllPdfVars', action='store_true', help='boolean')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument( '--doT2WS', action='store_true', help='boolean')
group.add_argument( '--dobestfit', action='store_true', help='boolean')
group.add_argument( '--dotoys', action='store_true', help='boolean')

parser.add_argument( '--nToys', type=int, default=1, help='default string' )
parser.add_argument( '--startseed', type=int, help='default string' )



group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument( '--dopt', action='store_true', help='boolean')
group2.add_argument( '--donjets', action='store_true', help='boolean')

parser.add_argument( '--newbins', action='store_true', default=False, help='boolean')

args = parser.parse_args()

if args.dopt:
    dopt    = True
    donjets = False
else:
    dopt    = False
    donjets = True


########################################
# Global
########################################

datestr = strftime( '%b%d' )

baseDir = abspath( join( os.environ['CMSSW_BASE'], 'src' ) )
os.chdir( baseDir )

# Make a directory to contain job output; clean if necessary
jobDir = 'jobs{0}'.format( datestr )
if args.clean:
    if isdir( join( baseDir, jobDir ) ):
        shutil.rmtree( jobDir )
    sys.exit()


# Determine platform
if 'lxplus' in os.environ['HOSTNAME']:
    host = 'lxplus'
elif 't3ui' in os.environ['HOSTNAME']:
    host = 'psi'
else:
    host = 'unknown'


def main():


    ########################################
    # Setting paths and variables for pt
    ########################################

    if args.dopt:

        if not args.newbins:

            chapter( 'Setting paths and variables for pt (old binning)' )

            nBins = 7
            nCats = 3

            ptDir = join( baseDir, 'pt_moriond17_Mar01_oldBins' )
            wsDir = join( baseDir, 'workspaces{0}'.format(datestr) )
            if not isdir( wsDir ): os.makedirs(wsDir)

            DataCard            = join( ptDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_oldBins.txt' )
            regularizedDataCard = join( ptDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_oldBins_regularized.txt' )

            genDataCardRoot         = join( wsDir, basename(DataCard.replace( '.txt', '_genMuFit.root' )) )
            recoDataCardRoot        = join( wsDir, basename(DataCard.replace( '.txt', '_recoMuFit.root' )) )
            reguralizedDataCardRoot = join( wsDir, basename(regularizedDataCard.replace( '.txt', '.root' )) )

            genPostFitRoot         = genDataCardRoot.replace( '.root', '_postfit.root' )
            recoPostFitRoot        = recoDataCardRoot.replace( '.root', '_postfit.root' )
            reguralizedPostFitRoot = reguralizedDataCardRoot.replace( '.root', '_postfit.root' )

            print 'nBins                    = {0}'.format( nBins )
            print 'nCats                    = {0}'.format( nCats )
            print 'ptDir                    = {0}'.format( ptDir )
            print 'wsDir                    = {0}'.format( wsDir )
            print 'DataCard                 = {0}'.format( DataCard )
            print 'regularizedDataCard      = {0}'.format( regularizedDataCard )
            print 'genDataCardRoot          = {0}'.format( genDataCardRoot )
            print 'recoDataCardRoot         = {0}'.format( recoDataCardRoot )
            print 'reguralizedDataCardRoot  = {0}'.format( reguralizedDataCardRoot )
            print 'genPostFitRoot           = {0}'.format( genPostFitRoot )
            print 'recoPostFitRoot          = {0}'.format( recoPostFitRoot )
            print 'reguralizedPostFitRoot   = {0}'.format( reguralizedPostFitRoot )

            process(
                nBins,
                nCats,
                DataCard,
                regularizedDataCard,
                genDataCardRoot,
                recoDataCardRoot,
                reguralizedDataCardRoot,
                genPostFitRoot,
                recoPostFitRoot,
                reguralizedPostFitRoot,
                )


        else:

            chapter( 'Setting paths and variables for pt with extra bin' )

            nBins = 8
            nCats = 3

            ptDir = join( baseDir, 'pt_moriond17_Mar01_ExtraBin' )
            wsDir = join( baseDir, 'workspaces{0}'.format(datestr) )
            if not isdir( wsDir ): os.makedirs(wsDir)

            DataCard            = join( ptDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin.txt' )
            regularizedDataCard = join( ptDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_regularized.txt' )

            genDataCardRoot         = join( wsDir, basename(DataCard.replace( '.txt', '_genMuFit.root' )) )
            recoDataCardRoot        = join( wsDir, basename(DataCard.replace( '.txt', '_recoMuFit.root' )) )
            reguralizedDataCardRoot = join( wsDir, basename(regularizedDataCard.replace( '.txt', '.root' )) )

            genPostFitRoot         = genDataCardRoot.replace( '.root', '_postfit.root' )
            recoPostFitRoot        = recoDataCardRoot.replace( '.root', '_postfit.root' )
            reguralizedPostFitRoot = reguralizedDataCardRoot.replace( '.root', '_postfit.root' )

            print 'nBins                    = {0}'.format( nBins )
            print 'nCats                    = {0}'.format( nCats )
            print 'ptDir                    = {0}'.format( ptDir )
            print 'wsDir                    = {0}'.format( wsDir )
            print 'DataCard                 = {0}'.format( DataCard )
            print 'regularizedDataCard      = {0}'.format( regularizedDataCard )
            print 'genDataCardRoot          = {0}'.format( genDataCardRoot )
            print 'recoDataCardRoot         = {0}'.format( recoDataCardRoot )
            print 'reguralizedDataCardRoot  = {0}'.format( reguralizedDataCardRoot )
            print 'genPostFitRoot           = {0}'.format( genPostFitRoot )
            print 'recoPostFitRoot          = {0}'.format( recoPostFitRoot )
            print 'reguralizedPostFitRoot   = {0}'.format( reguralizedPostFitRoot )

            process(
                nBins,
                nCats,
                DataCard,
                regularizedDataCard,
                genDataCardRoot,
                recoDataCardRoot,
                reguralizedDataCardRoot,
                genPostFitRoot,
                recoPostFitRoot,
                reguralizedPostFitRoot,
                )

    ########################################
    # Setting paths and variables for pt
    ########################################

    if args.donjets:

        chapter( 'Setting paths and variables for pt' )

        nBins = 5
        nCats = 3

        nJetsDir = join( baseDir, 'nJets_moriond17_Mar01' )

        wsDir = join( baseDir, 'workspaces{0}'.format(datestr) )
        if not isdir( wsDir ): os.makedirs(wsDir)

        DataCard            = join( nJetsDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod.txt' )
        regularizedDataCard = join( nJetsDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_regularized.txt' )

        genDataCardRoot         = join( wsDir, basename(DataCard.replace( '.txt', '_genMuFit.root' )) )
        recoDataCardRoot        = join( wsDir, basename(DataCard.replace( '.txt', '_recoMuFit.root' )) )
        reguralizedDataCardRoot = join( wsDir, basename(regularizedDataCard.replace( '.txt', '.root' )) )

        genPostFitRoot         = genDataCardRoot.replace( '.root', '_postfit.root' )
        recoPostFitRoot        = recoDataCardRoot.replace( '.root', '_postfit.root' )
        reguralizedPostFitRoot = reguralizedDataCardRoot.replace( '.root', '_postfit.root' )

        print 'nBins                    = {0}'.format( nBins )
        print 'nCats                    = {0}'.format( nCats )
        print 'nJetsDir                    = {0}'.format( nJetsDir )
        print 'wsDir                    = {0}'.format( wsDir )
        print 'DataCard                 = {0}'.format( DataCard )
        print 'regularizedDataCard      = {0}'.format( regularizedDataCard )
        print 'genDataCardRoot          = {0}'.format( genDataCardRoot )
        print 'recoDataCardRoot         = {0}'.format( recoDataCardRoot )
        print 'reguralizedDataCardRoot  = {0}'.format( reguralizedDataCardRoot )
        print 'genPostFitRoot           = {0}'.format( genPostFitRoot )
        print 'recoPostFitRoot          = {0}'.format( recoPostFitRoot )
        print 'reguralizedPostFitRoot   = {0}'.format( reguralizedPostFitRoot )

        process(
            nBins,
            nCats,
            DataCard,
            regularizedDataCard,
            genDataCardRoot,
            recoDataCardRoot,
            reguralizedDataCardRoot,
            genPostFitRoot,
            recoPostFitRoot,
            reguralizedPostFitRoot,
            )




def process(
    nBins,
    nCats,
    DataCard,
    regularizedDataCard,
    genDataCardRoot,
    recoDataCardRoot,
    reguralizedDataCardRoot,
    genPostFitRoot,
    recoPostFitRoot,
    reguralizedPostFitRoot,
    ):


    ########################################
    # text2workspace and best fit 
    ########################################

    if args.doT2WS:
        chapter( 'Running text2workspace' )
        t2w(
            DataCard,
            genDataCardRoot,
            isPtFit=dopt,
            isGenFit=True,
            )
        t2w(
            DataCard,
            recoDataCardRoot,
            isPtFit=dopt,
            isGenFit=False,
            )
        # t2w(
        #     regularizedDataCard,
        #     reguralizedDataCardRoot,
        #     isPtFit=True,
        #     isGenFit=True,
        #     )

    if args.dobestfit:
        chapter( 'Making the best asimov fit' )
        bestfit(
            genDataCardRoot,
            genPostFitRoot,
            nBins,
            0,
            )

        bestfit(
            recoDataCardRoot,
            recoPostFitRoot,
            nBins,
            nCats,
            )
        # bestfit(
        #     regularizedDataCardRoot,
        #     regularizedPostFitRoot,
        #     nBins,
        #     0,
        #     )


    ########################################
    # Run jobs
    ########################################

    if args.dotoys:
        chapter( 'Generating and fitting toys' )

        if not args.startseed:
            seed  = 100
        else:
            seed = args.startseed


        for iToy in xrange(args.nToys):

            try:

                runJob(
                    seed,
                    genDataCardRoot,
                    recoDataCardRoot,
                    reguralizedDataCardRoot,
                    genPostFitRoot,
                    recoPostFitRoot,
                    reguralizedPostFitRoot,
                    nBins,
                    nCats,
                    )

            finally:
                os.chdir( baseDir )

            seed += 1





def runJob(
    seed,
    genDataCardRoot,
    recoDataCardRoot,
    reguralizedDataCardRoot,
    genPostFitRoot,
    recoPostFitRoot,
    reguralizedPostFitRoot,
    nBins,
    nCats,
    ):

    if not args.test and not args.printCommands:
        if not isdir( jobDir ): os.makedirs( jobDir )
        tmpdir = tempfile.mkdtemp( prefix='tmptoy_seed{0}_'.format(seed), dir=jobDir )
        tmpdir = abspath( tmpdir )

        print '\nEntering {0}'.format( tmpdir )
        os.chdir( tmpdir )

        # # Get freezeVars; requires going back to baseDir temporarily
        # os.chdir( baseDir )
        # freezeVars = getCommaSeparatedListOfVariablesToFreeze( genPostFitRoot, freezeAll=args.freezeAllPdfVars )
        # os.chdir( tmpdir )


    # Make the toy command
    toycmd = ''
    toycmd += " combine  "
    toycmd += " {0} ".format( recoPostFitRoot )
    toycmd += " --seed {0} ".format(seed)
    toycmd += " -M GenerateOnly  "
    toycmd += " --snapshotName MultiDimFit "
    toycmd += " -n generatedToy "

    toycmd += " --saveToys  "
    toycmd += " -t 1 "
    toycmd += " --toysFrequentist --bypassFrequentistFit  "

    toycmd += " -m 125   "
    toycmd += " --setPhysicsModelParameters {0} ".format( getPhysicsModelParameterString(nBins,nCats) )

    toysFile = abspath( 'higgsCombinegeneratedToy.GenerateOnly.mH125.{0}.root'.format( seed ) )


    # Make the bestfit command
    bestfitcmd = ''
    bestfitcmd += "combine  "

    # bestfitcmd += " {0} ".format( recoDataCardRoot )
    bestfitcmd += " {0} ".format( recoPostFitRoot )

    # bestfitcmd += " --seed {0} ".format(seed)

    bestfitcmd += " -M MultiDimFit  "
    bestfitcmd += " -n bestfitToy "

    # bestfitcmd += " -t 1 "
    # bestfitcmd += " --toysFile {0} ".format( toysFile )
    bestfitcmd += " --dataset {0}:toys/toy_1 ".format( toysFile )

    bestfitcmd += " -m 125  "
    # bestfitcmd += " --setPhysicsModelParameters {0} ".format( getPhysicsModelParameterString(nBins,nCats) )

    bestfitcmd += " --floatOtherPOIs=1  "
    bestfitcmd += " --saveWorkspace  "
    bestfitcmd += " --minimizerStrategy 2 "

    # bestfitToyFile = 'higgsCombinebestfitToy.MultiDimFit.mH125.{0}.root'.format( seed )
    bestfitToyFile = 'higgsCombinebestfitToy.MultiDimFit.mH125.root'


    freezepdfvarscmd =  'python '
    freezepdfvarscmd += ' {0} '.format( join( baseDir, 'getCommaSeparatedListOfVariablesToFreeze.py' ) )
    freezepdfvarscmd += ' {0} '.format( bestfitToyFile )

    covmatcmd = ''
    covmatcmd += " combine "
    covmatcmd += " {0} ".format( bestfitToyFile )
    covmatcmd += " -M MultiDimFit "
    covmatcmd += " -n covmat "
    covmatcmd += " -m 125 "
    covmatcmd += " --saveWorkspace  "

    covmatcmd += " --dataset {0}:toys/toy_1 ".format( toysFile )
    # covmatcmd += " -t 1 "
    # covmatcmd += " --toysFrequentist --bypassFrequentistFit "

    covmatcmd += " --algo none "
    covmatcmd += " --snapshotName MultiDimFit "
    covmatcmd += " --setPhysicsModelParameters {0} ".format( getPhysicsModelParameterString(nBins,nCats) )
    covmatcmd += " --freezeNuisances \"$({0})\" ".format( freezepdfvarscmd )
    covmatcmd += " -v 2"


    if args.printCommands:
        print '\nToy generation command:'
        print toycmd
        print '\nBest fit command:'
        print bestfitcmd
        print '\nCorrelation matrix command:'
        print covmatcmd
        return


    shfile = 'runBestFitOnToy.sh'

    with open( shfile, 'w' ) as shfp:

        w = lambda text: shfp.write( text + ' \n' )


        if host == 'psi' or host == 'unknown':

            w( '#$ -q all.q' )
            # w( '#$ -o {0}'.format( tmpdir ) )
            # w( '#$ -e {0}'.format( tmpdir ) )
            w( 'cd {0}'.format(tmpdir) )
            w( '#$ -cwd' )
            w( 'OUTFILES=""' )
            w( 'JOBDIR=/scratch/`whoami`/job_${JOB_ID}' )
            w( 'mkdir -p $JOBDIR' )
            w( 'echo "JOBDIR = $JOBDIR"' )
            w( 'echo "HOSTNAME = $HOSTNAME"' )

            w( 'source $VO_CMS_SW_DIR/cmsset_default.sh' )
            w( 'cd {0}'.format(baseDir) )
            w( 'eval `scramv1 runtime -sh`' )
            w( 'cd $JOBDIR' )


            w( toycmd )
            w( bestfitcmd )
            w( covmatcmd )

            w( 'cp {0} {1}'.format( bestfitToyFile, join( tmpdir, bestfitToyFile ) ) )


        elif host == 'lxplus':

            w( 'cd {0}'.format( baseDir ) )
            w( 'eval `scramv1 runtime -sh`' )
            w( 'cd {0}'.format( tmpdir ) )
            w( toycmd )
            w( bestfitcmd )
            w( covmatcmd )


    if host == 'lxplus':
        os.system( 'chmod 777 {0}'.format(shfile) )


    if not args.local:
        if host == 'psi':
            if args.moreMemory:
                executeCommand( 'qsub -l h_vmem=4g {0}'.format( shfile ) )
            else:
                executeCommand( 'qsub {0}'.format( shfile ) )
        elif host == 'lxplus':
            executeCommand( 'bsub -q 1nh {0}'.format( shfile ) )
    else:
        executeCommand( toycmd )
        executeCommand( bestfitcmd )
        executeCommand( covmatcmd )






def bestfit(
    dataCardRoot,
    postFitRoot,
    nBins,
    nCats,
    docovmat=True,
    ):


    name, extension = splitext(basename(dataCardRoot))
    name = 'asimov_' + name

    bestfitcmd = ''
    bestfitcmd += " combine  "
    bestfitcmd += " {0} ".format(dataCardRoot)
    bestfitcmd += " -M MultiDimFit  "
    bestfitcmd += " -t -1 "
    bestfitcmd += " --saveWorkspace  "
    bestfitcmd += " -n {0} ".format( name )
    bestfitcmd += " --setPhysicsModelParameters {0} ".format( getPhysicsModelParameterString(nBins,nCats) )
    bestfitcmd += " -m 125  "
    bestfitcmd += " --minimizerStrategy 2 "
    bestfitcmd += " ; mv higgsCombine{0}.MultiDimFit.mH125.root {1} ".format( name, postFitRoot )

    executeCommand( bestfitcmd )
    

    if docovmat:

        freezeVars = getCommaSeparatedListOfVariablesToFreeze( postFitRoot, freezeAll=False )

        covmatRoot = dataCardRoot.replace( '.root', '_covmat.root' )
        covmatWSRoot = dataCardRoot.replace( '.root', '_covmatWS.root' )

        covmatcmd = ''
        covmatcmd += " combine "
        covmatcmd += " {0} ".format( postFitRoot )
        covmatcmd += " -M MultiDimFit "
        covmatcmd += " -n covmat "
        covmatcmd += " -m 125 "
        covmatcmd += " --saveWorkspace  "


        covmatcmd += " -t -1 "
        covmatcmd += " --toysFrequentist --bypassFrequentistFit "

        covmatcmd += " --algo none "
        covmatcmd += " --snapshotName MultiDimFit "
        covmatcmd += " --setPhysicsModelParameters {0} ".format( getPhysicsModelParameterString(nBins,nCats) )
        covmatcmd += " --freezeNuisances {0} ".format( freezeVars )
        covmatcmd += " -v 2"
        covmatcmd += " ; mv COVMATISHERE.root {0} ".format( covmatRoot )
        covmatcmd += " ; mv higgsCombinecovmat.MultiDimFit.mH125.root {0} ".format( covmatWSRoot )

        executeCommand( covmatcmd )


    printWSVars( postFitRoot, verbose=True )
    printWSVars( covmatWSRoot, verbose=True )



def t2w( DC, DCroot, isPtFit, isGenFit ):

    cmd = ''

    if isPtFit and isGenFit:
        cmd += "text2workspace.py "
        cmd += "{0} ".format( DC )
        cmd += "-o {0} ".format( DCroot )
        cmd += "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "
        cmd += "--PO verbose "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_0p0to15p0:r0[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_15p0to30p0:r1[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_30p0to45p0:r2[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_45p0to85p0:r3[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_85p0to125p0:r4[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_125p0to200p0:r5[1,0,2]' "
        cmd += "--PO 'map=.*/InsideAcceptance_genPt_200p0to10000p0:r6[1,0,2]' "

    elif isPtFit and not isGenFit:
        cmd += "text2workspace.py  "
        cmd += "{0} ".format( DC )
        cmd += "-o {0} ".format( DCroot )
        cmd += "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  "
        cmd += "--PO verbose  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_0p0to15p0.*/InsideAcceptance.*:r0_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_15p0to30p0.*/InsideAcceptance.*:r1_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_30p0to45p0.*/InsideAcceptance.*:r2_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_45p0to85p0.*/InsideAcceptance.*:r3_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_85p0to125p0.*/InsideAcceptance.*:r4_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_125p0to200p0.*/InsideAcceptance.*:r5_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoPt_200p0to10000p0.*/InsideAcceptance.*:r6_cat2[1,0,2]' "

    elif (not isPtFit ) and isGenFit:
        cmd += "text2workspace.py  "
        cmd += "{0} ".format( DC )
        cmd += "-o {0} ".format( DCroot )
        cmd += "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  "
        cmd += "--PO verbose  "
        cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_m0p5to0p5:r0[1,0,2]'  "
        cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_0p5to1p5:r1[1,0,2]'  "
        cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_1p5to2p5:r2[1,0,2]'  "
        cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_2p5to3p5:r3[1,0,2]'  "
        cmd += "--PO 'map=.*/InsideAcceptance_genNjets2p5_3p5to100p0:r4[1,0,2]' "

    elif (not isPtFit) and (not isGenFit):
        cmd += "text2workspace.py  "
        cmd += "{0} ".format( DC )
        cmd += "-o {0} ".format( DCroot )
        cmd += "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  "
        cmd += "--PO verbose  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_0_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat0[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_1_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat1[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoNjets2p5_m0p5to0p5.*/InsideAcceptance.*:r0_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoNjets2p5_0p5to1p5.*/InsideAcceptance.*:r1_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoNjets2p5_1p5to2p5.*/InsideAcceptance.*:r2_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoNjets2p5_2p5to3p5.*/InsideAcceptance.*:r3_cat2[1,0,2]'  "
        cmd += "--PO 'map=.*SigmaMpTTag_2_recoNjets2p5_3p5to100p0.*/InsideAcceptance.*:r4_cat2[1,0,2]' "

    executeCommand( cmd )



def getPhysicsModelParameterString( nBins, nCats=0 ):

    ret = []

    if nCats == 0:
        for iBin in xrange(nBins):
            ret.append( 'r{0}=1'.format(iBin) )
    else:
        for iCat in xrange(nCats):
            for iBin in xrange(nBins):
                ret.append( 'r{0}_cat{1}=1'.format(iBin,iCat) )

    return ','.join(ret)


def chapter( text ):
    print
    print '='*70
    print '=== {0}'.format( text )
    print '='*70
    print


def executeCommand( cmd, verbose=True ):
    if not args.test:
        if verbose: print '\nEXECUTING: {0}'.format(cmd)
        os.system( cmd )
    else:
        print '\nTESTMODE: {0}'.format(cmd)



def printWSVars( wsFile, verbose=True ):
    
    wsFp = ROOT.TFile.Open( wsFile )
    w = wsFp.Get('w')

    try:
        w.loadSnapshot('MultiDimFit')
    except:
        print 'Snapshot loading failed'
        return


    variables = []

    varArgSet = ROOT.RooArgList( w.allVars() )
    for i in xrange(varArgSet.getSize()):
        variables.append(( varArgSet[i].GetName(), varArgSet[i].getVal(), varArgSet[i].getError() ))

        # variable = varArgSet[i]
        # name = variable.GetName()
        # match = re.search( r'r\d', name )
        # if match:
        #     print '{0:20} = {1:20}     +-  {2}'.format( name, variable.getVal(), variable.getError() )


    catArgSet = ROOT.RooArgList( w.allCats() )
    for i in xrange(catArgSet.getSize()):
        variables.append(( catArgSet[i].GetName(), catArgSet[i].getIndex(), None ))


    variables.sort()
    variables.sort( key=lambda i: -2 * i[0].startswith('r') + -1*i[0].startswith('pdfindex') )

    outFile = wsFile.replace( '.root', '_variablePrintout.txt' )
    with open( outFile, 'w' ) as outFp:
        for name, val, err in variables:

            if not err:
                printstr = '{0} = {1}'.format( name, val )
            else:
                printstr = '{0} = {1}   +- {2}'.format( name, val, err )

            outFp.write( printstr + '\n' )
            if verbose: print printstr




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()