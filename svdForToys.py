#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, sys, re, json, glob, shutil
from os.path import *
from math import exp, sqrt, log, pi

import ROOT

from time import strftime
datestr = strftime( '%b%d' )

from reg2 import doSVDMathematics

# doSVDMathematics(
#     name,
#     nBins,
#     nCats,
#     Binning,
#     genMuFitFile,
#     recoMuFitFile,
#     responseMatrixTxtFile,
#     recoCovMatRootFile,
#     genCovMatRootFile,
#     plotdir,
#     seed=123,
#     makePlots=True,
#     )


########################################
# Main
########################################

def main():


    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument( '--string', type=str, default='default', help='default string' )
    parser.add_argument( '--ispt', action='store_true', help='boolean')
    parser.add_argument( '--newbins', action='store_true', help='boolean')
    # parser.add_argument( 'jobdirs', metavar='N', type=str, nargs='+', help='list of strings' )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument( '--jobdirs', metavar='N', type=str, nargs='+', help='list of strings' )
    group.add_argument( '--plotk', action='store_true', help='boolean')

    parser.add_argument( '--verbose', action='store_true', help='boolean')
    args = parser.parse_args()


    baseDir = abspath( join( os.environ['CMSSW_BASE'], 'src' ))
    wsDir =  join( baseDir, 'workspacesMar01' )


    if not args.ispt:

        genMuFitFile          = join( wsDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_genMuFit_postfit.root' )
        responseMatrixTxtFile = 'signalEv_transferMatrix_genNjets2p5_vs_recoNjets2p5_all_normToUnity_v2.txt'
        genCovMatRootFile     = join( wsDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_genMuFit_covmat.root' )

        overArchingPlotdir = 'plots_{0}_nJets_jobs'.format( datestr )

        nBins   = 5
        nCats   = 3
        Binning = [ -0.5, 0.5, 1.5, 2.5, 3.5, 5.5 ]

    else:

        if not args.newbins:
            genMuFitFile          = join( wsDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_oldBins_genMuFit_postfit.root' )
            responseMatrixTxtFile = 'signalEv_transferMatrix_genPt_vs_recoPt_all_normToUnity_v2.txt'
            genCovMatRootFile     = join( wsDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_oldBins_genMuFit_covmat.root' )

            overArchingPlotdir = 'plots_{0}_pt_oldbins_jobs'.format( datestr )

            nBins   = 7
            nCats   = 3
            Binning = [ 0., 15., 30., 45., 85., 125., 200., 300. ]

        else:
            genMuFitFile          = join( wsDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_genMuFit_postfit.root' )
            responseMatrixTxtFile = 'signalEv_transferMatrix_genPt_vs_recoPt_all_normToUnity_v2.txt'
            genCovMatRootFile     = join( wsDir, 'Datacard_13TeV_differential_pT_moriond17_reminiaod_extrabin_genMuFit_covmat.root' )

            overArchingPlotdir = 'plots_{0}_pt_extrabin_jobs'.format( datestr )

            nBins   = 8
            nCats   = 3
            Binning = [ 0., 15., 30., 45., 85., 125., 200., 350., 500. ]
        






    if not args.plotk:

        for jobdir in args.jobdirs:

            print '\nProcessing {0}'.format( jobdir )

            jobdir = abspath( jobdir )
            seed = re.search( r'seed(\d+)_', basename(jobdir) ).group(1)

            recoMuFitFile         = join( jobdir, 'higgsCombinebestfitToy.MultiDimFit.mH125.root' )
            recoCovMatRootFile    = join( jobdir, 'COVMATISHERE.root' )

            plotdir = '{0}/job_seed{1}'.format( overArchingPlotdir, seed )

            regResult = doSVDMathematics(
                'nJets',
                nBins,
                nCats,
                Binning,
                genMuFitFile,
                recoMuFitFile,
                responseMatrixTxtFile,
                recoCovMatRootFile,
                genCovMatRootFile,
                plotdir,
                forceVerbosity=args.verbose,
                )


            with open( join( plotdir, 'output.json' ), 'wb' ) as outputFp:

                dumpDict = {

                    'k' : regResult.k,
                    'S' : regResult.S.tolist(),
                    'correlationMat' : regResult.correlationMat.tolist(),
                    'covMat' : regResult.covMat.tolist(),
                    'd' : regResult.d.tolist(),

                    }

                json.dump( dumpDict, outputFp, sort_keys=True, indent=4)

    else:

        ks = []
        taus = []

        jsonFiles = glob.glob( '{0}/*/output.json'.format( overArchingPlotdir ) )
        for jsonFile in jsonFiles:

            with open( jsonFile, 'r' ) as jsonFp:
                D = json.load( jsonFp )

                k = D['k']

                if not k == -1:
                    ks.append(k)

                    if k > nBins-1:
                        tau = sqrt( D['S'][nBins-1][nBins-1] )
                    else:
                        tau = sqrt( D['S'][ int(round(k)) ][ int(round(k)) ] )

                    taus.append(tau)


        print ks
        print taus




        plotdirs = glob.glob( '{0}/job_seed*'.format( overArchingPlotdir ) )
        ddir = '{0}/allds'.format(overArchingPlotdir)
        if not isdir(ddir): os.makedirs(ddir)

        for plotdir in plotdirs:

            seed = re.search( r'seed(\d+)', basename(plotdir) ).group(1)

            src = join( plotdir, 'd.pdf' )

            dst = join( ddir, 'd_seed{0}.pdf'.format(seed) )

            shutil.copyfile( src, dst )


        ########################################
        # average k plot
        ########################################

        ck = ROOT.TCanvas( 'ck', 'ck', 800, 600 )

        k_TH1F = ROOT.TH1F( 'k', 'k', 20, 0., nBins+1 )
        for k in ks:
            k_TH1F.Fill( k )

        k_TH1F.Draw()


        # Checking at which bin we cross certain threshold
        ck.Update()
        ck.cd()
        yMax = ROOT.gPad.GetUymax()
        yMin = ROOT.gPad.GetUymin()

        k_average = sum(ks)/len(ks)
        l = ROOT.TLine( k_average, yMin, k_average, yMax )
        ROOT.SetOwnership( l, False )
        l.SetLineStyle(2)
        l.Draw()

        ck.SaveAs( '{0}/averagek.pdf'.format(overArchingPlotdir) )



        ########################################
        # average tau plot
        ########################################



        ctau = ROOT.TCanvas( 'ctau', 'ctau', 800, 600 )

        tau_TH1F = ROOT.TH1F( 'tau', 'tau', 20, 0., 1. )
        for tau in taus:
            tau_TH1F.Fill( tau )

        tau_TH1F.Draw()


        # Checking at which bin we cross certain threshold
        ctau.Update()
        ctau.cd()
        yMax = ROOT.gPad.GetUymax()
        yMin = ROOT.gPad.GetUymin()

        tau_average = sum(taus)/len(taus)
        l = ROOT.TLine( tau_average, yMin, tau_average, yMax )
        ROOT.SetOwnership( l, False )
        l.SetLineStyle(2)
        l.Draw()

        ctau.SaveAs( '{0}/averagetau.pdf'.format(overArchingPlotdir) )










########################################
# End of Main
########################################
if __name__ == "__main__":
    main()