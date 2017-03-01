#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
from os.path import *
import argparse
from array import array
import random
from math import exp, log, sqrt, pi

import numpy
numpy.set_printoptions( precision=5, linewidth=100 )

from time import strftime
datestr = strftime( '%b%d' )

import ROOT

from makeRegPlots import RegResult

ROOT.gROOT.SetBatch(True)
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
# ROOT.gSystem.Load('/mnt/t3nfs01/data01/shome/tklijnsm/Combination/RooUnfold-1.1.1/libRooUnfold')

Verbosity = 0


########################################
# Main
########################################

def main():

    baseDir = abspath( join( os.environ['CMSSW_BASE'], 'src' ))


    # do_pT = False
    do_pT = True

    do_nJets = False
    # do_nJets = True

    do_regularized = True


    global Verbosity
    Verbosity = 1


    if do_pT:

        dcDir =  join( baseDir, 'pt_moriond17' )

        genMuFitDir            = abspath(join( baseDir, 'results_Feb23_pt_moriond17' ))
        genMuFitFile          = join( genMuFitDir,  'Datacard_13TeV_differential_pT_moriond17_postfit.root' )

        recoMuFitDir           = abspath(join( baseDir, 'results_Feb23_pt_moriond17_recoMuFit' ))
        recoMuFitFile         = join( recoMuFitDir, 'Datacard_13TeV_differential_pT_moriond17_postfit_recoMuFit.root' )

        responseMatrixTxtFile = 'eff_times_acceptance.txt'
        recoCovMatRootFile    = join( recoMuFitDir, 'COVMATISHERE.root' )
        genCovMatRootFile     = join( genMuFitDir,  'COVMATISHERE.root' )

        nBins   = 7
        nCats   = 3
        Binning = [ 0., 15., 30., 45., 85., 125., 200., 300. ]

        plotdir = 'plots_{0}_pT'.format(datestr)

        pT_regResult = doSVDMathematics(
            'pt',
            nBins,
            nCats,
            Binning,
            genMuFitFile,
            recoMuFitFile,
            responseMatrixTxtFile,
            recoCovMatRootFile,
            genCovMatRootFile,
            plotdir,
            )


        if do_regularized:

            # ======================================
            # Get observed genMus

            nBins   = 7
            nCats   = 3
            Binning = [ 0., 15., 30., 45., 85., 125., 200., 300. ]

            regularizedGenMuFitDir = abspath(join( baseDir, 'results_Feb27_pt_moriond17_regularized' ))

            regularizedGenMuFitFile   = join( regularizedGenMuFitDir, 'Datacard_13TeV_differential_pT_moriond17_regularization_postfit.root' )
            regularizedGenMuFitRootFp = ROOT.TFile.Open( regularizedGenMuFitFile )
            regularizedGenMuFitW      = regularizedGenMuFitRootFp.Get('w')
            regularizedGenMuFitW.loadSnapshot('MultiDimFit')

            regularizedGenMuObserved      = numpy.zeros( (nBins,1) )
            regularizedGenMuObservedError = numpy.zeros( (nBins,1) )
            for iBin in xrange(nBins):
                r = regularizedGenMuFitW.var( 'r{0}'.format(iBin) )
                mu    = r.getVal()
                error = r.getError()

                # Simulate a toy
                # mu = Smear( mu, error )

                # Add to lists
                regularizedGenMuObserved[iBin] =  mu
                regularizedGenMuObservedError[iBin] =  error


            oneOverSqrtTau = regularizedGenMuFitW.var('OneOverSqrtTau').getVal()
            tau = 1./(oneOverSqrtTau**2)
            if Verbosity>0: print '1/sqrt(tau) = {0}'.format( oneOverSqrtTau )

            if Verbosity>0: print '\nPrinting regularizedGenMuObserved\n', regularizedGenMuObserved
            if Verbosity>0: print '\nPrinting regularizedGenMuObservedError\n', regularizedGenMuObservedError

            regularizedGenMuFitRootFp.Close()


            pT_regResult.compareRegularizedSpectrum( regularizedGenMuObserved, regularizedGenMuObservedError )

            if Verbosity>0: print '\nPrinting pT_regResult.genMuObserved\n', pT_regResult.genMuObserved
            if Verbosity>0: print '\nPrinting pT_regResult.genMuObservedError\n', pT_regResult.genMuObservedError


            # Draw correlation matrix
            regularizedCovMatFile = join( regularizedGenMuFitDir, 'COVMATISHERE.root' )
            regularizedCorrelationMat, regularizedCorrelationMat_byCat, regularizedCovMat = readCorrMatrix( regularizedCovMatFile, nBins, 1, regularizedGenMuObservedError )
            pT_regResult.drawMatrixInRoot( regularizedCorrelationMat, 'corrMat_GenMuRegularized', title='p_{{T}} correlation matrix after regularization (in gen bins, #tau={0})'.format(tau) )




    if do_nJets:

        dcDir =  join( baseDir, 'nJets_moriond17' )

        genMuFitDir  = abspath( join( baseDir, 'results_Feb23_nJets_moriond17' ))
        recoMuFitDir = abspath( join( baseDir, 'results_Feb23_nJets_moriond17_recoMuFit' ))

        responseMatrixTxtFile = 'eff_times_acceptance_nJets.txt'
        recoCovMatRootFile    = join( recoMuFitDir, 'COVMATISHERE.root' )
        genCovMatRootFile     = join( genMuFitDir,  'COVMATISHERE.root' )
        genMuFitFile          = join( genMuFitDir,  'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_postfit.root' )
        recoMuFitFile         = join( recoMuFitDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_postfit_recoMuFit.root' )

        nBins   = 5
        nCats   = 3
        Binning = [ -0.5, 0.5, 1.5, 2.5, 3.5, 5.5 ]

        plotdir = 'plots_{0}_nJets'.format(datestr)

        nJets_regResult = doSVDMathematics(
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
            )


        if do_regularized:

            nBins   = 5
            nCats   = 3
            Binning = [ -0.5, 0.5, 1.5, 2.5, 3.5, 5.5 ]

            # regularizedGenMuFitDir    = abspath( join( baseDir, 'results_Feb23_pt_moriond17_regularized' ))
            regularizedGenMuFitDir = abspath(join( baseDir, 'results_Feb24_nJets_moriond17_regularized' ))

            regularizedGenMuFitFile   = join( regularizedGenMuFitDir, 'Datacard_13TeV_differential_Njets_moriond17_skipAndDebug_reminiaod_regularized_postfit.root' )
            regularizedGenMuFitRootFp = ROOT.TFile.Open( regularizedGenMuFitFile )
            regularizedGenMuFitW      = regularizedGenMuFitRootFp.Get('w')
            regularizedGenMuFitW.loadSnapshot('MultiDimFit')

            regularizedGenMuObserved      = numpy.zeros( (nBins,1) )
            regularizedGenMuObservedError = numpy.zeros( (nBins,1) )
            for iBin in xrange(nBins):
                r = regularizedGenMuFitW.var( 'r{0}'.format(iBin) )
                mu    = r.getVal()
                error = r.getError()

                # Simulate a toy
                # mu = Smear( mu, error )

                # Add to lists
                regularizedGenMuObserved[iBin] =  mu
                regularizedGenMuObservedError[iBin] =  error

            if Verbosity>0: print '\nPrinting regularizedGenMuObserved\n', regularizedGenMuObserved
            if Verbosity>0: print '\nPrinting regularizedGenMuObservedError\n', regularizedGenMuObservedError

            oneOverSqrtTau = regularizedGenMuFitW.var('OneOverSqrtTau').getVal()
            tau = 1./(oneOverSqrtTau**2)
            if Verbosity>0: print '1/sqrt(tau) = {0}'.format( oneOverSqrtTau )

            regularizedGenMuFitRootFp.Close()


            nJets_regResult.compareRegularizedSpectrum( regularizedGenMuObserved, regularizedGenMuObservedError )

            if Verbosity>0: print '\nPrinting nJets_regResult.genMuObserved\n', nJets_regResult.genMuObserved
            if Verbosity>0: print '\nPrinting nJets_regResult.genMuObservedError\n', nJets_regResult.genMuObservedError


            # Draw correlation matrix
            regularizedCovMatFile = join( regularizedGenMuFitDir, 'COVMATISHERE.root' )
            regularizedCorrelationMat, regularizedCorrelationMat_byCat, regularizedCovMat = readCorrMatrix( regularizedCovMatFile, nBins, 1, regularizedGenMuObservedError )
            nJets_regResult.drawMatrixInRoot( regularizedCorrelationMat, 'corrMat_GenMuRegularized', title='nJets correlation matrix after regularization (in gen bins, #tau={0})'.format(tau) )






################################################################################
# SVD mathematics
################################################################################


def doSVDMathematics(
    name,
    nBins,
    nCats,
    Binning,
    genMuFitFile,
    recoMuFitFile,
    responseMatrixTxtFile,
    recoCovMatRootFile,
    genCovMatRootFile,
    plotdir,
    seed=123,
    makePlots=True,
    ):

    # baseDir = abspath( join( os.environ['CMSSW_BASE'], 'src' ))
    # dcDir   =  join( baseDir, 'pt_moriond17' )

    # genMuFitDir  = abspath( join( baseDir, 'results_Feb15_asimov' ))
    # recoMuFitDir = abspath( join( baseDir, 'results_Feb20_asimovPerReco' ))

    # nBins   = 7
    # nCats   = 3
    # Binning = [ 0., 15., 30., 45., 85., 125., 200., 300. ]


    ########################################
    # Read root files
    ########################################

    chapter( 'Input from root files' )

    # ======================================
    # Get observed genMus

    # genMuFitFile   = join( genMuFitDir, 'Datacard_13TeV_differential_pT_moriond17_post_asmov_fit.root' )
    genMuFitRootFp = ROOT.TFile.Open( genMuFitFile )
    genMuFitW      = genMuFitRootFp.Get('w')
    genMuFitW.loadSnapshot('MultiDimFit')

    random.seed(123)
    Smear = lambda x, err: x + random.gauss( 0, err )

    genMuObserved      = numpy.zeros( (nBins,1) )
    genMuObservedError = numpy.zeros( (nBins,1) )
    for iBin in xrange(nBins):
        r = genMuFitW.var( 'r{0}'.format(iBin) )
        mu    = r.getVal()
        error = r.getError()

        # Simulate a toy
        # mu = Smear( mu, error )

        # Add to lists
        genMuObserved[iBin] =  mu
        genMuObservedError[iBin] =  error

    if Verbosity>0:     print '\nPrinting genMuObserved\n', genMuObserved
    if Verbosity>0: print '\nPrinting genMuObservedError\n', genMuObservedError

    genMuFitRootFp.Close()


    # ======================================
    # Get observed recoMus

    # recoMuFitFile   = join( recoMuFitDir, 'Datacard_13TeV_differential_pT_moriond17_post_asmov_fit_perReco.root' )
    recoMuFitRootFp = ROOT.TFile.Open( recoMuFitFile )
    recoMuFitW      = recoMuFitRootFp.Get('w')
    recoMuFitW.loadSnapshot('MultiDimFit')

    random.seed(123) # Default seed

    random.seed(seed)
    Smear = lambda x, err: x + random.gauss( 0, err )

    recoMuObserved      = numpy.zeros( ( nCats, nBins ) )
    recoMuObservedError = numpy.zeros( ( nCats, nBins ) )
    for iBin in xrange(nBins):
        for iCat in xrange(nCats):
            r = recoMuFitW.var( 'r{0}_cat{1}'.format( iBin, iCat ) )
            mu    = r.getVal()
            error = r.getError()

            # Simulate a toy
            mu = Smear( mu, error )

            # Add to lists
            recoMuObserved[iCat][iBin]      = mu
            recoMuObservedError[iCat][iBin] = error


    recoMuObserved = recoMuObserved.reshape(nBins*nCats,1)
    recoMuObservedError = recoMuObservedError.reshape(nBins*nCats,1)

    if Verbosity>0: print '\nPrinting recoMuObserved\n', recoMuObserved
    if Verbosity>0: print '\nPrinting recoMuObservedError\n', recoMuObservedError

    recoMuFitRootFp.Close()


    # ======================================
    # Get response matrix (A)

    # eff = 1./(nCats) * numpy.array([ numpy.identity(nBins).ravel() for iCat in xrange(nCats) ]).reshape( nCats*nBins, nBins )
    # if Verbosity>0: print '\nPrinting eff\n', eff

    # acceptance = 1./nBins * numpy.ones((nBins,1))
    # if Verbosity>0: print '\nPrinting acceptance\n', acceptance


    eff_times_acceptance = numpy.transpose(numpy.genfromtxt( responseMatrixTxtFile ))

    # Last row is noTag
    eff_times_acceptance = eff_times_acceptance[:-1]

    numpy.set_printoptions( precision=2 )
    if Verbosity>0: print '\nPrinting eff_times_acceptance\n', eff_times_acceptance
    numpy.set_printoptions( precision=5 )


    ########################################
    # Create A*x (smeared fitted genMus) and b (separately fitted recoMus)
    ########################################

    chapter( 'SVD input' )

    A = eff_times_acceptance
    x = genMuObserved

    summed_eff_times_acceptance = numpy.sum( eff_times_acceptance, axis=1 )
    if Verbosity>0: print '\nPrinting summed_eff_times_acceptance\n', summed_eff_times_acceptance

    # element-wise product keeps failing, doing manually:
    b = numpy.zeros((nBins*nCats,1))
    for iCat in xrange(nCats):
        for iBin in xrange(nBins):
            iBinReco =  iCat*nBins + iBin 
            b[iBinReco] = recoMuObserved[iBinReco] * summed_eff_times_acceptance[iBinReco]

    if Verbosity>0: print '\nPrinting A\n', A
    if Verbosity>0: print '\nPrinting x\n', x
    if Verbosity>0: print '\nPrinting b\n', b
    
    Ax = A.dot(x)
    if Verbosity>0: print '\nPrinting Ax\n', Ax



    ########################################
    # Covariance matrix
    ########################################

    chapter( 'Decomposing covariance matrix' )
    numpy.set_printoptions( precision=2 )

    correlationMat, correlationMat_byCat, covMat = readCorrMatrix( recoCovMatRootFile, nBins, nCats, recoMuObservedError )
    genCorrelationMat, genCorrelationMat_byCat, genCovMat = readCorrMatrix( genCovMatRootFile, nBins, 1, genMuObservedError )

    if Verbosity>0: print '\nPrinting correlationMat\n', correlationMat

    # # covmatRootFile = join( recoMuFitDir, 'COVMATISHERE.root' )
    # covmatRootFp   = ROOT.TFile.Open( covmatRootFile )
    # fit            = covmatRootFp.Get('fit')

    # covMat               = numpy.zeros( (nBins*nCats, nBins*nCats) )
    # correlationMat       = numpy.zeros( (nBins*nCats, nBins*nCats) )
    # correlationMat_byCat = numpy.zeros( (nBins*nCats, nBins*nCats) )

    # for iBin1 in xrange(nBins):
    #     for iCat1 in xrange(nCats):
    #         for iBin2 in xrange(nBins):
    #             for iCat2 in xrange(nCats):

    #                 iBinReco1 = iCat1*nBins+iBin1
    #                 iBinReco2 = iCat2*nBins+iBin2

    #                 correlationMat[iBinReco1][iBinReco2] = fit.correlation( 'r{0}_cat{1}'.format(iBin1,iCat1), 'r{0}_cat{1}'.format(iBin2,iCat2) )

    #                 covMat[iBinReco1][iBinReco2] = (
    #                     correlationMat[iBinReco1][iBinReco2]
    #                     * recoMuObservedError[iBinReco1]
    #                     * recoMuObservedError[iBinReco2]
    #                     )

    #                 correlationMat_byCat[iBin1*nCats+iCat1][iBin2*nCats+iCat2] = fit.correlation( 'r{0}_cat{1}'.format(iBin1,iCat1), 'r{0}_cat{1}'.format(iBin2,iCat2) )

    # # # Test with some migrations
    # # smearing = 0.3
    # # correlationMat = [ [ 0 for i in xrange(nBins) ] for j in xrange(nBins) ]
    # # correlationMat[0]  = [ 1.-smearing, smearing ] + [ 0. for i in xrange(nBins-2) ]
    # # correlationMat[-1] = [ 0. for i in xrange(nBins-2) ] + [ smearing, 1.-smearing ]
    # # for i in xrange(1,nBins-1):
    # #     correlationMat[i] = [ 0. for j in xrange(0,i-1) ] + [ smearing, 1.-2*smearing, smearing ] + [ 0. for j in xrange(i+2,nBins) ]

    # if Verbosity>0: print '\nPrinting correlationMat\n', correlationMat
    # if Verbosity>0: print '\nPrinting correlationMat (by category)\n', correlationMat_byCat
    # if Verbosity>0: print '\nPrinting covMat\n', covMat


    # ======================================
    # Start decomposition

    decompCovMatObject = ROOT.TDecompSVD( numpy_to_TMatrix(covMat) )
    decompCovMatObject.Decompose()

    Q_TMat      = decompCovMatObject.GetU()
    Q_same_TMat = decompCovMatObject.GetV() # Should be exactly the same as Q
    R_TMat      = decompCovMatObject.GetSig()

    Q      = TMatrix_to_numpy( Q_TMat )
    Q_same = TMatrix_to_numpy( Q_same_TMat )
    R      = TMatrix_to_numpy( R_TMat )

    R = numpy.diag( R.ravel() )


    if Verbosity>0: print '\nPrinting Q\n', Q
    # if Verbosity>0: print '\nPrinting Q_same\n', Q_same
    if Verbosity>0: print '\nPrinting R\n', R

    should_be_zero = (Q.dot(R)).dot( numpy.transpose(Q) ) - covMat
    # if Verbosity>0: print '\nPrinting should_be_zero\n', should_be_zero


    ########################################
    # Computing Atilde and btilde
    ########################################

    chapter( 'Computing Atilde and btilde' )

    Atilde = numpy.zeros((nBins*nCats, nBins))
    btilde = numpy.zeros((nBins*nCats, 1))

    for iCat in xrange(nCats):
        for iBin in xrange(nBins):
            i = iCat*nBins + iBin

            for j in xrange(nBins):

                Atilde[i][j] = (
                    1./sqrt(R[i][i]) *
                    sum([ Q[i][m]*A[m][j] for m in xrange(nBins*nCats) ])
                    )

            btilde[i] = (
                1./sqrt(R[i][i]) *
                sum([ Q[i][m]*b[m] for m in xrange(nBins) ])
                )


    if Verbosity>0: print '\nPrinting Atilde\n', Atilde

    numpy.set_printoptions( precision=5 )
    if Verbosity>0: print '\nPrinting btilde\n', btilde


    # Define C (second derivative matrix)
    C = [ [ 0 for i in xrange(nBins) ] for j in xrange(nBins) ]
    C[0]  = [ -1., 1. ] + [ 0. for i in xrange(nBins-2) ]
    C[-1] = [ 0. for i in xrange(nBins-2) ] + [ 1., -1. ]
    for i in xrange(1,nBins-1):
        C[i] = [ 0. for j in xrange(0,i-1) ] + [ 1., -2., 1. ] + [ 0. for j in xrange(i+2,nBins) ]

    # Add small non-zero component
    for i in xrange(nBins):
        C[i][i] += 0.0001

    C = numpy.array(C)
    if Verbosity>0: print '\nPrinting C\n', C    


    decomposeC = False
    if decomposeC:

        decompCObject = ROOT.TDecompSVD( numpy_to_TMatrix(C) )
        decompCObject.Decompose()

        QC_TMat      = decompCObject.GetU()
        QC_same_TMat = decompCObject.GetV() # Should be exactly the same as Q
        SC_TMat      = decompCObject.GetSig()

        QC      = TMatrix_to_numpy( QC_TMat )
        QC_same = TMatrix_to_numpy( QC_same_TMat )
        SC      = TMatrix_to_numpy( SC_TMat )

        SC = numpy.diag( R.ravel() )


        numpy.set_printoptions( precision=1, linewidth=170 )
        if Verbosity>0: print '\nPrinting QC\n', QC
        if Verbosity>0: print '\nPrinting QC_same\n', QC_same
        if Verbosity>0: print '\nPrinting SC\n', SC
        numpy.set_printoptions( precision=5, linewidth=100 )


    Atilde_times_Cinv = Atilde.dot( numpy.linalg.inv(C) )
    if Verbosity>0: print '\nPrinting Atilde_times_Cinv\n', Atilde_times_Cinv        


    # ======================================
    # Start decomposition

    chapter( 'Decomposing Atilde_times_Cinv' )

    decompTildeObject = ROOT.TDecompSVD( numpy_to_TMatrix(Atilde_times_Cinv) )
    decompTildeObject.Decompose()

    U_TMat = decompTildeObject.GetU()
    V_TMat = decompTildeObject.GetV()
    S_TMat = decompTildeObject.GetSig()

    U                  = TMatrix_to_numpy( U_TMat )
    V                  = TMatrix_to_numpy( V_TMat )
    S_diagonalElements = TMatrix_to_numpy( S_TMat )

    S = numpy.zeros( ( nBins*nCats, nBins ) )
    for iBin in xrange(nBins):
        S[iBin][iBin] = S_diagonalElements[iBin]
    # S = numpy.diag( S.ravel() )


    if Verbosity>0: print '\nPrinting U\n', U
    if Verbosity>0: print '\nPrinting V\n', V
    if Verbosity>0: print '\nPrinting S\n', S


    AtildeCin_SVDtest = ( U.dot(S) ).dot( numpy.transpose(V) ) - Atilde_times_Cinv
    if Verbosity>0: print '\nPrinting AtildeCin_SVDtest\n', AtildeCin_SVDtest


    d = numpy.abs( numpy.transpose(U).dot(btilde) )
    if Verbosity>0: print '\nPrinting d\n', d


    ########################################
    # Try some unfolding
    ########################################

    # PROBLEMATIC!!

    # unfoldTest = RegResult()

    # # Load response matrix into TH2D
    # A_TH2D = ROOT.TH2D(
    #     'A_TH2D', 'A_TH2D',
    #     nCats*nBins, 0, nCats*nBins,
    #     nBins, 0, nBins,
    #     )

    # for iCat in xrange(nCats):
    #     for iBin in xrange(nBins):
    #         i = iCat*nBins + iBin
    #         for j in xrange(nBins):
    #             A_TH2D.SetBinContent( i+1, j+1, A[i][j] )

    # A_TH2D.SetDirectory(0)
    # ROOT.SetOwnership( A_TH2D, False )


    # # Load covariance matrix into TH2D
    # covMat_TH2D = ROOT.TH2D(
    #     'covMat_TH2D', 'covMat_TH2D',
    #     nCats*nBins, 0, nCats*nBins,
    #     nCats*nBins, 0, nCats*nBins,
    #     )

    # for iCat in xrange(nCats):
    #     for iBin in xrange(nBins):
    #         i = iCat*nBins + iBin
    #         for j in xrange(nBins):
    #             covMat_TH2D.SetBinContent( i+1, j+1, covMat[i][j] )

    # covMat_TH2D.SetDirectory(0)
    # ROOT.SetOwnership( covMat_TH2D, False )


    # genMuObserved_TH1D = unfoldTest.OneDimNumpy_to_TH1D( genMuObserved )
    # for iBin in xrange(nBins):
    #     genMuObserved_TH1D.SetBinError( iBin+1, genMuObservedError[iBin] )

    # recoMuObserved_TH1D = unfoldTest.OneDimNumpy_to_TH1D( recoMuObserved )
    # for iBin in xrange(nBins):
    #     recoMuObserved_TH1D.SetBinError( iBin+1, recoMuObservedError[iBin] )

    # genMuExpected_TH1D  = unfoldTest.OneDimNumpy_to_TH1D( numpy.ones(nBins) )
    # recoMuExpected_TH1D = unfoldTest.OneDimNumpy_to_TH1D( numpy.ones(nBins*nCats) )


    # unfoldObject = ROOT.TSVDUnfold(
    #     recoMuObserved_TH1D,
    #     covMat_TH2D,
    #     recoMuExpected_TH1D,
    #     genMuObserved_TH1D,
    #     A_TH2D
    #     )
    # ROOT.SetOwnership( unfoldObject, False )

    # muUnfolded = unfoldObject.Unfold(5)
    # print muUnfolded


    # ########################################
    # # C times mu
    # ########################################

    # C_times_muRecoObserved = C.dot( recoMuObserved )
    # squarednorm_C_times_muRecoObserved = ( numpy.linalg.norm( C_times_muRecoObserved ) )**2

    # if Verbosity>0: print '\nPrinting squarednorm_C_times_muRecoObserved\n', squarednorm_C_times_muRecoObserved


    ########################################
    # Load into simple class (for easy plotting)
    ########################################

    if not plotdir:
        regResult = RegResult()
    else:
        regResult = RegResult( plotdir )

    # regResult.baseDir                     = baseDir
    # regResult.dcDir                       = dcDir
    # regResult.genMuFitDir                 = genMuFitDir
    # regResult.genMuFitFile                = genMuFitFile
    # regResult.recoMuFitDir                = recoMuFitDir
    # regResult.recoMuFitFile               = recoMuFitFile
    regResult.name                        = name
    regResult.nBins                       = nBins
    regResult.nCats                       = nCats
    regResult.Binning                     = Binning
    regResult.genMuObserved               = genMuObserved
    regResult.genMuObservedError          = genMuObservedError
    regResult.recoMuObserved              = recoMuObserved
    regResult.recoMuObservedError         = recoMuObservedError
    regResult.eff_times_acceptance        = eff_times_acceptance
    regResult.summed_eff_times_acceptance = summed_eff_times_acceptance
    regResult.A                           = A
    regResult.x                           = x
    regResult.b                           = b
    regResult.Ax                          = Ax
    regResult.correlationMat              = correlationMat
    regResult.correlationMat_byCat        = correlationMat_byCat
    regResult.covMat                      = covMat
    regResult.genCorrelationMat           = genCorrelationMat
    regResult.genCorrelationMat_byCat     = genCorrelationMat_byCat
    regResult.genCovMat                   = genCovMat
    regResult.Q                           = Q
    regResult.Q_same                      = Q_same
    regResult.R                           = R
    regResult.should_be_zero              = should_be_zero
    regResult.Atilde                      = Atilde
    regResult.btilde                      = btilde
    regResult.C                           = C    
    regResult.Atilde_times_Cinv           = Atilde_times_Cinv        
    regResult.U                           = U
    regResult.V                           = V
    regResult.S                           = S
    regResult.d                           = d


    if name == 'pt':
        regResult.label = 'p_{T} [GeV]'
    else:
        regResult.label = 'nJets'


    if makePlots:

        regResult.plot_d()
        regResult.plot_spectra()
        # regResult.drawCorrMatrix()

        regResult.drawMatrixInRoot( S, 'S' )

        regResult.drawMatrixInRoot( regResult.correlationMat, 'corrMat_RecoMu', overlayCatLines=True, title='{0} correlation matrix (in reco bins)'.format(regResult.label) )
        # regResult.drawMatrixInRoot( regResult.correlationMat_byCat, 'corrMat_RecoMu_ByCat' )
        regResult.drawMatrixInRoot( regResult.genCorrelationMat, 'corrMat_GenMu', title='{0} correlation matrix (in gen bins)'.format(regResult.label) )


    return regResult







def readCorrMatrix( covmatRootFile, nBins, nCats=1, observedErrors=[] ):

    # covmatRootFile = join( recoMuFitDir, 'COVMATISHERE.root' )
    covmatRootFp   = ROOT.TFile.Open( covmatRootFile )
    fit            = covmatRootFp.Get('fit')

    covMat               = numpy.zeros( (nBins*nCats, nBins*nCats) )
    correlationMat       = numpy.zeros( (nBins*nCats, nBins*nCats) )
    correlationMat_byCat = numpy.zeros( (nBins*nCats, nBins*nCats) )

    for iBin1 in xrange(nBins):
        for iCat1 in xrange(nCats):
            for iBin2 in xrange(nBins):
                for iCat2 in xrange(nCats):

                    iBinReco1 = iCat1*nBins+iBin1
                    iBinReco2 = iCat2*nBins+iBin2

                    if nCats == 1:
                        varStr1 = 'r{0}'.format(iBin1)
                        varStr2 = 'r{0}'.format(iBin2)
                    else:
                        varStr1 = 'r{0}_cat{1}'.format(iBin1,iCat1)
                        varStr2 = 'r{0}_cat{1}'.format(iBin2,iCat2)

                    correlationMat[iBinReco1][iBinReco2] = fit.correlation( varStr1, varStr2 )

                    if len(observedErrors) > 0:
                        covMat[iBinReco1][iBinReco2] = (
                            correlationMat[iBinReco1][iBinReco2]
                            * observedErrors[iBinReco1]
                            * observedErrors[iBinReco2]
                            )

                    correlationMat_byCat[iBin1*nCats+iCat1][iBin2*nCats+iCat2] = fit.correlation( varStr1, varStr2 )

    # # Test with some migrations
    # smearing = 0.3
    # correlationMat = [ [ 0 for i in xrange(nBins) ] for j in xrange(nBins) ]
    # correlationMat[0]  = [ 1.-smearing, smearing ] + [ 0. for i in xrange(nBins-2) ]
    # correlationMat[-1] = [ 0. for i in xrange(nBins-2) ] + [ smearing, 1.-smearing ]
    # for i in xrange(1,nBins-1):
    #     correlationMat[i] = [ 0. for j in xrange(0,i-1) ] + [ smearing, 1.-2*smearing, smearing ] + [ 0. for j in xrange(i+2,nBins) ]

    # if Verbosity>0: print '\nPrinting correlationMat\n', correlationMat
    # if Verbosity>0: print '\nPrinting correlationMat (by category)\n', correlationMat_byCat
    
    if len(observedErrors) > 0:
        # if Verbosity>0: print '\nPrinting covMat\n', covMat
        return correlationMat, correlationMat_byCat, covMat
    else:
        return correlationMat, correlationMat_byCat







def chapter( text ):
    if Verbosity>0: print '\n\n=================================================='
    if Verbosity>0: print text
    if Verbosity>0: print



def OneDimNumpy_to_TH1D( self, numpy_hist, Binning=None ):

    dimension = [ i for i in numpy_hist.shape if i!=1 ][0]

    if Binning:
        H = ROOT.TH1D(
            self.rootname(), self.rootname(),
            dimension, array( 'd', Binning )
            )
    else:
        H = ROOT.TH1D(
            'temp', 'temp',
            dimension, 0, dimension,
            )

    for iBin in xrange(dimension):
        H.SetBinContent( iBin+1, numpy_hist[iBin] )

    ROOT.SetOwnership( H, False )

    return H

def numpy_to_TMatrix( M ):
    xdim, ydim = M.shape
    M_TMat = ROOT.TMatrixD( xdim, ydim, array( 'd', M.ravel() ) )
    return M_TMat

def TMatrix_to_numpy( M_TMat ):

    nRows = M_TMat.GetNrows()

    try:
        nCols = M_TMat.GetNcols()
        dim = 2
    except AttributeError:
        nCols = 1
        dim = 1

    M = numpy.zeros((nRows, nCols))

    for iRow in xrange(nRows):
        for iCol in xrange(nCols):
            if dim == 2:
                M[iRow][iCol] = M_TMat[iRow][iCol]
            elif dim == 1:
                M[iRow][iCol] = M_TMat[iRow]

    return M


########################################
# End of Main
########################################
if __name__ == "__main__":
    main()