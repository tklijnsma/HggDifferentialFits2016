#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import argparse
import ROOT

DEBUG = False


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( 'wsFile', type=str, default='default', help='default string' )
    # parser.add_argument( '--boolean', action='store_true', help='boolean')
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )
    args = parser.parse_args()
    wsFile = args.wsFile

    print getCommaSeparatedListOfVariablesToFreeze( wsFile )




def getCommaSeparatedListOfVariablesToFreeze( wsFile ):


    if not os.path.isfile( wsFile ):
        print 'ERROR: File {0} does not exist'.format( wsFile )
        return

    wsFp = ROOT.TFile.Open( wsFile )
    ws   = wsFp.Get('w')

    loadedSnapshot = ws.loadSnapshot('MultiDimFit')

    if not loadedSnapshot:
        print 'ERROR: Could not load MultiDimFit snapshot - are you passing the post fit workspace?'
        return


    # Make a guess for nJets or pT
    if 'njets' in wsFile.lower():
        allCats = allCats_nJets
        allBins = allBins_nJets
        allVars = allVars_nJets
        alwaysFreeze = alwaysFreeze_nJets
    elif 'pt' in wsFile.lower():
        allCats = allCats_pT
        allBins = allBins_pT
        allVars = allVars_pT
        alwaysFreeze = alwaysFreeze_pT


    varsToFloat = []

    for Bin in allBins:
        for Cat in allCats:

            shapeName = 'shapeBkg_bkg_mass_SigmaMpTTag_{0}_{1}_13TeV'.format( Cat, Bin )

            envl = ws.pdf(shapeName)

            dependingVars = ROOT.RooArgList( envl.getPdf( envl.getCurrentIndex() ).getObservables( ws.allVars() ))
            nVars = dependingVars.getSize()

            if DEBUG: print '\n{0} variables:'.format( shapeName )

            for i in xrange(nVars):
                if dependingVars[i].GetName() == 'CMS_hgg_mass': continue
                if DEBUG: print '    ' + dependingVars[i].GetName()
                varsToFloat.append( dependingVars[i].GetName() )


    varsToFreeze = list(set(allVars) - set(varsToFloat))
    varsToFreeze.extend( alwaysFreeze )

    # For some reason, the first variable is never frozen; simply append it again at end of list
    varsToFreeze.append( varsToFreeze[0] )


    if DEBUG:
        print '\n\nFreezing the following variables:'
        for i in varsToFreeze:
            print i

    wsFp.Close()

    return ','.join(varsToFreeze)





allCats_nJets = [
    '0', '1', '2'
    ]

allBins_nJets = [
    'recoNjets2p5_m0p5to0p5',
    'recoNjets2p5_0p5to1p5',
    'recoNjets2p5_1p5to2p5',
    'recoNjets2p5_2p5to3p5',
    'recoNjets2p5_3p5to100p0',
    ]

allCats_pT = [
    '0', '1', '2'
    ]

allBins_pT = [
    'recoPt_0p0to15p0',
    'recoPt_15p0to30p0',
    'recoPt_30p0to45p0',
    'recoPt_45p0to85p0',
    'recoPt_85p0to125p0',
    'recoPt_125p0to200p0',
    'recoPt_200p0to10000p0',
    ]


allVars_nJets = [
    'env_pdf_0_13TeV_bern3_p0',
    'env_pdf_0_13TeV_bern3_p1',
    'env_pdf_0_13TeV_bern3_p2',
    'env_pdf_0_13TeV_bern4_p0',
    'env_pdf_0_13TeV_bern4_p1',
    'env_pdf_0_13TeV_bern4_p2',
    'env_pdf_0_13TeV_bern4_p3',
    'env_pdf_0_13TeV_exp1_p1',
    'env_pdf_0_13TeV_exp3_p1',
    'env_pdf_0_13TeV_exp3_f1',
    'env_pdf_0_13TeV_exp3_p2',
    'env_pdf_0_13TeV_pow1_p1',
    'env_pdf_0_13TeV_lau1_l1',
    'env_pdf_1_13TeV_bern2_p0',
    'env_pdf_1_13TeV_bern2_p1',
    'env_pdf_1_13TeV_bern3_p0',
    'env_pdf_1_13TeV_bern3_p1',
    'env_pdf_1_13TeV_bern3_p2',
    'env_pdf_1_13TeV_exp1_p1',
    'env_pdf_1_13TeV_exp3_p1',
    'env_pdf_1_13TeV_exp3_f1',
    'env_pdf_1_13TeV_exp3_p2',
    'env_pdf_1_13TeV_pow1_p1',
    'env_pdf_1_13TeV_lau1_l1',
    'env_pdf_1_13TeV_lau2_l1',
    'env_pdf_1_13TeV_lau2_h1',
    'env_pdf_2_13TeV_bern1_p0',
    'env_pdf_2_13TeV_bern2_p0',
    'env_pdf_2_13TeV_bern2_p1',
    'env_pdf_2_13TeV_exp1_p1',
    'env_pdf_2_13TeV_exp3_p1',
    'env_pdf_2_13TeV_exp3_f1',
    'env_pdf_2_13TeV_exp3_p2',
    'env_pdf_2_13TeV_pow1_p1',
    'env_pdf_2_13TeV_lau1_l1',
    'env_pdf_2_13TeV_lau2_l1',
    'env_pdf_2_13TeV_lau2_h1',
    'env_pdf_3_13TeV_bern1_p0',
    'env_pdf_3_13TeV_bern2_p0',
    'env_pdf_3_13TeV_bern2_p1',
    'env_pdf_3_13TeV_exp1_p1',
    'env_pdf_3_13TeV_pow1_p1',
    'env_pdf_3_13TeV_lau2_l1',
    'env_pdf_3_13TeV_lau2_h1',
    'env_pdf_4_13TeV_bern1_p0',
    'env_pdf_4_13TeV_bern2_p0',
    'env_pdf_4_13TeV_bern2_p1',
    'env_pdf_4_13TeV_exp1_p1',
    'env_pdf_4_13TeV_pow1_p1',
    'env_pdf_4_13TeV_lau1_l1',
    'env_pdf_4_13TeV_lau2_l1',
    'env_pdf_4_13TeV_lau2_h1',
    'env_pdf_5_13TeV_bern4_p0',
    'env_pdf_5_13TeV_bern4_p1',
    'env_pdf_5_13TeV_bern4_p2',
    'env_pdf_5_13TeV_bern4_p3',
    'env_pdf_5_13TeV_bern5_p0',
    'env_pdf_5_13TeV_bern5_p1',
    'env_pdf_5_13TeV_bern5_p2',
    'env_pdf_5_13TeV_bern5_p3',
    'env_pdf_5_13TeV_bern5_p4',
    'env_pdf_5_13TeV_bern6_p0',
    'env_pdf_5_13TeV_bern6_p1',
    'env_pdf_5_13TeV_bern6_p2',
    'env_pdf_5_13TeV_bern6_p3',
    'env_pdf_5_13TeV_bern6_p4',
    'env_pdf_5_13TeV_bern6_p5',
    'env_pdf_5_13TeV_exp3_p1',
    'env_pdf_5_13TeV_exp3_f1',
    'env_pdf_5_13TeV_exp3_p2',
    'env_pdf_5_13TeV_pow1_p1',
    'env_pdf_5_13TeV_lau1_l1',
    'env_pdf_6_13TeV_bern2_p0',
    'env_pdf_6_13TeV_bern2_p1',
    'env_pdf_6_13TeV_bern3_p0',
    'env_pdf_6_13TeV_bern3_p1',
    'env_pdf_6_13TeV_bern3_p2',
    'env_pdf_6_13TeV_bern4_p0',
    'env_pdf_6_13TeV_bern4_p1',
    'env_pdf_6_13TeV_bern4_p2',
    'env_pdf_6_13TeV_bern4_p3',
    'env_pdf_6_13TeV_exp1_p1',
    'env_pdf_6_13TeV_exp3_p1',
    'env_pdf_6_13TeV_exp3_f1',
    'env_pdf_6_13TeV_exp3_p2',
    'env_pdf_6_13TeV_pow1_p1',
    'env_pdf_6_13TeV_lau1_l1',
    'env_pdf_7_13TeV_bern2_p0',
    'env_pdf_7_13TeV_bern2_p1',
    'env_pdf_7_13TeV_bern3_p0',
    'env_pdf_7_13TeV_bern3_p1',
    'env_pdf_7_13TeV_bern3_p2',
    'env_pdf_7_13TeV_exp1_p1',
    'env_pdf_7_13TeV_exp3_p1',
    'env_pdf_7_13TeV_exp3_f1',
    'env_pdf_7_13TeV_exp3_p2',
    'env_pdf_7_13TeV_pow1_p1',
    'env_pdf_7_13TeV_lau1_l1',
    'env_pdf_7_13TeV_lau2_l1',
    'env_pdf_7_13TeV_lau2_h1',
    'env_pdf_8_13TeV_bern1_p0',
    'env_pdf_8_13TeV_bern2_p0',
    'env_pdf_8_13TeV_bern2_p1',
    'env_pdf_8_13TeV_exp1_p1',
    'env_pdf_8_13TeV_pow1_p1',
    'env_pdf_8_13TeV_lau1_l1',
    'env_pdf_8_13TeV_lau2_l1',
    'env_pdf_8_13TeV_lau2_h1',
    'env_pdf_9_13TeV_bern1_p0',
    'env_pdf_9_13TeV_bern2_p0',
    'env_pdf_9_13TeV_bern2_p1',
    'env_pdf_9_13TeV_exp1_p1',
    'env_pdf_9_13TeV_pow1_p1',
    'env_pdf_9_13TeV_lau2_l1',
    'env_pdf_9_13TeV_lau2_h1',
    'env_pdf_10_13TeV_bern4_p0',
    'env_pdf_10_13TeV_bern4_p1',
    'env_pdf_10_13TeV_bern4_p2',
    'env_pdf_10_13TeV_bern4_p3',
    'env_pdf_10_13TeV_bern5_p0',
    'env_pdf_10_13TeV_bern5_p1',
    'env_pdf_10_13TeV_bern5_p2',
    'env_pdf_10_13TeV_bern5_p3',
    'env_pdf_10_13TeV_bern5_p4',
    'env_pdf_10_13TeV_bern6_p0',
    'env_pdf_10_13TeV_bern6_p1',
    'env_pdf_10_13TeV_bern6_p2',
    'env_pdf_10_13TeV_bern6_p3',
    'env_pdf_10_13TeV_bern6_p4',
    'env_pdf_10_13TeV_bern6_p5',
    'env_pdf_10_13TeV_exp3_p1',
    'env_pdf_10_13TeV_exp3_f1',
    'env_pdf_10_13TeV_exp3_p2',
    'env_pdf_10_13TeV_pow1_p1',
    'env_pdf_10_13TeV_pow3_p1',
    'env_pdf_10_13TeV_pow3_f1',
    'env_pdf_10_13TeV_pow3_p2',
    'env_pdf_10_13TeV_lau1_l1',
    'env_pdf_10_13TeV_lau2_l1',
    'env_pdf_10_13TeV_lau2_h1',
    'env_pdf_10_13TeV_lau3_l1',
    'env_pdf_10_13TeV_lau3_l2',
    'env_pdf_10_13TeV_lau3_h1',
    'env_pdf_11_13TeV_bern3_p0',
    'env_pdf_11_13TeV_bern3_p1',
    'env_pdf_11_13TeV_bern3_p2',
    'env_pdf_11_13TeV_bern4_p0',
    'env_pdf_11_13TeV_bern4_p1',
    'env_pdf_11_13TeV_bern4_p2',
    'env_pdf_11_13TeV_bern4_p3',
    'env_pdf_11_13TeV_bern5_p0',
    'env_pdf_11_13TeV_bern5_p1',
    'env_pdf_11_13TeV_bern5_p2',
    'env_pdf_11_13TeV_bern5_p3',
    'env_pdf_11_13TeV_bern5_p4',
    'env_pdf_11_13TeV_bern6_p0',
    'env_pdf_11_13TeV_bern6_p1',
    'env_pdf_11_13TeV_bern6_p2',
    'env_pdf_11_13TeV_bern6_p3',
    'env_pdf_11_13TeV_bern6_p4',
    'env_pdf_11_13TeV_bern6_p5',
    'env_pdf_11_13TeV_exp3_p1',
    'env_pdf_11_13TeV_exp3_f1',
    'env_pdf_11_13TeV_exp3_p2',
    'env_pdf_11_13TeV_pow1_p1',
    'env_pdf_11_13TeV_lau1_l1',
    'env_pdf_12_13TeV_bern2_p0',
    'env_pdf_12_13TeV_bern2_p1',
    'env_pdf_12_13TeV_bern3_p0',
    'env_pdf_12_13TeV_bern3_p1',
    'env_pdf_12_13TeV_bern3_p2',
    'env_pdf_12_13TeV_exp1_p1',
    'env_pdf_12_13TeV_exp3_p1',
    'env_pdf_12_13TeV_exp3_f1',
    'env_pdf_12_13TeV_exp3_p2',
    'env_pdf_12_13TeV_pow1_p1',
    'env_pdf_12_13TeV_lau1_l1',
    'env_pdf_12_13TeV_lau2_l1',
    'env_pdf_12_13TeV_lau2_h1',
    'env_pdf_13_13TeV_bern2_p0',
    'env_pdf_13_13TeV_bern2_p1',
    'env_pdf_13_13TeV_bern3_p0',
    'env_pdf_13_13TeV_bern3_p1',
    'env_pdf_13_13TeV_bern3_p2',
    'env_pdf_13_13TeV_bern4_p0',
    'env_pdf_13_13TeV_bern4_p1',
    'env_pdf_13_13TeV_bern4_p2',
    'env_pdf_13_13TeV_bern4_p3',
    'env_pdf_13_13TeV_exp1_p1',
    'env_pdf_13_13TeV_exp3_p1',
    'env_pdf_13_13TeV_exp3_f1',
    'env_pdf_13_13TeV_exp3_p2',
    'env_pdf_13_13TeV_pow1_p1',
    'env_pdf_13_13TeV_lau1_l1',
    'env_pdf_13_13TeV_lau2_l1',
    'env_pdf_13_13TeV_lau2_h1',
    'env_pdf_14_13TeV_bern1_p0',
    'env_pdf_14_13TeV_bern2_p0',
    'env_pdf_14_13TeV_bern2_p1',
    'env_pdf_14_13TeV_exp1_p1',
    'env_pdf_14_13TeV_pow1_p1',
    'env_pdf_14_13TeV_lau1_l1',
    'env_pdf_14_13TeV_lau2_l1',
    'env_pdf_14_13TeV_lau2_h1',
    ]

alwaysFreeze_nJets = [
    'pdfindex_SigmaMpTTag_0_recoNjets2p5_m0p5to0p5_13TeV',
    'pdfindex_SigmaMpTTag_0_recoNjets2p5_0p5to1p5_13TeV',
    'pdfindex_SigmaMpTTag_0_recoNjets2p5_1p5to2p5_13TeV',
    'pdfindex_SigmaMpTTag_0_recoNjets2p5_2p5to3p5_13TeV',
    'pdfindex_SigmaMpTTag_0_recoNjets2p5_3p5to100p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoNjets2p5_m0p5to0p5_13TeV',
    'pdfindex_SigmaMpTTag_1_recoNjets2p5_0p5to1p5_13TeV',
    'pdfindex_SigmaMpTTag_1_recoNjets2p5_1p5to2p5_13TeV',
    'pdfindex_SigmaMpTTag_1_recoNjets2p5_2p5to3p5_13TeV',
    'pdfindex_SigmaMpTTag_1_recoNjets2p5_3p5to100p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoNjets2p5_m0p5to0p5_13TeV',
    'pdfindex_SigmaMpTTag_2_recoNjets2p5_0p5to1p5_13TeV',
    'pdfindex_SigmaMpTTag_2_recoNjets2p5_1p5to2p5_13TeV',
    'pdfindex_SigmaMpTTag_2_recoNjets2p5_2p5to3p5_13TeV',
    'pdfindex_SigmaMpTTag_2_recoNjets2p5_3p5to100p0_13TeV',
    ]


allVars_pT = [
    'env_pdf_0_13TeV_bern3_p0',
    'env_pdf_0_13TeV_bern3_p1',
    'env_pdf_0_13TeV_bern3_p2',
    'env_pdf_0_13TeV_bern4_p0',
    'env_pdf_0_13TeV_bern4_p1',
    'env_pdf_0_13TeV_bern4_p2',
    'env_pdf_0_13TeV_bern4_p3',
    'env_pdf_0_13TeV_bern5_p0',
    'env_pdf_0_13TeV_bern5_p1',
    'env_pdf_0_13TeV_bern5_p2',
    'env_pdf_0_13TeV_bern5_p3',
    'env_pdf_0_13TeV_bern5_p4',
    'env_pdf_0_13TeV_bern6_p0',
    'env_pdf_0_13TeV_bern6_p1',
    'env_pdf_0_13TeV_bern6_p2',
    'env_pdf_0_13TeV_bern6_p3',
    'env_pdf_0_13TeV_bern6_p4',
    'env_pdf_0_13TeV_bern6_p5',
    'env_pdf_0_13TeV_exp1_p1',
    'env_pdf_0_13TeV_exp3_p1',
    'env_pdf_0_13TeV_exp3_f1',
    'env_pdf_0_13TeV_exp3_p2',
    'env_pdf_0_13TeV_pow1_p1',
    'env_pdf_0_13TeV_lau1_l1',
    'env_pdf_1_13TeV_bern3_p0',
    'env_pdf_1_13TeV_bern3_p1',
    'env_pdf_1_13TeV_bern3_p2',
    'env_pdf_1_13TeV_bern4_p0',
    'env_pdf_1_13TeV_bern4_p1',
    'env_pdf_1_13TeV_bern4_p2',
    'env_pdf_1_13TeV_bern4_p3',
    'env_pdf_1_13TeV_exp1_p1',
    'env_pdf_1_13TeV_exp3_p1',
    'env_pdf_1_13TeV_exp3_f1',
    'env_pdf_1_13TeV_exp3_p2',
    'env_pdf_1_13TeV_pow1_p1',
    'env_pdf_1_13TeV_lau1_l1',
    'env_pdf_2_13TeV_bern2_p0',
    'env_pdf_2_13TeV_bern2_p1',
    'env_pdf_2_13TeV_bern3_p0',
    'env_pdf_2_13TeV_bern3_p1',
    'env_pdf_2_13TeV_bern3_p2',
    'env_pdf_2_13TeV_bern4_p0',
    'env_pdf_2_13TeV_bern4_p1',
    'env_pdf_2_13TeV_bern4_p2',
    'env_pdf_2_13TeV_bern4_p3',
    'env_pdf_2_13TeV_bern5_p0',
    'env_pdf_2_13TeV_bern5_p1',
    'env_pdf_2_13TeV_bern5_p2',
    'env_pdf_2_13TeV_bern5_p3',
    'env_pdf_2_13TeV_bern5_p4',
    'env_pdf_2_13TeV_bern6_p0',
    'env_pdf_2_13TeV_bern6_p1',
    'env_pdf_2_13TeV_bern6_p2',
    'env_pdf_2_13TeV_bern6_p3',
    'env_pdf_2_13TeV_bern6_p4',
    'env_pdf_2_13TeV_bern6_p5',
    'env_pdf_2_13TeV_exp1_p1',
    'env_pdf_2_13TeV_exp3_p1',
    'env_pdf_2_13TeV_exp3_f1',
    'env_pdf_2_13TeV_exp3_p2',
    'env_pdf_2_13TeV_pow1_p1',
    'env_pdf_2_13TeV_lau1_l1',
    'env_pdf_3_13TeV_bern2_p0',
    'env_pdf_3_13TeV_bern2_p1',
    'env_pdf_3_13TeV_bern3_p0',
    'env_pdf_3_13TeV_bern3_p1',
    'env_pdf_3_13TeV_bern3_p2',
    'env_pdf_3_13TeV_exp1_p1',
    'env_pdf_3_13TeV_pow1_p1',
    'env_pdf_3_13TeV_lau1_l1',
    'env_pdf_4_13TeV_bern1_p0',
    'env_pdf_4_13TeV_bern2_p0',
    'env_pdf_4_13TeV_bern2_p1',
    'env_pdf_4_13TeV_exp1_p1',
    'env_pdf_4_13TeV_pow1_p1',
    'env_pdf_4_13TeV_lau1_l1',
    'env_pdf_5_13TeV_bern1_p0',
    'env_pdf_5_13TeV_bern2_p0',
    'env_pdf_5_13TeV_bern2_p1',
    'env_pdf_5_13TeV_exp1_p1',
    'env_pdf_5_13TeV_pow1_p1',
    'env_pdf_5_13TeV_lau1_l1',
    'env_pdf_6_13TeV_bern1_p0',
    'env_pdf_6_13TeV_exp1_p1',
    'env_pdf_6_13TeV_pow1_p1',
    'env_pdf_6_13TeV_lau1_l1',
    'env_pdf_6_13TeV_lau2_l1',
    'env_pdf_6_13TeV_lau2_h1',
    'env_pdf_7_13TeV_bern3_p0',
    'env_pdf_7_13TeV_bern3_p1',
    'env_pdf_7_13TeV_bern3_p2',
    'env_pdf_7_13TeV_bern4_p0',
    'env_pdf_7_13TeV_bern4_p1',
    'env_pdf_7_13TeV_bern4_p2',
    'env_pdf_7_13TeV_bern4_p3',
    'env_pdf_7_13TeV_exp3_p1',
    'env_pdf_7_13TeV_exp3_f1',
    'env_pdf_7_13TeV_exp3_p2',
    'env_pdf_7_13TeV_pow1_p1',
    'env_pdf_7_13TeV_lau1_l1',
    'env_pdf_8_13TeV_bern3_p0',
    'env_pdf_8_13TeV_bern3_p1',
    'env_pdf_8_13TeV_bern3_p2',
    'env_pdf_8_13TeV_bern4_p0',
    'env_pdf_8_13TeV_bern4_p1',
    'env_pdf_8_13TeV_bern4_p2',
    'env_pdf_8_13TeV_bern4_p3',
    'env_pdf_8_13TeV_exp3_p1',
    'env_pdf_8_13TeV_exp3_f1',
    'env_pdf_8_13TeV_exp3_p2',
    'env_pdf_8_13TeV_pow1_p1',
    'env_pdf_8_13TeV_lau1_l1',
    'env_pdf_9_13TeV_bern3_p0',
    'env_pdf_9_13TeV_bern3_p1',
    'env_pdf_9_13TeV_bern3_p2',
    'env_pdf_9_13TeV_bern4_p0',
    'env_pdf_9_13TeV_bern4_p1',
    'env_pdf_9_13TeV_bern4_p2',
    'env_pdf_9_13TeV_bern4_p3',
    'env_pdf_9_13TeV_exp1_p1',
    'env_pdf_9_13TeV_exp3_p1',
    'env_pdf_9_13TeV_exp3_f1',
    'env_pdf_9_13TeV_exp3_p2',
    'env_pdf_9_13TeV_pow1_p1',
    'env_pdf_9_13TeV_lau1_l1',
    'env_pdf_10_13TeV_bern2_p0',
    'env_pdf_10_13TeV_bern2_p1',
    'env_pdf_10_13TeV_bern3_p0',
    'env_pdf_10_13TeV_bern3_p1',
    'env_pdf_10_13TeV_bern3_p2',
    'env_pdf_10_13TeV_exp1_p1',
    'env_pdf_10_13TeV_exp3_p1',
    'env_pdf_10_13TeV_exp3_f1',
    'env_pdf_10_13TeV_exp3_p2',
    'env_pdf_10_13TeV_pow1_p1',
    'env_pdf_10_13TeV_lau1_l1',
    'env_pdf_10_13TeV_lau2_l1',
    'env_pdf_10_13TeV_lau2_h1',
    'env_pdf_11_13TeV_bern1_p0',
    'env_pdf_11_13TeV_bern2_p0',
    'env_pdf_11_13TeV_bern2_p1',
    'env_pdf_11_13TeV_exp1_p1',
    'env_pdf_11_13TeV_pow1_p1',
    'env_pdf_11_13TeV_lau2_l1',
    'env_pdf_11_13TeV_lau2_h1',
    'env_pdf_12_13TeV_bern1_p0',
    'env_pdf_12_13TeV_exp1_p1',
    'env_pdf_12_13TeV_pow1_p1',
    'env_pdf_12_13TeV_lau2_l1',
    'env_pdf_12_13TeV_lau2_h1',
    'env_pdf_13_13TeV_bern1_p0',
    'env_pdf_13_13TeV_exp1_p1',
    'env_pdf_13_13TeV_pow1_p1',
    'env_pdf_13_13TeV_lau2_l1',
    'env_pdf_13_13TeV_lau2_h1',
    'env_pdf_14_13TeV_bern4_p0',
    'env_pdf_14_13TeV_bern4_p1',
    'env_pdf_14_13TeV_bern4_p2',
    'env_pdf_14_13TeV_bern4_p3',
    'env_pdf_14_13TeV_bern5_p0',
    'env_pdf_14_13TeV_bern5_p1',
    'env_pdf_14_13TeV_bern5_p2',
    'env_pdf_14_13TeV_bern5_p3',
    'env_pdf_14_13TeV_bern5_p4',
    'env_pdf_14_13TeV_bern6_p0',
    'env_pdf_14_13TeV_bern6_p1',
    'env_pdf_14_13TeV_bern6_p2',
    'env_pdf_14_13TeV_bern6_p3',
    'env_pdf_14_13TeV_bern6_p4',
    'env_pdf_14_13TeV_bern6_p5',
    'env_pdf_14_13TeV_exp3_p1',
    'env_pdf_14_13TeV_exp3_f1',
    'env_pdf_14_13TeV_exp3_p2',
    'env_pdf_14_13TeV_pow1_p1',
    'env_pdf_14_13TeV_pow3_p1',
    'env_pdf_14_13TeV_pow3_f1',
    'env_pdf_14_13TeV_pow3_p2',
    'env_pdf_14_13TeV_lau1_l1',
    'env_pdf_15_13TeV_bern3_p0',
    'env_pdf_15_13TeV_bern3_p1',
    'env_pdf_15_13TeV_bern3_p2',
    'env_pdf_15_13TeV_bern4_p0',
    'env_pdf_15_13TeV_bern4_p1',
    'env_pdf_15_13TeV_bern4_p2',
    'env_pdf_15_13TeV_bern4_p3',
    'env_pdf_15_13TeV_bern5_p0',
    'env_pdf_15_13TeV_bern5_p1',
    'env_pdf_15_13TeV_bern5_p2',
    'env_pdf_15_13TeV_bern5_p3',
    'env_pdf_15_13TeV_bern5_p4',
    'env_pdf_15_13TeV_exp3_p1',
    'env_pdf_15_13TeV_exp3_f1',
    'env_pdf_15_13TeV_exp3_p2',
    'env_pdf_15_13TeV_pow1_p1',
    'env_pdf_15_13TeV_lau1_l1',
    'env_pdf_16_13TeV_bern3_p0',
    'env_pdf_16_13TeV_bern3_p1',
    'env_pdf_16_13TeV_bern3_p2',
    'env_pdf_16_13TeV_bern4_p0',
    'env_pdf_16_13TeV_bern4_p1',
    'env_pdf_16_13TeV_bern4_p2',
    'env_pdf_16_13TeV_bern4_p3',
    'env_pdf_16_13TeV_bern5_p0',
    'env_pdf_16_13TeV_bern5_p1',
    'env_pdf_16_13TeV_bern5_p2',
    'env_pdf_16_13TeV_bern5_p3',
    'env_pdf_16_13TeV_bern5_p4',
    'env_pdf_16_13TeV_exp3_p1',
    'env_pdf_16_13TeV_exp3_f1',
    'env_pdf_16_13TeV_exp3_p2',
    'env_pdf_16_13TeV_pow1_p1',
    'env_pdf_16_13TeV_lau1_l1',
    'env_pdf_17_13TeV_bern3_p0',
    'env_pdf_17_13TeV_bern3_p1',
    'env_pdf_17_13TeV_bern3_p2',
    'env_pdf_17_13TeV_bern4_p0',
    'env_pdf_17_13TeV_bern4_p1',
    'env_pdf_17_13TeV_bern4_p2',
    'env_pdf_17_13TeV_bern4_p3',
    'env_pdf_17_13TeV_bern5_p0',
    'env_pdf_17_13TeV_bern5_p1',
    'env_pdf_17_13TeV_bern5_p2',
    'env_pdf_17_13TeV_bern5_p3',
    'env_pdf_17_13TeV_bern5_p4',
    'env_pdf_17_13TeV_bern6_p0',
    'env_pdf_17_13TeV_bern6_p1',
    'env_pdf_17_13TeV_bern6_p2',
    'env_pdf_17_13TeV_bern6_p3',
    'env_pdf_17_13TeV_bern6_p4',
    'env_pdf_17_13TeV_bern6_p5',
    'env_pdf_17_13TeV_exp3_p1',
    'env_pdf_17_13TeV_exp3_f1',
    'env_pdf_17_13TeV_exp3_p2',
    'env_pdf_17_13TeV_pow1_p1',
    'env_pdf_17_13TeV_lau1_l1',
    'env_pdf_18_13TeV_bern1_p0',
    'env_pdf_18_13TeV_bern2_p0',
    'env_pdf_18_13TeV_bern2_p1',
    'env_pdf_18_13TeV_exp1_p1',
    'env_pdf_18_13TeV_pow1_p1',
    'env_pdf_18_13TeV_lau2_l1',
    'env_pdf_18_13TeV_lau2_h1',
    'env_pdf_19_13TeV_bern1_p0',
    'env_pdf_19_13TeV_exp1_p1',
    'env_pdf_19_13TeV_pow1_p1',
    'env_pdf_19_13TeV_lau2_l1',
    'env_pdf_19_13TeV_lau2_h1',
    'env_pdf_20_13TeV_bern1_p0',
    'env_pdf_20_13TeV_exp1_p1',
    'env_pdf_20_13TeV_pow1_p1',
    'env_pdf_20_13TeV_lau2_l1',
    'env_pdf_20_13TeV_lau2_h1',
    ]


alwaysFreeze_pT = [
    'pdfindex_SigmaMpTTag_0_recoPt_0p0to15p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_15p0to30p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_30p0to45p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_45p0to85p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_85p0to125p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_125p0to200p0_13TeV',
    'pdfindex_SigmaMpTTag_0_recoPt_200p0to10000p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_0p0to15p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_15p0to30p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_30p0to45p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_45p0to85p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_85p0to125p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_125p0to200p0_13TeV',
    'pdfindex_SigmaMpTTag_1_recoPt_200p0to10000p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_0p0to15p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_15p0to30p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_30p0to45p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_45p0to85p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_85p0to125p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_125p0to200p0_13TeV',
    'pdfindex_SigmaMpTTag_2_recoPt_200p0to10000p0_13TeV',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_0p0to15p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_15p0to30p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_30p0to45p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_45p0to85p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_85p0to125p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_125p0to200p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_0_recoPt_200p0to10000p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_0p0to15p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_15p0to30p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_30p0to45p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_45p0to85p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_85p0to125p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_125p0to200p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_1_recoPt_200p0to10000p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_0p0to15p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_15p0to30p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_30p0to45p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_45p0to85p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_85p0to125p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_125p0to200p0_13TeV__norm',
    # 'shapeBkg_bkg_mass_SigmaMpTTag_2_recoPt_200p0to10000p0_13TeV__norm',
    # 'n_exp_binSigmaMpTTag_0_recoPt_0p0to15p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_15p0to30p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_30p0to45p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_45p0to85p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_85p0to125p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_125p0to200p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_0_recoPt_200p0to10000p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_0p0to15p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_15p0to30p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_30p0to45p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_45p0to85p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_85p0to125p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_125p0to200p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_1_recoPt_200p0to10000p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_0p0to15p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_15p0to30p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_30p0to45p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_45p0to85p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_85p0to125p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_125p0to200p0_13TeV_proc_bkg_mass',
    # 'n_exp_binSigmaMpTTag_2_recoPt_200p0to10000p0_13TeV_proc_bkg_mass',
    ]



########################################
# End of Main
########################################
if __name__ == "__main__":
    main()