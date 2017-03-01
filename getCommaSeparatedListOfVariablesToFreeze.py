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

    varsToFreeze = []
    
    # Find all pdf indexes 
    allCats = ws.allCats()
    catitr = allCats.createIterator()
    cat = catitr.Next()
    while cat:
        if cat.GetName().startswith("pdfindex"):
            varsToFreeze.append( cat.GetName() )
        cat = catitr.Next()
        
    # Find all background pdfs
    allPdfs = ws.allPdfs()
    pdfitr = allPdfs.createIterator()
    pdf = pdfitr.Next()
    while pdf:
        if pdf.GetName().startswith("shapeBkg_bkg"):
            # Loop over all shapes in the envelope
            for ishape in xrange(pdf.getNumPdfs()):
                # Freeze all pdf parameters except those from the best fit function
                if ishape == pdf.getCurrentIndex(): continue
                shape = pdf.getPdf( ishape )
                observables = ROOT.RooArgList(shape.getObservables( ws.allVars() ))
                observables = filter(lambda x: not x.startswith('CMS_hgg_mass'), map(lambda x: observables[x].GetName(), xrange(observables.getSize()) ) )
                varsToFreeze.extend( observables )
        pdf = pdfitr.Next()

    # For some reason, the first variable is never frozen; simply append it again at end of list
    varsToFreeze.append( varsToFreeze[0] )


    if DEBUG:
        print '\n\nFreezing the following variables:'
        for i in varsToFreeze:
            print i

    wsFp.Close()

    return ','.join(varsToFreeze)



########################################
# End of Main
########################################
if __name__ == "__main__":
    main()
