#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os
import numpy
from array import array
from math import log, exp, sqrt, pi

from time import strftime
datestr = strftime( '%b%d' )

import matplotlib.pyplot as plt

import ROOT
ROOT.gROOT.SetBatch(True)
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
# ROOT.gSystem.Load('/mnt/t3nfs01/data01/shome/tklijnsm/Combination/RooUnfold-1.1.1/libRooUnfold')


########################################
# Main
########################################

def main():

    regTest = RegResult( plotdir = 'plots_{0}_test'.format(datestr) )
    regTest.drawMatrixInRoot()



# Container class for all the results
class RegResult():

    c = ROOT.TCanvas( 'c', 'c', 1000, 800 )
    ROOT.SetOwnership( c, False )
    ROOTCOUNTER = 1000


    def __init__( self, plotdir = 'plots_{0}'.format(datestr) ):
        self.plotdir = plotdir
        self.SetColorPalette()



    def linSpace( self, left, right, nPoints ):
        dx = (right-left)/float(nPoints-1)
        return [ left + i*dx for i in xrange(nPoints) ]

    def SetColorPalette( self ):
        # n_stops = 20
        # stops  = [ 0.0 ] + self.linSpace( 0.3, 1.0, n_stops-1 )
        # greens = [ 1.0 ] + self.linSpace( 0.5, 0.3, n_stops-1 )
        # blues  = [ 1.0 ] + self.linSpace( 0.5, 0.3, n_stops-1 )
        # reds   = [ 1.0 ] + self.linSpace( 1.0, 0.7, n_stops-1 )

        # n_stops = 50
        # stops  = self.linSpace( 0.0, 0.5, n_stops ) + self.linSpace( 0.5, 1.0, n_stops )[1:]
        # greens = self.linSpace( 0.0, 1.0, n_stops ) + self.linSpace( 1.0, 0.0, n_stops )[1:]
        # blues  = self.linSpace( 1.0, 1.0, n_stops ) + self.linSpace( 0.0, 0.0, n_stops )[1:]
        # reds   = self.linSpace( 0.0, 0.0, n_stops ) + self.linSpace( 1.0, 1.0, n_stops )[1:]

        # l = lambda x, y, z: self.linSpace( x, y, z )

        # n_stops = 20
        # stops   = l( 0.0, 1.0, 20 )
        # greens  = l( 0.0, 0.0, 5 ) + l( 0.0, 1.0, 3 )  + l( 1.0, 1.0, 4 ) + l( 1.0, 0.0, 3 ) + l( 0.0, 0.0, 5 )
        # reds    = l( 0.0, 0.0, 5 )  + l( 0.5, 1.0, 5 )  + l( 1.0, 1.0, 10 )
        # blues   = l( 1.0, 1.0, 10 ) + l( 1.0, 0.5, 10 ) + l( 0.0, 0.0, 5 )

        # ROOT.TColor.CreateGradientColorTable(
        #     n_stops,
        #     array('d', stops ),
        #     array('d', reds ),
        #     array('d', greens ),
        #     array('d', blues ),
        #     109 )

        # ROOT.gStyle.SetPalette( ROOT.TColor.kLightTemperature )
        ROOT.gStyle.SetPalette( 1 )



    def get_k( self ):

        d_TH1D = self.OneDimNumpy_to_TH1D( self.d )

        fitFunction = ROOT.TF1(
            'dfit',
            '[0]*exp(-[1]*x) + [2]',
            0., float( self.nCats * self.nBins )
            )
        ROOT.SetOwnership( fitFunction, False )

        d_TH1D.Fit( fitFunction, 'q' )

        p0 = fitFunction.GetParameter(0)
        p1 = fitFunction.GetParameter(1)
        p2 = fitFunction.GetParameter(2)

        threshold = 2. * p2
        xThreshold = -log(( threshold - p2 )/p0)/p1

        return xThreshold


    def plot_d( self ):

        self.c.Clear()
        self.c.SetLogy()

        d_TH1D = self.OneDimNumpy_to_TH1D( self.d )
        d_TH1D.SetTitle( 'log|d|' )
        d_TH1D.Draw()

        self.c.Update()
        
        
        print '\n\nDoing fit with exponential on d plot'

        fitFunction = ROOT.TF1(
            'dfit',
            '[0]*exp(-[1]*x) + [2]',
            0., float( self.nCats * self.nBins )
            )
        ROOT.SetOwnership( fitFunction, False )

        d_TH1D.Fit( fitFunction )


        fitFunction.SetLineColor(9)
        fitFunction.SetLineWidth(3)
        fitFunction.Draw('SAMEL')




        p0 = fitFunction.GetParameter(0)
        p1 = fitFunction.GetParameter(1)
        p2 = fitFunction.GetParameter(2)


        # expFunction = ROOT.TF1(
        #     'expfn',
        #     '[0]*exp(-[1]*x)',
        #     0., float( self.nCats * self.nBins )
        #     )
        # ROOT.SetOwnership( expFunction, False )
        # expFunction.SetParameter( 0, p0 )
        # expFunction.SetParameter( 1, p1 )
        # expFunction.Draw('LSAME')


        l = ROOT.TLine( 0., p2, self.nCats * self.nBins, p2 )
        l.SetLineWidth(2)
        l.SetLineColor(2)
        l.Draw()


        # Checking at which bin we cross certain threshold

        self.c.Update()
        self.c.cd()
        yMax = 10**( ROOT.gPad.GetUymax() )
        yMin = 10**( ROOT.gPad.GetUymin() )

        # yMin = 0.
        # # yMax = 100000.
        # yMax = d_TH1D.GetYaxis().GetXmax()

        # print 'Drawing line from {0} to {1}'.format( yMin, yMax )

        threshold = 2. * p2
        xThreshold = -log(( threshold - p2 )/p0)/p1

        lThreshold = ROOT.TLine( xThreshold, yMin,xThreshold, yMax )
        ROOT.SetOwnership( lThreshold, False )
        lThreshold.SetLineStyle(2)
        lThreshold.Draw()


        self.Save( 'd' )
        self.c.SetLogy(False)
        




    def plot_spectra( self ):

        self.c.Clear()

        base = self.getBase()

        xMin = self.Binning[0]
        xMax = self.Binning[-1]
        yMin = 0.
        yMax = 1.5

        base.SetXTitle(self.label)
        base.GetXaxis().SetTitleSize(0.04)
        base.GetXaxis().SetTitleOffset(1.07)
        base.SetYTitle('#mu')
        base.GetYaxis().SetTitleSize(0.04)
        base.GetYaxis().SetTitleOffset(1.25)
        base.SetMinimum(yMin)
        base.SetMaximum(yMax)
        base.GetXaxis().SetLimits( xMin, xMax )

        base.Draw('P')


        for bound in self.Binning[1:-1]:
            l = ROOT.TLine( bound, yMin, bound, yMax )
            ROOT.SetOwnership( l, False )
            l.SetLineColor(17)
            l.Draw()


        H_truth = self.OneDimNumpy_to_TH1D( numpy.ones(self.nBins), Binning=self.Binning )
        H_truth.Draw('SAME')
        H_truth.SetLineStyle(1)
        H_truth.SetLineWidth(3)
        H_truth.SetLineColor(1)


        # H_obs = self.OneDimNumpy_to_TH1D( self.genMuObserved, Binning=self.Binning )
        # for iBin in xrange(self.nBins):
        #     H_obs.SetBinError( iBin+1, self.genMuObservedError[iBin] )
        # H_obs.Draw('SAME')
        # H_obs.SetLineStyle(2)
        # H_obs.SetLineWidth(2)
        # H_obs.SetLineColor(2)


        recoMus = numpy.split( self.recoMuObserved, self.nCats )
        recoMuErrors = numpy.split( self.recoMuObservedError, self.nCats )
        # print 'recoMus       = ', recoMus
        # print 'recoMuErrors  = ', recoMuErrors        

        leg = ROOT.TLegend( 0.6, 0.1, 0.9, 0.4 )
        leg.AddEntry( H_truth.GetName(), 'Truth' , 'l' )
        # leg.AddEntry( H_obs.GetName(),   'Data'  , 'l' )
        leg.Draw()




        ROOT.gStyle.SetEndErrorSize(5)


        colors = [ 9, 8, 40, 42 ]
        for iCat in xrange(self.nCats):
            recoMu = recoMus[iCat]
            recoMuError = recoMuErrors[iCat]

            H = self.OneDimNumpy_to_TH1D( recoMu, Binning=self.Binning )

            for iBin in xrange(self.nBins):
                H.SetBinError( iBin+1, recoMuError[iBin] )

            H.SetLineWidth(2)
            H.SetLineColor( colors.pop(0) )
            H.SetMarkerSize(0.6)

            H.Draw('SAMEE1')
            leg.AddEntry( H.GetName(), 'Cat {0}'.format(iCat), 'l' )


        self.Save( 'spectra' )

        ROOT.gStyle.SetEndErrorSize()







    def compareRegularizedSpectrum( self, regularizedGenMuObserved, regularizedGenMuObservedError ):



        self.c.Clear()

        base = self.getBase()

        xMin = self.Binning[0]
        xMax = self.Binning[-1]
        yMin = 0.
        yMax = 1.5

        base.SetXTitle(self.label)
        base.GetXaxis().SetTitleSize(0.04)
        base.GetXaxis().SetTitleOffset(1.07)
        base.SetYTitle('#mu')
        base.GetYaxis().SetTitleSize(0.04)
        base.GetYaxis().SetTitleOffset(1.25)
        base.SetMinimum(yMin)
        base.SetMaximum(yMax)
        base.GetXaxis().SetLimits( xMin, xMax )

        base.Draw('P')


        for bound in self.Binning[1:-1]:
            l = ROOT.TLine( bound, yMin, bound, yMax )
            ROOT.SetOwnership( l, False )
            l.SetLineColor(17)
            l.Draw()


        H_truth = self.OneDimNumpy_to_TH1D( numpy.ones(self.nBins), Binning=self.Binning )
        H_truth.Draw('SAME')
        H_truth.SetLineStyle(1)
        H_truth.SetLineWidth(3)
        H_truth.SetLineColor(1)


        H_obs = self.OneDimNumpy_to_TH1D( self.genMuObserved, Binning=self.Binning )
        for iBin in xrange(self.nBins):
            H_obs.SetBinError( iBin+1, self.genMuObservedError[iBin] )
        H_obs.Draw('SAMEE')
        H_obs.SetLineStyle(1)
        H_obs.SetLineWidth(2)
        H_obs.SetLineColor(2)


        H_reg = self.OneDimNumpy_to_TH1D( regularizedGenMuObserved, Binning=self.Binning )
        for iBin in xrange(self.nBins):
            H_reg.SetBinError( iBin+1, regularizedGenMuObservedError[iBin] )
        H_reg.Draw('SAMEE2')
        H_reg.SetLineStyle(1)
        H_reg.SetLineWidth(2)
        H_reg.SetLineColor(9)
        H_reg.SetFillStyle(3004)
        H_reg.SetFillColor(9)


        leg = ROOT.TLegend( 0.6, 0.1, 0.9, 0.4 )
        leg.AddEntry( H_truth.GetName(), 'Truth' , 'l' )
        leg.AddEntry( H_obs.GetName(),   'Data'  , 'l' )
        leg.AddEntry( H_reg.GetName(),   'Data (regul.)'  , 'f' )
        leg.Draw()


        self.Save( 'spectraCompared' )


            



    def drawMatrix( self, M, name='corrMat', majorTickDistance=None ):

        if not majorTickDistance:
            majorTickDistance = self.nBins

        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1)                     

        plt.rcParams['image.cmap'] = 'bwr' 
        plt.matshow( M, vmin=-1, vmax=1,fignum=0)
        plt.colorbar()

        major_ticks = numpy.arange(0, self.nBins*self.nCats, majorTickDistance )
        minor_ticks = numpy.arange(0, self.nBins*self.nCats, 1 )

        # plt.xticks( np.arange(len(nBins*nCats)), nBins*nCats, rotation='vertical' )
        # plt.yticks( np.arange(len(nBins*nCats)), nBins*nCats );

        ax.set_xticks( major_ticks )                                                       
        ax.set_xticks( minor_ticks, minor=True)                                           
        ax.set_yticks( major_ticks )                                                       
        ax.set_yticks( minor_ticks, minor=True) 

        if not os.path.isdir(self.plotdir): os.makedirs(self.plotdir)
        plt.savefig( os.path.join( self.plotdir, name + '.png' ) )
        plt.clf()


    def drawCorrMatrix( self ):
        # self.drawMatrix( self.correlationMat )
        # self.drawMatrix( self.correlationMat_byCat, 'corrMatt_byCat', self.nCats )

        self.drawMatrixInRoot( self.correlationMat, 'corrMat_root', overlayCatLines=True )
        self.drawMatrixInRoot( self.correlationMat_byCat, 'corrMatByCat_root' )



    def OneDimNumpy_to_TH1D( self, numpy_hist, Binning=None ):

        dimension = [ i for i in numpy_hist.shape if i!=1 ][0]

        if Binning:
            H = ROOT.TH1D(
                self.rootname(), self.rootname(),
                dimension, array( 'd', Binning )
                )
        else:
            H = ROOT.TH1D(
                self.rootname(), self.rootname(),
                dimension, 0, dimension,
                )

        for iBin in xrange(dimension):
            H.SetBinContent( iBin+1, numpy_hist[iBin] )

        ROOT.SetOwnership( H, False )

        return H








    def drawMatrixInRoot( self, M, outname, overlayCatLines=False, title='' ):

        # M = numpy.arange(9).reshape(3,3)

        nRows, nCols = M.shape


        self.c.Clear()

        # p = self.getPad(
        #     User_x_min = 0.0,
        #     User_x_max = nCols,
        #     User_y_min = 0.0,
        #     User_y_max = nRows,
        #     )

        T = ROOT.TH2D(
            self.rootname(), title,
            nCols, 0., nCols,
            nRows, 0., nRows
            )
        ROOT.SetOwnership( T, False )
        T.SetContour(100)

        for iRow in xrange(nRows):
            for iCol in xrange(nCols):
                T.SetBinContent( iCol+1, iRow+1, M[iRow][iCol] )


        if 'cor' in outname.lower():
            T.GetZaxis().SetRangeUser(-1.0,1.0)


        # TH2D *h_z = new TH2D(*h); // create a "copy"
        # h_z->GetZaxis()->SetRangeUser(0.89, 1.11); // ... set the range ...
        # h_z->Draw("z"); // draw "axes", "color palette", "statistics box"

        # zScale = ROOT.TH2D( T )
        # # zScale.GetZaxis()
        # zScale.Draw('Z')

        ROOT.gStyle.SetPaintTextFormat('4.2g')

        T.Draw('COLZ TEXT')

        # Draw grid
        # ls = []
        # for iRow in xrange(nRows+1):
        #     l = ROOT.TLine( 0., float(iRow)/nRows, 1., float(iRow)/nRows )
        #     ROOT.SetOwnership( l, False )
        #     l.Draw()
        #     ls.append( l )
        # for iCol in xrange(nCols+1):
        #     l = ROOT.TLine( float(iCol)/nCols, 0., float(iCol)/nCols, 1. )
        #     ROOT.SetOwnership( l, False )
        #     l.Draw()
        #     ls.append( l )


        if overlayCatLines:

            # Draw grid
            ls = []
            for iRow in xrange(nRows+1):
                if not iRow % self.nBins :
                    l = ROOT.TLine( 0., float(iRow), nRows, float(iRow) )
                    l.SetLineWidth(2)
                    ROOT.SetOwnership( l, False )
                    l.Draw()
                    ls.append( l )
            for iCol in xrange(nCols+1):
                if not iCol % self.nBins :
                    l = ROOT.TLine( float(iCol), 0., float(iCol), nCols )
                    l.SetLineWidth(2)
                    ROOT.SetOwnership( l, False )
                    l.Draw()
                    ls.append( l )



        self.c.cd()
        self.c.Update()

        self.Save( outname )


    def getBase( self ):
        base = ROOT.TH1F()
        base.SetName( self.rootname() )
        base.SetTitle('')
        base.SetMarkerSize(0)
        base.SetMarkerColor(0)
        ROOT.SetOwnership( base, False )
        return base

    def getPad( self,
        pad_name = None,

        NDC_x_min = 0.1,
        NDC_y_min = 0.1,
        NDC_x_max = 0.9,
        NDC_y_max = 0.9,

        User_x_min = 0.0,
        User_x_max = 1.0,
        User_y_min = 0.0,
        User_y_max = 1.0,
        ):

        if not pad_name: pad_name = self.rootname()

        # print '    Drawing {0} from ( {1}, {2} ) to ( {3}, {4} )'.format(
        #     pad_name,
        #     NDC_x_min,
        #     NDC_y_min,
        #     NDC_x_max,
        #     NDC_y_max,
        #     )

        # Define the actual pad
        TPad = ROOT.TPad(
            pad_name, pad_name,
            NDC_x_min,
            NDC_y_min,
            NDC_x_max,
            NDC_y_max,
            )

        # print '    Setting ranges {0} from ( {1}, {2} ) to ( {3}, {4} )'.format(
        #     pad_name,
        #     User_x_min,
        #     User_y_min,
        #     User_x_max,
        #     User_y_max,
        #     )

        TPad.Range(
            User_x_min,
            User_y_min,
            User_x_max,
            User_y_max,
            )

        TPad.Draw()
        ROOT.SetOwnership( TPad, False )
        self.c.Update()
        TPad.cd()

        # Set horizontal margins to zero
        TPad.SetLeftMargin(0.0)
        TPad.SetRightMargin(0.0)
        TPad.SetTopMargin(0.0)
        TPad.SetBottomMargin(0.0)
        TPad.SetFillStyle(4000)

        return TPad


    def rootname( self ):
        newname = 'robj{0}'.format(self.ROOTCOUNTER)
        self.ROOTCOUNTER += 1
        return newname


    def Save( self, outname ):

        outname = os.path.basename(outname).replace('.pdf','').replace('.png','')

        if not os.path.isdir(self.plotdir): os.makedirs(self.plotdir)

        outname = os.path.join( self.plotdir, outname )

        self.c.SaveAs( outname + '.pdf' )



########################################
# End of Main
########################################
if __name__ == "__main__":
    main()