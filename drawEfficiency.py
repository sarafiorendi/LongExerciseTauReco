from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i"  , "--input"     , dest = "input"     ,  help = "input file"       , default = ''                        )
parser.add_argument("-c"  , "--compare"   , dest = "compfile"  ,  help = "file to compare"  , default = ''                        )
parser.add_argument("-d"  , "--diff"      , dest = "diff"      ,  help = "plot differences" , default = False, action='store_true')
parser.add_argument("-m"  , "--mc"        , dest = "mc"        ,  help = "comparison is mc" , default = False, action='store_true')
parser.add_argument("-l"  , "--leg"       , dest = "leg"       ,  help = "legend labels"    , default = ''                        )

options = parser.parse_args()
if not options.input:   
  parser.error('Input filename not given')

import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis, TPad, TGraphErrors
from   array import array
from   math  import sqrt, isnan
from   copy  import deepcopy as dc
import numpy as np

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

namefiles = options.input.split(',')
nfiles   = len(namefiles)
files    = []

print 'number of input files is ' + str(nfiles)

for i in range(0, nfiles):
  print 'opening file ' + str(i) + ': ' + namefiles[i]
  files.append(TFile.Open(namefiles[i]   , 'r') )


pt_bins  = [  0, 15, 16, 17, 18, 20, 22, 25, 30, 40, 50, 60, 80, 100] 
eta_bins = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
lxy_bins  = [  0, 1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 30, 100] 

colorlist = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kAzure+1, ROOT.kViolet, ROOT.kOrange, ROOT.kBlue, ROOT.kGray, ROOT.kMagenta+3, ROOT.kOrange+2, ROOT.kYellow]

sel = ''

def doHisto(file, var, thecolor, i, lmin, lmax, sel):

  tt = file.Get('tree')
  if var[9] == -1:
    hnum = ROOT.TH1F('hnum', 'hnum', var[3], var[4][0], var[4][1])
    hden = ROOT.TH1F('hden', 'hden', var[3], var[4][0], var[4][1])
  else:  
    print len(var[9])-1
    print np.asarray(var[9])
    hnum = ROOT.TH1F('hnum', 'hnum', len(var[9])-1, array('d',var[9]))
    hden = ROOT.TH1F('hden', 'hden', len(var[9])-1, array('d',var[9]))
  tt.Draw('tau_gen_%s>>+hnum'%(var[0]), 'tau_gen_lxy>-100 && tau_gen_vis_pt>15 && abs(tau_gen_vis_eta)<2.1 && tau_gen_vis_pt>20 && tau_gen_vis_pt<1000 && tau_gen_lxy>%s && tau_gen_lxy<%s && tau_l1_iso==1 && tau_l1_pt>20 '%(lmin,lmax) + sel, 'goff')
  tt.Draw('tau_gen_%s>>+hden'%(var[0]), 'tau_gen_lxy>-100 && tau_gen_vis_pt>15 && abs(tau_gen_vis_eta)<2.1 && tau_gen_vis_pt>20 && tau_gen_vis_pt<1000 && tau_gen_lxy>%s && tau_gen_lxy<%s'%(lmin,lmax) + sel, 'goff')
  pEff = dc(TEfficiency(hnum,hden));
    
  pEff.SetLineColor  (thecolor)
  pEff.SetMarkerColor(thecolor)
  pEff.SetMarkerStyle(8  )
  pEff.SetMarkerSize(0.8)

  pEff.SetTitle(";" + var[1] + ";" + var[2])
  return pEff

ytitle = 'efficiency'
variables = [
#  numerator          # x axis title            # y title   # nbins    # x range      # y range      # pdf name        # legend position         #y range ratio          
 ('lxy'               , 'L_{xy}'               , ytitle,   30   , (  0   , 100   ), (0.  , 1.4),  'efficiency_lxy',  (0.2 , 0.58, 0.7, 0.88),  (0.501, 1.05 ) , lxy_bins), 
 ('vis_pt'            , 'vis p_{T} [GeV]'      , ytitle,   100  , (  0.  , 100  ), (0.   , 1.4),  'efficiency_pt',   (0.2 , 0.58, 0.7, 0.88),  (0.5  , 1.05  ), pt_bins),
 ('vis_eta'           , '#eta'                 , ytitle,    50  , ( -2.2 , 2.2  ), (0.5  , 1.4),  'efficiency_eta',  (0.2 , 0.58, 0.7, 0.88),  (0.5  , 1.05  ), -1),
 ('cosxy'             , 'vis cosine'           , ytitle,    50  , ( -1.  ,  1.  ), (0.   , 1.4),  'efficiency_cos',  (0.2 , 0.58, 0.7, 0.88),  (0.5  , 1.05  ), -1),
 ('momct'             , 'ct [cm]'              , ytitle,    50  , ( 0   , 0.50  ), (0.   , 1.4),  'efficiency_ct',  (0.2 , 0.58, 0.72, 0.88),  (0.5  , 1.05  ), -1),
] 

import pdb

c2 = TCanvas('c2', 'c2', 600,600)
c1 = TCanvas('c1', 'c1', 600,600)
if options.diff:
  stackPad = ROOT.TPad('stackPad', 'stackPad', 0.,  .25, 1., 1.  , 0, 0)  
  ratioPad = ROOT.TPad('ratioPad', 'ratioPad', 0., 0. , 1.,  .3, 0, 0)  
else:
  stackPad = ROOT.TPad('stackPad', 'stackPad', 0,  0. , 1., 1.  , 0, 0)  
  ratioPad = ROOT.TPad('ratioPad', 'ratioPad', 0., 0. , 0., 0.  , 0, 0)  


for var in variables:

  c1.cd()
  stackPad.Draw()
  ratioPad.Draw()
  
  l = TLegend(var[7][0], var[7][2], var[7][1], var[7][3])
  l.SetBorderSize(0)
  l.SetTextSize(0.026)
  
  eff_list = [] 
  sel = ''
#   sel = '&& (tau_gen_decaymode>0 && tau_gen_decaymode<=4)'
  for i, ifile in enumerate(files):
    c2.cd() 
    eff_list.append(dc(doHisto(ifile , var, colorlist[i], i, 0, 100000, sel)))

  c1.cd()
  stackPad.cd()
  for i,k in enumerate(eff_list):
    if (i == 0):
      k.Draw('AP')
      c1.Update()
      c1.Modified()

      k.GetPaintedGraph().GetXaxis().SetLimits(var[4][0], var[4][1])
      k.GetPaintedGraph().GetHistogram().SetMinimum(var[5][0])             
      k.GetPaintedGraph().GetHistogram().SetMaximum(var[5][1])        
  
      k.GetPaintedGraph().GetXaxis().SetLabelSize(0.04)
      k.GetPaintedGraph().GetYaxis().SetLabelSize(0.04)   
      k.GetPaintedGraph().GetXaxis().SetTitleSize(0.04)
      k.GetPaintedGraph().GetYaxis().SetTitleSize(0.04)
      k.GetPaintedGraph().GetYaxis().SetTitleOffset(1.5)
      k.GetPaintedGraph().GetXaxis().SetTitleOffset(1.2)

    else:
      k.Draw('P same')

    if options.leg:
      l.AddEntry(k , options.leg.split(',')[i]  , "pel")

    
  l.Draw()
#   if 'eff_tkIP' in var[0]:
#     gPad.SetLogx()  
#   else:
#     gPad.SetLogx(False)  


  gPad.SetLogx()  
  gPad.SetGridx(True)
  gPad.SetGridy(True)
  c1.SaveAs("" +  var[6] + "_isotau_hnl_logL.pdf")



