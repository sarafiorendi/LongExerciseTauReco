import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis, TPad, TGraphErrors
import pdb

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

ROOT.gROOT.SetBatch(True)


thefile = TFile.Open('ntuples/ditau_gmsb_200Gev_100mm_forDiTau.root'   , 'r')
thetree = thefile.Get('tree')

min_tau_gen_pt = 20
min_l1_tau_pt = 32

min_gen_lxy_vec = [0,  2,  10]
max_gen_lxy_vec = [2, 10, 1000]

range_tau_pt = [i for i in range(32,45)]

for min_gen_lxy,max_gen_lxy in zip(min_gen_lxy_vec, max_gen_lxy_vec):

  h_2d_eff_pt = ROOT.TH2F('h_2d_eff_pt', ';p_{T}(#tau_{1}) [GeV] ;p_{T}(#tau_{2}) [GeV]', len(range_tau_pt)-1, range_tau_pt[0], range_tau_pt[-1], len(range_tau_pt)-1, range_tau_pt[0], range_tau_pt[-1]  )
  
  
  has_gen_tau1 = '(abs(tau_gen_vis_eta_t1)<2.1 && tau_gen_vis_pt_t1>%s)' %(min_tau_gen_pt) 
  has_gen_tau2 = '(abs(tau_gen_vis_eta_t2)<2.1 && tau_gen_vis_pt_t2>%s)' %(min_tau_gen_pt) 
  dxy_gen_tau = '(abs(tau_gen_lxy_t1) > %s && abs(tau_gen_lxy_t1) < %s) && (abs(tau_gen_lxy_t2) > %s && abs(tau_gen_lxy_t2) < %s)'%(min_gen_lxy, max_gen_lxy, min_gen_lxy, max_gen_lxy) 
  denominator_string = has_gen_tau1 + ' && ' + has_gen_tau2 + ' && ' + dxy_gen_tau
  nden = thetree.Draw('tau_gen_pt_t1>>hden', denominator_string, 'goff')
  
  l1_string = '(tau_l1_pt_t1 > %s && abs(tau_l1_eta_t1)<2.1) && (tau_l1_pt_t2 > %s && abs(tau_l1_eta_t2)<2.1) '%(min_l1_tau_pt, min_l1_tau_pt)
  
  for hlt_tau1_pt in range_tau_pt[:-1]:
    hlt_tau1_string = '(tau_hltPFdispltau_pt_t1 > %s && abs(tau_hltPFdispltau_eta_t1)<2.1 && tau_hltPFdispltau_passChargedIso_t1 > 0 && tau_hltPFdispltau_passFilters_t1 > 0 )'%hlt_tau1_pt
    for hlt_tau2_pt in range_tau_pt[:-1]:
  
      hlt_tau2_string = '(tau_hltPFdispltau_pt_t2 > %s && abs(tau_hltPFdispltau_eta_t2)<2.1 && tau_hltPFdispltau_passChargedIso_t2 > 0 && tau_hltPFdispltau_passFilters_t2 > 0 ) '%hlt_tau2_pt
      
      num_string = l1_string + ' && ' + hlt_tau1_string + ' && ' + hlt_tau2_string + ' && ' + denominator_string
      nnum = thetree.Draw('tau_gen_pt_t1>>hnum', num_string, 'goff')
  
      print ('eff for tau %s and mu %s: %s (n passing = %s)'%(hlt_tau1_pt, hlt_tau2_pt, float(nnum)/nden, nnum))
      h_2d_eff_pt.Fill(hlt_tau1_pt+0.01, hlt_tau2_pt+0.01, float(nnum)/nden)
      
  c1 = TCanvas()
  h_2d_eff_pt.SetStats(0)
  h_2d_eff_pt.GetZaxis().SetRangeUser(0.02, .14)
  h_2d_eff_pt.SetMinimum(0.02)
  h_2d_eff_pt.SetMaximum(0.14)

  gStyle.SetPaintTextFormat("0.3f");
  h_2d_eff_pt.Draw('colz text')
  c1.SaveAs('2d_eff_pt_ditau_genLxy%s_%s.pdf'%(min_gen_lxy, max_gen_lxy)   ) 




h_1d_eff_lxy = ROOT.TH1F('h_1d_eff_lxy', ';L_{xy} [cm] ;#varepsilon_{evt}', len(range_tau_dxy)-1, range_tau_dxy[0], range_tau_dxy[-1])






