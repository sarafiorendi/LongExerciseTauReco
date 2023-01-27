import ROOT
from   ROOT  import TFile, TTree, gDirectory, TH1F, TCanvas, TLegend, TEfficiency, gPad, gStyle, TGaxis, TPad, TGraphErrors
import pdb
import numpy as np

gStyle.SetOptStat('emr')
gStyle.SetTitleAlign(23)
gStyle.SetPadLeftMargin(0.16)
gStyle.SetPadBottomMargin(0.16)
TGaxis.SetMaxDigits(3)

ROOT.gROOT.SetBatch(True)


thefile = TFile.Open('ntuples/mutau_gmsb_200Gev_100mm_muAndTau_v18_l1tau_l1mu18.root'   , 'r')
thetree = thefile.Get('tree')

min_tau_gen_pt = 20
min_l1_mu_pt = 18
min_l1_singlemu_pt = 22
min_l1_tau_pt = 24

min_dxy = 0;
max_dxy = 0.01;
step_dxy = 0.001;

min_gen_lxy_vec = [0,  2,  10]
max_gen_lxy_vec = [2, 10, 1000]
# min_gen_lxy_vec = [0,  1,  3,  10, 100]
# max_gen_lxy_vec = [1,  3, 10, 100, 1000]

range_tau_pt = [i for i in range(24,33)]
range_mu_pt = [i for i in range(18,28)]
range_mu_dxy = [i for i in np.arange(min_dxy, max_dxy, step_dxy)]
range_singlemu_pt = [i for i in range(22,40)]
range_singlemu_eta = [i for i in np.arange(1.4, 2.1, 0.1)]


has_gen_mu  = '(abs(mu_gen_vis_eta)<2.1 && mu_gen_vis_pt>%s ) ' %(min_tau_gen_pt)
has_gen_tau = '(abs(tau_gen_vis_eta)<2.1 && tau_gen_vis_pt>%s)' %(min_tau_gen_pt) 

for min_gen_lxy,max_gen_lxy in zip(min_gen_lxy_vec, max_gen_lxy_vec):

  h_2d_eff_singlemu =     ROOT.TH2F('h_2d_eff_singlemu',     ';p_{T}(#mu) [GeV]; muon d_{xy}', len(range_singlemu_pt)-1, range_singlemu_pt[0], range_singlemu_pt[-1], len(range_mu_dxy)-1, range_mu_dxy[0], range_mu_dxy[-1]  )
  h_2d_eta_eff_singlemu = ROOT.TH2F('h_2d_eta_eff_singlemu', ';p_{T}(#mu) [GeV]; |#eta(#mu)|', len(range_singlemu_pt)-1, range_singlemu_pt[0], range_singlemu_pt[-1], len(range_singlemu_eta)-1, range_singlemu_eta[0], range_singlemu_eta[-1]  )

  h_2d_eff_pt =    ROOT.TH2F('h_2d_eff_pt',    ';p_{T}(#tau) [GeV];p_{T}(#mu) [GeV]', len(range_tau_pt)-1, range_tau_pt[0], range_tau_pt[-1], len(range_mu_pt)-1, range_mu_pt[0], range_mu_pt[-1]  )
  h_2d_eff_pt_mu = ROOT.TH2F('h_2d_eff_pt_mu', ';p_{T}(#tau) [GeV];p_{T}(#mu) [GeV]', len(range_tau_pt)-1, range_tau_pt[0], range_tau_pt[-1], len(range_mu_pt)-1, range_mu_pt[0], range_mu_pt[-1]  )
  
  dxy_gen_tau = '(abs(tau_gen_lxy) > %s && abs(tau_gen_lxy) < %s) && (abs(mu_gen_lxy) > %s && abs(mu_gen_lxy) < %s)'%(min_gen_lxy, max_gen_lxy, min_gen_lxy, max_gen_lxy) 
  denominator_string = has_gen_mu + ' && ' + has_gen_tau + ' && ' + dxy_gen_tau 
  nden = thetree.Draw('tau_gen_pt>>hden', denominator_string, 'goff')
  print (min_gen_lxy,max_gen_lxy, ' ->', nden)
  
  l1_string_mutau = '(tau_l1mu18_pt > %s && abs(tau_l1mu18_eta)<2.1) && (tau_l1_pt > %s && abs(tau_l1_eta)<2.1) '%(min_l1_mu_pt, min_l1_tau_pt)
  l1_string_mu = '(tau_l1mu_pt > %s && abs(tau_l1mu_eta)<2.1) '%(min_l1_singlemu_pt)


  ## single muon trigger
  for cut_mu_dxy in range_mu_dxy:
    hlt_cut_dxy_string = '(abs(tau_glb_dxy) > %s )'%cut_mu_dxy
    for hlt_mu_pt in range_singlemu_pt[:-1]:
  
      hlt_mu_string = '(tau_glb_pt > %s && abs(tau_glb_eta)<2.1 ) '%hlt_mu_pt
      num_string = l1_string_mu + ' && ' + hlt_mu_string + ' && ' + hlt_cut_dxy_string + ' && ' + denominator_string
      nnum = thetree.Draw('tau_gen_pt>>hnum', num_string, 'goff')

#       print (min_gen_lxy,max_gen_lxy, ' + ' , cut_mu_dxy , ' -> num = ', nnum, ' -> eff = ', float(nnum)/nden)
  
#       print ('eff for dxy %s and mu %s: %s'%(cut_mu_dxy, hlt_mu_pt, float(nnum)/nden))
      h_2d_eff_singlemu.Fill(hlt_mu_pt+0.01, cut_mu_dxy+0.0001, float(nnum)/nden)
  #     h_2d_eff_pt.Fill(cut_mu_dxy+0.0001, hlt_mu_pt+0.01, float(nnum)/nden)
      
  c1 = TCanvas()
  ROOT.gStyle.SetPaintTextFormat("4.3f ");
  h_2d_eff_singlemu.SetStats(0)
  h_2d_eff_singlemu.GetZaxis().SetRangeUser(0.4, .8)
  h_2d_eff_singlemu.SetMinimum(0.4)
  h_2d_eff_singlemu.SetMaximum(0.8)
  h_2d_eff_singlemu.SetTitle('GEN Lxy %s-%s cm'%(min_gen_lxy, max_gen_lxy))
  h_2d_eff_singlemu.Draw('colz text')
  c1.SaveAs('2d_eff_dxy_pt_singlemu_genLxy%s_%s.pdf'%(min_gen_lxy, max_gen_lxy))    

  ##  **********************************************
  for cut_mu_eta in range_singlemu_eta:
    hlt_cut_eta_string = '(abs(tau_glb_eta) < %s )'%cut_mu_eta
    for hlt_mu_pt in range_singlemu_pt[:-1]:
  
      hlt_mu_string = '(tau_glb_pt > %s && abs(tau_glb_eta)<2.1 ) '%hlt_mu_pt
      num_string = l1_string_mu + ' && ' + hlt_mu_string + ' && ' + hlt_cut_eta_string + ' && ' + denominator_string
      nnum = thetree.Draw('tau_gen_pt>>hnum', num_string, 'goff')
      h_2d_eta_eff_singlemu.Fill(hlt_mu_pt+0.01, cut_mu_eta-0.0001, float(nnum)/nden)
  #     h_2d_eff_pt.Fill(cut_mu_dxy+0.0001, hlt_mu_pt+0.01, float(nnum)/nden)
      
  c1 = TCanvas()
  ROOT.gStyle.SetPaintTextFormat("4.3f ");
  h_2d_eta_eff_singlemu.SetStats(0)
  h_2d_eta_eff_singlemu.GetZaxis().SetRangeUser(0.3, .8)
  h_2d_eta_eff_singlemu.SetMinimum(0.3)
  h_2d_eta_eff_singlemu.SetMaximum(0.8)
  h_2d_eta_eff_singlemu.SetTitle('GEN Lxy %s-%s cm'%(min_gen_lxy, max_gen_lxy))
  h_2d_eta_eff_singlemu.Draw('colz text')
  c1.SaveAs('2d_eff_eta_pt_singlemu_genLxy%s_%s.pdf'%(min_gen_lxy, max_gen_lxy))    
  


  ##  **********************************************
  #  mu + tau trigger + dxy 0.005 for muon
  for hlt_tau_pt in range_tau_pt[:-1]:
    hlt_tau_string = '(tau_hltPFdispltau_pt > %s && abs(tau_hltPFdispltau_eta)<2.1 && tau_hltPFdispltau_passChargedIso > 0 && tau_hltPFdispltau_passFilters > 0 )'%hlt_tau_pt
    for hlt_mu_pt in range_mu_pt[:-1]:
  
      hlt_mu_string = '(tau_glb_pt > %s && abs(tau_glb_eta)<2.1 && abs(tau_glb_dxy) > 0.005) '%hlt_mu_pt
      num_string = l1_string_mutau + ' && ' + hlt_tau_string + ' && ' + hlt_mu_string + ' && ' + denominator_string
      nnum = thetree.Draw('tau_gen_pt>>hnum', num_string, 'goff')
      print ('n ev passing for pt tau %s and mu %s: %s'%(hlt_tau_pt, hlt_mu_pt, float(nnum)))
  
      h_2d_eff_pt.Fill(hlt_tau_pt+0.01, hlt_mu_pt+0.01, float(nnum)/nden)
      
  c1 = TCanvas()
  h_2d_eff_pt.SetStats(0)
  h_2d_eff_pt.GetZaxis().SetRangeUser(0.07, .35)
  h_2d_eff_pt.SetMinimum(0.07)
  h_2d_eff_pt.SetMaximum(0.35)
  h_2d_eff_pt.SetTitle('GEN Lxy %s-%s cm'%(min_gen_lxy, max_gen_lxy))
  h_2d_eff_pt.Draw('colz text')
  c1.SaveAs('2d_eff_pt_mu0p005_genLxy%s_%s.pdf'%(min_gen_lxy, max_gen_lxy))    
  
  
  ##  **********************************************
  #  mu + tau trigger
  for hlt_tau_pt in range_tau_pt[:-1]:
    hlt_tau_string = '(tau_hltPFdispltau_pt > %s && abs(tau_hltPFdispltau_eta)<2.1 && tau_hltPFdispltau_passChargedIso > 0 && tau_hltPFdispltau_passFilters > 0 )'%hlt_tau_pt
    for hlt_mu_pt in range_mu_pt[:-1]:
  
      hlt_mu_string = '(tau_glb_pt > %s && abs(tau_glb_eta)<2.1) '%hlt_mu_pt
      num_string = l1_string_mutau + ' && ' + hlt_tau_string + ' && ' + hlt_mu_string + ' && ' + denominator_string
      nnum = thetree.Draw('tau_gen_pt>>hnum', num_string, 'goff')
      h_2d_eff_pt_mu.Fill(hlt_tau_pt+0.01, hlt_mu_pt+0.01, float(nnum)/nden)
      
  c1 = TCanvas()
  h_2d_eff_pt_mu.SetStats(0)
  h_2d_eff_pt_mu.GetZaxis().SetRangeUser(0.07, .3)
  h_2d_eff_pt_mu.SetMinimum(0.07)
  h_2d_eff_pt_mu.SetMaximum(0.3)
#   h_2d_eff_pt_mu.GetZaxis().SetRangeUser(0, .35)
#   h_2d_eff_pt_mu.SetMinimum(0)
#   h_2d_eff_pt_mu.SetMaximum(0.35)
  h_2d_eff_pt_mu.SetTitle('GEN Lxy %s-%s cm'%(min_gen_lxy, max_gen_lxy))
  h_2d_eff_pt_mu.Draw('colz text')
  c1.SaveAs('2d_eff_pt_genLxy%s_%s.pdf'%(min_gen_lxy, max_gen_lxy))    
  
  