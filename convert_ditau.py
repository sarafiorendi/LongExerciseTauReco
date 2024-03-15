import os, sys, pdb
import numpy as np
import pandas
from copy import deepcopy
from pdb import set_trace
from array import array
from math  import sqrt 
import uproot

# sig_list = [ 'tau_ditau_tuple_gmsb_m100_ctau100mm_addL2isoVals_eos_2_withDispl_compareToPromptTrigger.root:tree']
sig_list = [ 'tau_ditau_tuple_gmsb_m100_ctau100mm_norehlt_summer22_fnal_1_withDispl_newsample.root:tree']

sig_uproot = uproot.concatenate(sig_list,  library="np")

sig = pandas.DataFrame(sig_uproot)
todrop = [
'tau_reco',
'tau_hltfo',
'tau_hlt_',
'tau_hps',
'tau_hltPFjet',
'tau_hltPFtau_',
'tau_hltPFp',
'tau_hltPF_',
'jet_',
'tau_offtk',
'tau_hltMuMTk',
'tau_hltIt',
'tau_hltTk',
# 'Reg',
'_0_',
'_1_',
'_2_',
'l1_emu',
'mu_l3',
'2fo_',
]

for istring in todrop:
  sig = sig[sig.columns.drop(list(sig.filter(regex=istring)))]

sig = sig.sort_values(['tau_gen_vis_pt'],ascending=False)

sig['counter'] = sig.groupby(['event']).cumcount() + 1 
tau_one = sig[sig.counter == 1]
tau_two = sig[sig.counter == 2]

tokeep = [
'run',
'lumi',
'event',
'trueNI',
'nvtx',
'PV_x',
'PV_y',
'PV_z',
'pass_met_hlt',
'pass_lep_hlt',
'pass_PFMET120_PFMHT120_IDTight', 
'pass_PFMET130_PFMHT130_IDTight', 
'pass_PFMET140_PFMHT140_IDTight', 
'pass_PFMETNoMu110_PFMHTNoMu110_IDTight', 
'pass_PFMETNoMu120_PFMHTNoMu120_IDTight', 
'pass_PFMETNoMu130_PFMHTNoMu130_IDTight', 
'pass_PFMETNoMu140_PFMHTNoMu140_IDTight', 
'pass_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg', 
'pass_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1',
'pass_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100', 
'pass_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110', 
'pass_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120', 
'pass_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130', 
'pass_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140', 
'pass_MET105_IsoTrk50', 
'pass_MET120_IsoTrk50',
]

tau_one = tau_one.rename(columns={col: col+'_t1' 
                          for col in tau_one.columns if col not in tokeep})
                          
tau_two = tau_two.rename(columns={col: col+'_t2' 
                          for col in tau_two.columns if col not in tokeep})                          
                          
                          
all =  pandas.merge(tau_one,tau_two, on='event')        

data = {key: all[key].values for key in all.columns}
import ROOT
rdf = ROOT.RDF.MakeNumpyDataFrame(data)
rdf.Snapshot("tree", sig_list[0].replace('.root:tree','_forDiTau_addhltinfo.root'))


                 
