'''
In this file you can define the branch to be included in the flat ntuple.
You can find some basic quantities here, expand with more specific observables,
such as isolation etc...
'''
branches = [
    'run',
    'lumi',
    'event',
    'trueNI',
    'nvtx',
    'gen_met',
    'nano_met',
    'nano_sumet',
    'miniaod_met_pt',

    'pass_met_hlt',
    'pass_lep_hlt',
    'pass_PFMET120_PFMHT120_IDTight', 
    'pass_PFMET130_PFMHT130_IDTight', 
    'pass_PFMET140_PFMHT140_IDTight', 
    'pass_PFMETNoMu110_PFMHTNoMu110_IDTight', 
    'pass_PFMETNoMu120_PFMHTNoMu120_IDTight', 
    'pass_PFMETNoMu130_PFMHTNoMu130_IDTight', 
    'pass_PFMETNoMu140_PFMHTNoMu140_IDTight', 
    'pass_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1',    
    'pass_MET105_IsoTrk50', 
    'pass_MET120_IsoTrk50',
    'pass_Ele30_WPTight_Gsf',
    'pass_Photon200',
    'pass_Photon20',
    'pass_Photon30_R9Id90_CaloIdL_IsoL_DisplacedIdL',
    'pass_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1',
    'pass_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1',

    'tau_reco_mass',
    'tau_reco_pt',
    'tau_reco_eta',
    'tau_reco_phi',
    'tau_reco_charge',
    'tau_reco_decaymode',
    'tau_reco_npix',

    'hlt_tau_mass'    ,
    'hlt_tau_pt'      ,
    'hlt_tau_eta'     ,
    'hlt_tau_phi'     ,
    'hlt_tau_charge'  ,
    'hlt_tau_decaymode'  ,

    'hlt_tau_leadChargedCandPdgId',
    'hlt_tau_leadChargedCandPt'   ,
    'hlt_tau_leadNeutralCandPdgId',
    'hlt_tau_leadNeutralCandPt'   ,
    'hlt_tau_leadPFCandPdgId'     ,
    'hlt_tau_leadPFCandPt'        ,
    'hlt_tau_maxHCALPFClusterEt',
    'hlt_tau_nChargedHad',
    'hlt_tau_nGamma'     ,
    'hlt_tau_sumPtCharged',
    'hlt_tau_sumPtNeutral',

    'hlt_tau_dxy'  ,
    'hlt_tau_dxyerr'  ,
    'hlt_tau_ip3d'  ,
    'hlt_tau_ip3derr'  ,

    'hlt_tau_passChargedIso'  ,
    'hlt_tau_passRelChargedIso'  ,
    'hlt_tau_passAbsChargedIso'  ,
    'hlt_tau_isoval',
    'hlt_tau_passFilters',

    'tau_hlttaumatchmom_pt'       ,
    'tau_hlttaumatchmom_eta'      ,
    'tau_hlttaumatchmom_phi'      ,
    'tau_hlttaumatchmom_charge'   ,
    'tau_hlttaumatchmom_decaymode',

    'tau_offtk_pt'     ,
    'tau_offtk_eta'    ,
    'tau_offtk_phi'    ,
    'tau_offtk_charge' ,
    'tau_offtk_algo',
    'tau_reco_tau_mass'     ,
    'tau_reco_tau_pt'       ,
    'tau_reco_tau_eta'      ,
    'tau_reco_tau_phi'      ,
    'tau_reco_tau_charge'   ,
    'tau_reco_tau_decaymode',

#     'tau_hltTk_dxy' ,
    'l1_tau_pt',
    'l1_tau_eta',
    'l1_tau_phi',
    'l1_tau_charge',
    'l1_tau_iso',

    'l1_mu_pt'    ,
    'l1_mu_eta'   ,
    'l1_mu_phi'   ,
    'l1_mu_charge',
    'l1_mu18_pt'    ,
    'l1_mu18_eta'   ,
    'l1_mu18_phi'   ,
    'l1_mu18_charge',
    'l3mu_pt'    ,
    'l3mu_eta'   ,
    'l3mu_phi'   ,
    'l3mu_charge',
    'tau_isomu_pt'    ,
    'tau_isomu_eta'   ,
    'tau_isomu_phi'   ,
    'tau_isomu_charge',

    'tau_l2disp_pt'    ,
    'tau_l2disp_eta'   ,
    'tau_l2disp_phi'   ,
    'tau_l2disp_charge',
    'tau_cascade_pt'    ,
    'tau_cascade_eta'   ,
    'tau_cascade_phi'   ,
    'tau_cascade_charge',
    'tau_cascade_dxy',
    'hlt_mu_pt'    ,
    'hlt_mu_eta'   ,
    'hlt_mu_phi'   ,
    'hlt_mu_charge',
    'hlt_mu_dxy'   ,

    'gen_mu_ndau',
    'gen_mu_pt',
    'gen_mu_eta',
    'gen_mu_phi',
    'gen_mu_charge',
    'gen_mu_lxy', 
    'gen_mu_lxy_mu',
    'gen_mu_dxy',  
    'gen_mu_vx',  
    'gen_mu_vy',  
    'gen_mu_vz',  
    'gen_mu_cosxy',  
    'gen_mu_mom_mass', 
    'gen_mu_mom_pt', 
    'gen_mu_mom_eta', 
    'gen_mu_mom_phi', 
    'gen_mu_momct2d',  
    'gen_mu_vis_mass',
    'gen_mu_vis_pt',
    'gen_mu_vis_eta',
    'gen_mu_vis_phi',

    'gen_tau_vis_mass',
    'gen_tau_vis_pt',
    'gen_tau_vis_eta',
    'gen_tau_vis_phi',
    'gen_tau_dau_pt' ,
    'gen_tau_dau_eta',
    'gen_tau_dau_phi',
    'gen_tau_lxy', 
    'gen_tau_dxy',  
    'gen_tau_visdxy',  
    'gen_tau_vx',  
    'gen_tau_vy',  
    'gen_tau_vz',  
    'gen_tau_cosxy',  
    'gen_tau_momct',  
    'gen_tau_momct2d',  
    'gen_tau_mom_mass', 
    'gen_tau_mom_pt', 
    'gen_tau_mom_eta', 
    'gen_tau_mom_phi', 

    'gen_tau_ndau',
    'gen_tau_pt',
    'gen_tau_eta',
    'gen_tau_phi',
    'gen_tau_charge',
    'gen_tau_decaymode',

    'PV_x',
    'PV_y',
    'PV_z',
]
