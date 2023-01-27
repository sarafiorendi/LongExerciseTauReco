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

    'l1_tau_pt',
    'l1_tau_eta',
    'l1_tau_phi',
    'l1_tau_charge',
    'l1_tau_iso',

     ## egamma variables
    'l1_ele_pt'    ,
    'l1_ele_eta'   ,
    'l1_ele_phi'   ,
    'l1_ele_charge',
    'hlt_ele_pt'    ,
    'hlt_ele_eta'   ,
    'hlt_ele_phi'   ,
    'hlt_ele_charge',
    'filter1_ele_pt'    ,
    'filter1_ele_eta'   ,
    'filter1_ele_phi'   ,
    'filter1_ele_charge',

    'filter2_ele_pt'    ,
    'filter2_ele_eta'   ,
    'filter2_ele_phi'   ,
    'filter2_ele_charge',

    'filter3_ele_pt'    ,
    'filter3_ele_eta'   ,
    'filter3_ele_phi'   ,
    'filter3_ele_charge',

    'filter4_ele_pt'    ,
    'filter4_ele_eta'   ,
    'filter4_ele_phi'   ,
    'filter4_ele_charge',

    'filter5_ele_pt'    ,
    'filter5_ele_eta'   ,
    'filter5_ele_phi'   ,
    'filter5_ele_charge',

    'filter6_ele_pt'    ,
    'filter6_ele_eta'   ,
    'filter6_ele_phi'   ,
    'filter6_ele_charge',
    'filter7_ele_pt'    ,
    'filter7_ele_eta'   ,
    'filter7_ele_phi'   ,
    'filter7_ele_charge',
    'filter8_ele_pt'    ,
    'filter8_ele_eta'   ,
    'filter8_ele_phi'   ,
    'filter8_ele_charge',

    'gen_ele_ndau',
    'gen_ele_pt',
    'gen_ele_eta',
    'gen_ele_phi',
    'gen_ele_charge',
    'gen_ele_lxy', 
    'gen_ele_dxy',  
    'gen_ele_vx',  
    'gen_ele_vy',  
    'gen_ele_vz',  
    'gen_ele_cosxy',  
    'gen_ele_mom_mass', 
    'gen_ele_mom_pt', 
    'gen_ele_mom_eta', 
    'gen_ele_mom_phi', 
    'gen_ele_momct2d',  
    'gen_ele_vis_mass',
    'gen_ele_vis_pt',
    'gen_ele_vis_eta',
    'gen_ele_vis_phi',

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
]