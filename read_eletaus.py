'''
Loops on the events and operates the matching between reconstructed and generated taus.
It produces two flat ntuples:
    - one with an entry for each gen tau (useful for efficiencies)
    - one with an entry for each reconstructed tau (useful for fake studies)
'''
import ROOT, sys, os, pdb, argparse
import numpy as np
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from treeVariablesHLT_ele import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer , genDecayModeGEANT, isAncestor, isGenLepTau# utility functions
from math import sqrt, pow
from copy import deepcopy as dc
from itertools import product as itertools_product

ROOT.gSystem.Load('libRecoTauTagRecoTau')
from ROOT import EcalClusterTools

# from RecoEcal.EgammaCoreTools/interface/EcalClusterTools.h
# EcalClusterTools


parser = argparse.ArgumentParser(description="Convert MiniAOD to flat ntuples!")
parser.add_argument(
	"--sample",
# 	choices=['ZTT','SMS', 'test'],
	required=True,
	help='Specify the sample you want to flatten')

parser.add_argument(
	"--file",
	required=False,
	default=0,
	type=int,
	help='Specify the file number in the list from das')

args = parser.parse_args()
ifile = args.file
sample = args.sample

##########################################################################################
feat_list = ['lxy', 'dxy', 'visdxy', 'cosxy', 'momct', 'momct2d', 'mom_mass', 'pi_lxy', 'pi_cosxy', 'lxy_ele']
pftau_feat_list = ['leadChargedCandPt', 'leadChargedCandPdgId', 'leadCandPt', 'leadCandPdgId', 'maxHCALPFClusterEt', 'nChargedHad', 'nGamma', 'sum_pt_charged', 'sum_pt_neutral', 'dxy', 'dxyerr', 'ip3d', 'ip3derr', 'isoVal', 'passChargedIso', 'passRelChargedIso', 'passAbsChargedIso', 'passFilters']
ele_feat_list = ['dxy', 'r9', 'clsh', 'ecaliso', 'hcaliso', 'trkiso', 'passFilters', 'sMin', 'sMaj', 'passSingleEGFilters','passSingleL1Filter']
l1_feat_list = ['passSingleL1Filter', 'passEleTauL1Filter', 'passGsfEleL1Filter']
c_const = 299.792458
dR_cone = 0.5
dPt_cone = 0.2
good_gen_status = 2
mom_pdgId = [1000015]
if 'HNL' in sample:  mom_pdgId = [9900012]
if 'HNL' in sample and 'Dirac' in sample:  mom_pdgId = [9990012]
if 'gmsb' in sample:  
    mom_pdgId = [2000015, 1000015]
    good_gen_status = 2 #@ was 8 in previous samples
if 'DY' in sample:  mom_pdgId = [23]

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_ele_tuple_{}_{}.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(list(zip(branches, [-9999.]*len(branches)))) # initialise all branches to unphysical -99       
print((outfile_gen.GetName())) 

##########################################################################################
# Get ahold of the events
file_min = 1
file_max = 200
f = open('%s_filelist.txt'%args.sample)
infile = f.readlines()[file_min:file_max]

os.environ['X509_USER_PROXY'] = '/afs/cern.ch/user/f/fiorendi/x509up_u58808' 
redirector = 'root://cms-xrd-global.cern.ch//'
if 'fnal' in sample:  redirector = 'root://cmseos.fnal.gov//'
if 'eos' in sample:   redirector = ''

# redirector = ''
events = Events(redirector+infile.strip())
# events = Events([
#                  '/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_6_2/src/HLTrigger/Configuration/test/hltOutput.root'
#                  '/eos/cms/store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC_Summer21/crab_ntuples_mutau_singlephoton_v18_gmsb_M200_100mm_summer21/220627_083854/0000/outputHLT_111.root',
#                  '/eos/cms/store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC_Summer21/crab_ntuples_mutau_singlephoton_v18_gmsb_M200_100mm_summer21/220628_134709/0000/outputHLT_69.root'
#                 ])

# print('using a local file!!!!!!!!!!!!!!!!!!!!!')

print(infile)
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files


##########################################################################################
def findMatchToGen(gen_taus, hlt_taus, hlt_tau, dR_cone_ = dR_cone): 
 
    gen_taus_match = gen_taus 
    for gg in gen_taus_match: 
        setattr(gg, hlt_tau, None)
        bestcom = bestMatch(gg.visp4, hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone_ :
            setattr(gg,hlt_tau,bestcom[0])
    return gen_taus_match    


def findMatchToGenMother(gen_taus, hlt_taus, hlt_tau, dR_cone_ = dR_cone): 
 
    gen_taus_match = gen_taus 
    for gg in gen_taus_match: 
        setattr(gg, hlt_tau, None)
#         pdb.set_trace()
        bestcom = bestMatch(gg.bestmom.p4(), hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone_ :
            setattr(gg,hlt_tau,bestcom[0])
    return gen_taus_match    




def printForDebug(cc, gg):
    print(('    pt: %.2f'%( cc.pt()), \
          '\t eta: %.2f'%( cc.eta()), \
          '\t phi: %.2f'%( cc.phi()), \
          '\t pdgID: %d'%( cc.pdgId()), \
          '\t dpT/pT = %.2f,  dR = %.2f '%(abs(cc.pt() - gg.vispt())/gg.vispt(), deltaR(cc.eta(),cc.phi(),gg.viseta(), gg.visphi()))))

def printGEN(gg):
    print(('GEN pt: %.2f'%(gg.vispt(),), \
          '\t eta: %.2f'%(gg.viseta()), \
          '\t phi: %.2f'%(gg.visphi()), '\n'))

   
##########################################################################################
# instantiate the handles to the relevant collections.

process_name = 'RECO'
if 'UL' or 'HNL' in sample:  process_name = 'PAT'

## gen particles
handle_gen              = Handle('std::vector<reco::GenParticle>')
## offline objects
handle_reco_taus        = Handle('std::vector<pat::Tau>')
handle_vtx              = Handle('std::vector<reco::Vertex>')
# L1 taus
handle_l1_tau = Handle('BXVector<l1t::Tau>')
handle_l1_ele = Handle('BXVector<l1t::EGamma>')
## HLT objects 
handle_hlt_taus = Handle('std::vector<reco::PFTau>')
handle_hlt_ele  = Handle('std::vector<reco::RecoEcalCandidate>')
handle_hlt_gsfele  = Handle('std::vector<reco::Electron>')

handle_beamspot  = Handle('reco::BeamSpot') 
## HLT ancillary variables
handle_IP_displ     = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_iso_displ    = Handle('reco::PFTauDiscriminator')

## HLT filters
handle_hlt_tau_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_hlt_overlap_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l1ele_filter    = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter1   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter2   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter3   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter4   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter5   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter6   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter7   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter8   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_isoele_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter9  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_lastsingleeg_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 


handle_l1gsfele_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter1  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter2  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter3  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter4  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter5  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter6  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter7  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter8  = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_tightele_filter9  = Handle('trigger::TriggerFilterObjectWithRefs') 


handle_r9           = Handle('reco::RecoEcalCandidateIsolationMap')
handle_clustershape = Handle('reco::RecoEcalCandidateIsolationMap')
handle_ecaliso      = Handle('reco::RecoEcalCandidateIsolationMap')
handle_hcaliso      = Handle('reco::RecoEcalCandidateIsolationMap')
handle_trkiso       = Handle('reco::RecoEcalCandidateIsolationMap')

handle_ebhits       = Handle('EcalRecHitCollection')
handle_eehits       = Handle('EcalRecHitCollection')
triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
triggerBitsMyHLT, triggerBitLabelMyHLT = Handle("edm::TriggerResults"), ("TriggerResults","","MYHLT")

handles = {}

handles['hlt_r9']      = [('hltEgammaR9ID', '', 'MYHLT'), handle_r9, False]
handles['hlt_clsh']    = [('hltEgammaClusterShape', '', 'MYHLT'), handle_clustershape, False]
handles['hlt_ecaliso'] = [('hltEgammaEcalPFClusterIso', '', 'MYHLT'), handle_ecaliso, False]
handles['hlt_hcaliso'] = [('hltEgammaHcalPFClusterIso', '', 'MYHLT'), handle_hcaliso, False]
handles['hlt_trkiso']  = [('hltEgammaHollowTrackIso', '', 'MYHLT'), handle_trkiso, False]
handles['hlt_ebhits']  = [('hltEcalRecHit','EcalRecHitsEB', 'MYHLT'), handle_ebhits, False]
handles['hlt_eehits']  = [('hltEcalRecHit','EcalRecHitsEE', 'MYHLT'), handle_eehits, False]


handles['gen_particles']     = [('prunedGenParticles', '', process_name) , handle_gen, False]
if 'gmsb' in sample:
    handles['gen_particles'] = [('genParticlePlusGeant', '', 'SIM') , handle_gen, False]


## L1 objects
handles['l1_eles'] = [('hltGtStage2Digis','EGamma','MYHLT'), handle_l1_ele, False]
handles['l1_taus'] = [('hltGtStage2Digis','Tau','MYHLT'), handle_l1_tau, False]
# handles['l1_taus'] = [('caloStage2Digis','Tau','RECO'), handle_l1_tau, False]
# handles['l1_taus'] = [('caloStage2Digis','Tau','RECO'), handle_l1_tau, False]
## HLT objects 
handles['hlt_taus'] = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_taus, False]
handles['hlt_eles'] = [('hltEgammaCandidates', '', 'MYHLT'), handle_hlt_ele, False]
handles['hlt_gsfeles'] = [('hltEgammaGsfElectrons', '', 'MYHLT'), handle_hlt_gsfele, False]

handles['beamspot'] = [('hltOnlineBeamSpot', '', 'MYHLT'), handle_beamspot, False]
## HLT ancillary variables
handles['hlt_IP_displ']   = [('hltHpsPFTauTransverseImpactParameters', '', 'MYHLT'), handle_IP_displ, False]
handles['hlt_iso_displ']  = [('hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator', '', 'MYHLT'), handle_iso_displ, False]

## HLT filters
handles['hlt_tau_filter']  = [('hltHpsDisplacedPhotonMediumChargedIsoDisplPFTau26TrackPt1L1HLTMatchedGlob', '', 'MYHLT'), handle_hlt_tau_filter, False]
handles['hlt_overlap_filter']  = [('hltHpsOverlapFilterDisplacedEle22DisplPFTau26', '', 'MYHLT'), handle_hlt_overlap_filter, False]

handles['hlt_l1ele_filter']   = [('hltL1sBigORLooseIsoEGXXerIsoTauYYerdRMin0p3'   , '', 'MYHLT'), handle_l1ele_filter, False]
handles['hlt_eg_filter1']     = [('hltEG22EtFilterForEGTau'                               , '', 'MYHLT'), handle_l3ele_filter1, False]
handles['hlt_eg_filter2']     = [('hltEG22HEFilterForEGTau'                               , '', 'MYHLT'), handle_l3ele_filter2, False]
handles['hlt_eg_filter3']     = [('hltEG22R9Id90CaloIdLIsoLR9IdForEGTauFilter'            , '', 'MYHLT'), handle_l3ele_filter3, False]
handles['hlt_eg_filter4']     = [('hltEG22R9Id90CaloIdLIsoLClusterShapeForEGTauFilter'    , '', 'MYHLT'), handle_l3ele_filter4, False]
handles['hlt_eg_filter5']     = [('hltEG22R9Id90CaloIdLIsoLEcalPFClusterIsoForEGTauFilter', '', 'MYHLT'), handle_l3ele_filter5, False]
handles['hlt_eg_filter6']     = [('hltEG22R9Id90CaloIdLIsoLHcalPFClusterIsoForEGTauFilter', '', 'MYHLT'), handle_l3ele_filter6, False]
handles['hlt_eg_filter7']     = [('hltEG22R9Id90CaloIdLIsoLHollowTrackIsoForEGTauFilter'  , '', 'MYHLT'), handle_l3ele_filter7, False]
handles['hlt_eg_filter8']     = [('hltEG22R9Id90CaloIdLIsoLHollowTrackIsoForEGTauFilter'  , '', 'MYHLT'), handle_l3ele_filter8, False]
handles['hlt_eg_filter9']     = [('hltEG22R9Id90CaloIdLIsoLDisplacedIdForEGTauFilter',      '', 'MYHLT'), handle_l3ele_filter9, False]

handles['hlt_l1singleeg_filt']  = [('hltL1sSingleEGNonIsoOrWithJetAndTauNoPS' , '', 'MYHLT'), handle_l3ele_filter1, False]
handles['hlt_lastsingleeg_filt']= [('hltEG30R9Id90CaloIdLIsoLDisplacedIdFilter', '', 'MYHLT'), handle_lastsingleeg_filter, False]

handles['l1_gsfele_filter']         = [('hltL1sSingleEGor'       , '', 'MYHLT'), handle_l1gsfele_filter, False]
handles['hlt_tightele_filter1']     = [('hltEGL1SingleEGOrFilter'                , '', 'MYHLT'), handle_tightele_filter1, False]
handles['hlt_tightele_filter2']     = [('hltEle30WPTightClusterShapeFilter'      , '', 'MYHLT'), handle_tightele_filter2, False]
handles['hlt_tightele_filter3']     = [('hltEle30WPTightHEFilter'                , '', 'MYHLT'), handle_tightele_filter3, False]
handles['hlt_tightele_filter4']     = [('hltEle30WPTightEcalIsoFilter'           , '', 'MYHLT'), handle_tightele_filter4, False]
handles['hlt_tightele_filter5']     = [('hltEle30WPTightHcalIsoFilter'           , '', 'MYHLT'), handle_tightele_filter5, False]
handles['hlt_tightele_filter6']     = [('hltEle30WPTightPixelMatchFilter'        , '', 'MYHLT'), handle_tightele_filter6, False]
handles['hlt_tightele_filter7']     = [('hltEle30WPTightPMS2Filter'              , '', 'MYHLT'), handle_tightele_filter7, False]
handles['hlt_tightele_filter8']     = [('hltEle30WPTightGsfOneOEMinusOneOPFilter', '', 'MYHLT'), handle_tightele_filter8, False]
handles['hlt_tightele_filter9']    = [('hltEle30WPTightGsfTrackIsoFilter', '', 'MYHLT'), handle_tightele_filter9, False]

## offline objects
handles['reco_taus']         = [('slimmedTaus', '', process_name),        handle_reco_taus, False]
handles['vtx']               = [('offlineSlimmedPrimaryVertices','','PAT'), handle_vtx, False]

met_paths = [
  "HLT_PFMET120_PFMHT120_IDTight", 
  "HLT_PFMET130_PFMHT130_IDTight", 
#   "HLT_PFMET140_PFMHT140_IDTight", 
#   "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", 
  "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", 
  "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", 
#   "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", 
#   "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg", 
#   'HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1',
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", 
  "HLT_MET105_IsoTrk50", 
  "HLT_MET120_IsoTrk50",
]
lep_paths = [
  "HLT_Ele30_WPTight_Gsf",
#   "HLT_Photon20"
  "HLT_Photon200",
#   "HLT_Photon30_R9Id90_CaloIdL_IsoL_DisplacedIdL"
]

##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print(('===> processing %d / %d event' %(i, totevents)))
    
    for k, v in list(handles.items()):
        setattr(ev, k, None)
        v[2] = False
        try:
            ev.getByLabel(v[0], v[1])
            setattr(ev, k, v[1].product())
            v[2] = True
        except:    
            v[2] = False
    
    gen_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and pp.status()==good_gen_status]
    # select only hadronically decaying taus
    gen_had_taus = [pp for pp in gen_taus if isGenHadTau(pp)]
    # select only taus decaying to electrons
    gen_ele_taus = [pp for pp in gen_taus if isGenLepTau(pp, 11)]
    # skip events where we do not have a tau_h + tau_ele
    if len(gen_had_taus) < 1 or len(gen_ele_taus) < 1 :  continue

    # select accepted tau mothers
    gen_moms = [imom for imom in ev.gen_particles if abs(imom.pdgId()) in mom_pdgId]

    ## calculate MET 
    gen_neutrinos = [pp for pp in ev.gen_particles if (abs(pp.pdgId())==12 or abs(pp.pdgId())==16)]
    good_gen_neutrinos = []
    for ineu in gen_neutrinos:
#       if ineu.numberOfMothers() > 1 : print ('number of neu moms: ' , ineu.numberOfMothers() )
      if abs(ineu.mother(0).pdgId()) == 15:
        good_gen_neutrinos.append(ineu)
                    
    good_gen_lsps = []
    gen_lsps = [pp for pp in ev.gen_particles if abs(pp.pdgId())==1000022 and pp.status()==1]
    for ilsp in gen_lsps:
#       if ilsp.numberOfMothers() > 1 : print ('number of lsp moms: ' , ilsp.numberOfMothers() )
      if abs(ilsp.mother(0).pdgId()) in mom_pdgId:
        good_gen_lsps.append(ilsp)
    
    gen_met = 0
    for pp in good_gen_lsps+good_gen_neutrinos:
      gen_met = gen_met + pp.pt()



#     print (len(gen_had_taus), len(gen_ele_taus))
    # info about tau_h
    for gg in gen_had_taus:
        gg.bestmom = None
        tau_moms = [imom for imom in gen_moms if isAncestor(imom, gg) ]
        if len(tau_moms) > 0:  gg.bestmom = tau_moms[0]

    # info about tau_ele
    for gg in gen_ele_taus:
        gg.bestmom = None
        tau_moms = [imom for imom in gen_moms if isAncestor(imom, gg)]
        if len(tau_moms) > 0:  gg.bestmom = tau_moms[0]

    gen_ele_taus = [gg for gg in gen_ele_taus if gg.bestmom !=None]
    gen_had_taus = [gg for gg in gen_had_taus if gg.bestmom !=None]
    
#     for gg in gen_ele_taus: print(('ele mom: ', gg.bestmom.pdgId()))
#     for gg in gen_had_taus: print(('had mom: ', gg.bestmom.pdgId()))

    if len(gen_ele_taus) == 0 or  len(gen_had_taus) == 0 : 
      continue
    
    ## temporary!
    if len(gen_ele_taus) > 1 : 
      print('more than one gen ele tau')
      continue
    if len(gen_had_taus) > 1 : 
      print('more than one gen had tau')
      continue
      
    
    if gen_ele_taus[0].bestmom.pdgId() != -gen_had_taus[0].bestmom.pdgId() :
      print('tau_ele and tau_had coming from different moms')

    ######################################################################
    ## save variables for tau->ele
    for gg in gen_ele_taus:
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)
            
        ### find the electron, child of the tau 
        gg.dau = None
        for idau in range(gg.numberOfDaughters()):
            if abs(gg.daughter(idau).pdgId()) == 11: 
                gg.dau = gg.daughter(idau)
                break
        if gg.dau == None:  print ('is none')

        if gg.bestmom == None or gg.dau == None:
            continue
        
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

        if 'taugun' in sample:
            gg.lxy = sqrt(pow(gg.dau.vx(),2)+pow(gg.dau.vy(),2))
        else:    
            gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
            gg.lxy_ele = sqrt(pow(gg.dau.vx()-gg.bestmom.vx(),2)+pow(gg.dau.vy()-gg.bestmom.vy(),2))

        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))

    ######################################################################
    ## save variables for tau->had
    for gg in gen_had_taus:

        ## reset other attributes
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])

        ### find first dau to be used for vertices 
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print ('is none')
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()
        gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))


#     ######################################################################

    if handles['hlt_l1ele_filter'][2]:  
        ele_filter_product = ev.hlt_l1ele_filter.l1tegammaRefs()
        gen_ele_taus = findMatchToGen(gen_ele_taus, ele_filter_product, 'l1_ele_filter', 0.3)

    for ifilter in range(1,10):
        if handles['hlt_eg_filter%s'%ifilter][2]:  
          ele_filter_product = getattr(ev, 'hlt_eg_filter%s'%ifilter).photonRefs()
          gen_ele_taus = findMatchToGen(gen_ele_taus, ele_filter_product, 'egamma%s'%ifilter)

    for ifilter in range(1,10):
        if handles['hlt_tightele_filter%s'%ifilter][2]:  
          tightele_filter_product = getattr(ev, 'hlt_tightele_filter%s'%ifilter).photonRefs()
          gen_ele_taus = findMatchToGen(gen_ele_taus, tightele_filter_product, 'tightele%s'%ifilter)

   ### add loop on l1 electrons and emulate isolation 
    if handles['l1_eles'][2]:  
        all_l1eles = []
        for i in range(ev.l1_eles.size(0)):
            all_l1eles.append(ev.l1_eles.at(0, i))
        for thel1ele,ifeat in itertools_product(all_l1eles, l1_feat_list):  
            setattr(thel1ele, ifeat, -9999.)

        ### this l1 passes my single photon
        if handles['hlt_l1singleeg_filt'][2]:    
            filter_product = ev.hlt_l1singleeg_filt.l1tegammaRefs()
            l1ele_pt_filter_list = [k.pt() for k in filter_product]
            for thel1ele in all_l1eles:
                thel1ele.passSingleL1Filter = 1 if thel1ele.pt() in l1ele_pt_filter_list else 0

        ### this l1 passes ele-tau l1 seed
        if handles['hlt_l1ele_filter'][2]:    
            filter_product = ev.hlt_l1ele_filter.l1tegammaRefs()
            l1eletau_pt_filter_list = [k.pt() for k in filter_product]
            for thel1ele in all_l1eles:
                thel1ele.passEleTauL1Filter = 1 if thel1ele.pt() in l1eletau_pt_filter_list else 0

        if handles['l1_gsfele_filter'][2]:    
            filter_product = ev.l1_gsfele_filter.l1tegammaRefs()
            l1gsfele_pt_filter_list = [k.pt() for k in filter_product]
            for thel1ele in all_l1eles:
                thel1ele.passGsfEleL1Filter = 1 if thel1ele.pt() in l1gsfele_pt_filter_list else 0

        gen_ele_taus = findMatchToGen(gen_ele_taus, all_l1eles, 'l1_ele', 0.3)

    # match hlt displaced photons to gen taus
    counter_ele = 0
    if handles['hlt_eles'][2]:  
        all_eles = [pp for pp in ev.hlt_eles ]
        for theele,ifeat in itertools_product(all_eles, ele_feat_list):  
            setattr(theele, ifeat, -9999.)

        if handles['hlt_r9'][2]:
            r9_product = ev.hlt_r9
            pt_r9_dict = {}
            for k in range(len(r9_product)):
                pt_r9_dict[r9_product.keys()[k].pt()] = r9_product.values()[k]

        if handles['hlt_clsh'][2]:
            clsh_product = ev.hlt_clsh
            pt_clsh_dict = {}
            for k in range(len(clsh_product)):
                pt_clsh_dict[clsh_product.keys()[k].pt()] = clsh_product.values()[k]

        if handles['hlt_ecaliso'][2]:
            eiso_product = ev.hlt_ecaliso
            pt_eiso_dict = {}
            for k in range(len(eiso_product)):
                pt_eiso_dict[eiso_product.keys()[k].pt()] = eiso_product.values()[k]

        if handles['hlt_hcaliso'][2]:
            hiso_product = ev.hlt_hcaliso
            pt_hiso_dict = {}
            for k in range(len(hiso_product)):
                pt_hiso_dict[hiso_product.keys()[k].pt()] = hiso_product.values()[k]

        if handles['hlt_trkiso'][2]:
            tiso_product = ev.hlt_trkiso
            pt_tiso_dict = {}
            for k in range(len(tiso_product)):
                pt_tiso_dict[tiso_product.keys()[k].pt()] = tiso_product.values()[k]

        if handles['hlt_ebhits'][2] and handles['hlt_eehits'][2]:
            ebhits_product = ev.hlt_ebhits
            eehits_product = ev.hlt_eehits
            pt_smin_smaj_dict = {}
            for iele,theele in enumerate(all_eles):
              ## from https://github.com/cms-sw/cmssw/blob/c291814dc4517d36c7cb452b8c5e2a8850920815/HLTrigger/Egamma/plugins/HLTDisplacedEgammaFilter.cc#L115
              SCseed =   theele.superCluster().seed()
              rechits = ebhits_product
              if abs(theele.eta()) > 1.479: 
                rechits = eehits_product
              
              pt_smin_smaj_dict[theele.pt()] = [ROOT.EcalClusterTools.cluster2ndMoments(SCseed.get(),rechits).sMin, ROOT.EcalClusterTools.cluster2ndMoments(SCseed.get(),rechits).sMaj]
    
        for theele in all_eles:
            try:    theele.r9    = pt_r9_dict[theele.pt()]
            except: theele.r9    =  -9999.
    
            try:    theele.clsh = pt_clsh_dict[theele.pt()]
            except: theele.clsh =  -9999.
    
            try:    theele.ecaliso = pt_eiso_dict[theele.pt()]
            except: theele.ecaliso =  -9999.
    
            try:    theele.hcaliso = pt_hiso_dict[theele.pt()]
            except: theele.hcaliso =  -9999.

            try:    theele.trkiso  = pt_tiso_dict[theele.pt()]
            except: theele.trkiso  =  -9999.

            try:    
              theele.smin  = pt_smin_smaj_dict[theele.pt()][0]
              theele.smaj  = pt_smin_smaj_dict[theele.pt()][1]
            except: 
              theele.smin  =  -9999.
              theele.smaj  =  -9999.

        if handles['hlt_eg_filter9'][2]:
            filter_product = ev.hlt_eg_filter9.photonRefs()
            ele_pt_filter_list = [k.pt() for k in filter_product]
            for theele in all_eles:
                if theele.pt() in ele_pt_filter_list:
                    counter_ele = counter_ele + 1
                    theele.passFilters = 1  
                else:
                    theele.passFilters = 0
                
        if handles['hlt_lastsingleeg_filt'][2]:
            filter_product = ev.hlt_lastsingleeg_filt.photonRefs()
            ele_pt_filter_list = [k.pt() for k in filter_product]
            for theele in all_eles:
                theele.passSingleEGFilters = 1 if theele.pt() in ele_pt_filter_list else 0
                
        if handles['hlt_l1singleeg_filt'][2]:
            filter_product = ev.hlt_l1singleeg_filt.photonRefs()
            ele_pt_filter_list = [k.pt() for k in filter_product]
            for theele in all_eles:
                theele.passSingleL1Filter = 1 if theele.pt() in ele_pt_filter_list else 0

        gen_ele_taus = findMatchToGen(gen_ele_taus, all_eles, 'hlt_ele')
        gen_ele_taus = findMatchToGenMother(gen_ele_taus, all_eles, 'hlt_mom_to_ele')


    # match hlt GSF electrons to gen taus
    if handles['hlt_eles'][2]:  
        all_gsfeles = [pp for pp in ev.hlt_eles ]
#     if handles['hlt_gsfeles'][2]:  
#         all_gsfeles = [pp for pp in ev.hlt_gsfeles ]
        for theele,ifeat in itertools_product(all_gsfeles, ele_feat_list):  
            setattr(theele, ifeat, -9999.)

        if handles['hlt_tightele_filter9'][2]:
            filter_product = ev.hlt_tightele_filter9.photonRefs()
#             pdb.set_trace()
            gsfele_pt_filter_list = [k.pt() for k in filter_product]
#             filter_product2 = ev.hlt_tightele_filter9.electronRefs()
#             gsfele_pt_filter_list2 = [k.pt() for k in filter_product2]
            
            for iele,theele in enumerate(all_gsfeles):
                if theele.pt() in gsfele_pt_filter_list: 
                  theele.passFilters = 1                  
                else: theele.passFilters = 0
#             pdb.set_trace()

        gen_ele_taus = findMatchToGen(gen_ele_taus, all_gsfeles, 'hlt_gsfele')
        gen_ele_taus = findMatchToGenMother(gen_ele_taus, all_gsfeles, 'hlt_mom_to_gsfele')
#         pdb.set_trace()


    ######################################
    if handles['l1_taus'][2]:  
        all_l1taus = []
        for i in range(ev.l1_taus.size(0)):
            all_l1taus.append(ev.l1_taus.at(0, i))
        for thel1tau,ifeat in itertools_product(all_l1taus, l1_feat_list):  
            setattr(thel1tau, ifeat, -9999.)
        gen_had_taus = findMatchToGen(gen_had_taus, all_l1taus, 'l1_tau', 0.3)
        
    counter_taus = 0      
    if handles['hlt_taus'][2]:  

        all_taus = [pp for pp in ev.hlt_taus ]
        for thetau,ifeat in itertools_product(all_taus, pftau_feat_list):  
            setattr(thetau, ifeat, -9999.)
        pt_filter_list = []
        
        if handles['hlt_tau_filter'][2]:
            filter_product = ev.hlt_tau_filter.pftauRefs()
            pt_filter_list = [k.pt() for k in filter_product]

        if handles['hlt_IP_displ'][2]:
            ip_product = ev.hlt_IP_displ
            pt_ip_dict = {}
            for k in range(len(ip_product)):
                pt_ip_dict[ip_product.key(k).pt()] = [ip_product.value(k).dxy(), ip_product.value(k).dxy_error(), 
                                                      ip_product.value(k).ip3d(), ip_product.value(k).ip3d_error()] 

        for itau,thetau in enumerate(all_taus):
            
            if thetau.pt() in pt_filter_list: 
                thetau.passFilters = 1
                counter_taus = counter_taus + 1
            else:
                thetau.passFilters = 0

            try:  
                thetau.dxy     = pt_ip_dict[thetau.pt()][0]
                thetau.dxyerr  = pt_ip_dict[thetau.pt()][1]
                thetau.ip3d    = pt_ip_dict[thetau.pt()][2]
                thetau.ip3derr = pt_ip_dict[thetau.pt()][3]
#                 print thetau.dxy 
            except:  pass

        gen_had_taus = findMatchToGen(gen_had_taus, all_taus, 'hlt_tau')

#     ######################################################################################

    pass_overlap = False
    if handles['hlt_overlap_filter'][2]:
        pass_overlap = len(ev.hlt_overlap_filter.pftauRefs()) > 0 and len(ev.hlt_overlap_filter.photonRefs()) > 0 

    pass_emu_hlt = 0
    if counter_ele > 0 and counter_taus > 0 :  pass_emu_hlt = 1
    pass_emu_hlt_overlap = 0
    if counter_ele > 0 and counter_taus > 0 and pass_overlap:  pass_emu_hlt_overlap = 1
    pass_eleonly_hlt = 0
    if counter_ele > 0 :  pass_eleonly_hlt = 1

#     ######################################################################################
    pass_met_hlt = 0
    pass_hlts = {}
    pass_lep_hlt = 0
    pass_lep_hlts = {}
    ev.getByLabel(triggerBitLabel, triggerBits)
    names = ev.object().triggerNames(triggerBits.product())
    for i in range(triggerBits.product().size()):   
        if names.triggerName(i).split('_v')[0] in met_paths:
#           print ("Trigger ", names.triggerName(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)"))
          if triggerBits.product().accept(i) :
#             print ("Trigger ", names.triggerName(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)"))
            pass_hlts[names.triggerName(i).split('_v')[0]] = 1
            if names.triggerName(i).split('_v')[0] != 'HLT_Photon30_R9Id90_CaloIdL_IsoL_DisplacedIdL':  
              pass_met_hlt = 1
          else:  
            pass_hlts[names.triggerName(i).split('_v')[0]] = 0

        if names.triggerName(i).split('_v')[0] in lep_paths:
          if triggerBits.product().accept(i) :
            pass_lep_hlts[names.triggerName(i).split('_v')[0]] = 1
            if names.triggerName(i).split('_v')[0] != 'HLT_Photon30_R9Id90_CaloIdL_IsoL_DisplacedIdL':  
              pass_lep_hlt = 1
          else:  
            pass_lep_hlts[names.triggerName(i).split('_v')[0]] = 0

#     ######################################################################################
    pass_EGTau_hlt = 0
    ev.getByLabel(triggerBitLabelMyHLT, triggerBitsMyHLT)
    names = ev.object().triggerNames(triggerBitsMyHLT.product())
    for i in range(triggerBitsMyHLT.product().size()):   
        if names.triggerName(i).split('_v')[0] == 'HLT_Photon30_R9Id90_CaloIdL_IsoL_DisplacedIdL_displacedTau':
          if triggerBitsMyHLT.product().accept(i) :
            pass_EGTau_hlt = 1
          else:  
            pass_EGTau_hlt = 0
          break  

    
    
#     print ('len taus: ', len(gen_had_taus))
#     print ('len eles: ', len(gen_ele_taus))

    # fill the ntuple: each gen tau makes an entry
    # since only one ele and one tau per event, remove the loop

    for k, v in tofill_gen.items(): tofill_gen[k] = -9999. # initialise before filling

    tofill_gen['run'               ] = ev.eventAuxiliary().run()
    tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
    tofill_gen['event'             ] = ev.eventAuxiliary().event()
    tofill_gen['gen_met'           ] = gen_met

    gg = gen_ele_taus[0]
    if gg.bestmom == None or gg.dau == None:  continue
    tofill_gen['gen_ele_ndau'      ] = gg.numberOfDaughters()
    tofill_gen['gen_ele_pt'        ] = gg.pt()
    tofill_gen['gen_ele_eta'       ] = gg.eta()
    tofill_gen['gen_ele_phi'       ] = gg.phi()
    tofill_gen['gen_ele_charge'    ] = gg.charge()
    tofill_gen['gen_ele_lxy'       ] = gg.lxy
    tofill_gen['gen_ele_lxy_ele'   ] = gg.lxy_ele
    tofill_gen['gen_ele_dxy'       ] = gg.dxy
    tofill_gen['gen_ele_vx'        ] = gg.bestmom.vx()  ## user defined ones (mother prod. vertex)
    tofill_gen['gen_ele_vy'        ] = gg.bestmom.vy()
    tofill_gen['gen_ele_vz'        ] = gg.bestmom.vz()
    tofill_gen['gen_ele_cosxy'     ] = gg.cosxy
    tofill_gen['gen_ele_mom_mass'  ] = gg.bestmom.mass()
    tofill_gen['gen_ele_mom_pt'    ] = gg.bestmom.pt()
    tofill_gen['gen_ele_mom_eta'   ] = gg.bestmom.eta()
    tofill_gen['gen_ele_mom_phi'   ] = gg.bestmom.phi()
    tofill_gen['gen_ele_momct2d'   ] = gg.momct2d
    tofill_gen['gen_ele_vis_mass'  ] = gg.vismass()
    tofill_gen['gen_ele_vis_pt'    ] = gg.vispt()
    tofill_gen['gen_ele_vis_eta'   ] = gg.viseta()
    tofill_gen['gen_ele_vis_phi'   ] = gg.visphi()

    if hasattr(gg, 'hlt_ele') and gg.hlt_ele:
        tofill_gen['hlt_ele_pt'] = gg.hlt_ele.pt()
        tofill_gen['hlt_ele_eta'] = gg.hlt_ele.eta()
        tofill_gen['hlt_ele_phi'] = gg.hlt_ele.phi()
        tofill_gen['hlt_ele_r9'] = gg.hlt_ele.r9 
        tofill_gen['hlt_ele_clsh'] = gg.hlt_ele.clsh
        tofill_gen['hlt_ele_eiso'] = gg.hlt_ele.ecaliso
        tofill_gen['hlt_ele_hiso'] = gg.hlt_ele.hcaliso
        tofill_gen['hlt_ele_tiso'] = gg.hlt_ele.trkiso
        tofill_gen['hlt_ele_smin'] = gg.hlt_ele.smin
        tofill_gen['hlt_ele_smaj'] = gg.hlt_ele.smaj
        tofill_gen['hlt_ele_passFilters'] = gg.hlt_ele.passFilters 
        tofill_gen['hlt_ele_passSingleEGFilters'] = gg.hlt_ele.passSingleEGFilters 
        tofill_gen['hlt_ele_passSingleL1Filter'] = gg.hlt_ele.passSingleL1Filter
        

    if hasattr(gg, 'hlt_gsfele') and gg.hlt_gsfele:
        tofill_gen['hlt_gsfele_pt']  = gg.hlt_gsfele.pt()
        tofill_gen['hlt_gsfele_eta'] = gg.hlt_gsfele.eta()
        tofill_gen['hlt_gsfele_phi'] = gg.hlt_gsfele.phi()
        tofill_gen['hlt_gsfele_passFilters'] = gg.hlt_gsfele.passFilters 
        
        
    if hasattr(gg, 'l1_ele') and gg.l1_ele:
        tofill_gen['l1_ele_pt'       ] = gg.l1_ele.pt()
        tofill_gen['l1_ele_eta'      ] = gg.l1_ele.eta()
        tofill_gen['l1_ele_phi'      ] = gg.l1_ele.phi()
        tofill_gen['l1_ele_charge'   ] = gg.l1_ele.charge()
        tofill_gen['l1_ele_hwiso'    ] = gg.l1_ele.hwIso()

        tofill_gen['l1_ele_passselel1'   ] = gg.l1_ele.passSingleL1Filter
        tofill_gen['l1_ele_passeletaul1' ] = gg.l1_ele.passEleTauL1Filter
        tofill_gen['l1_ele_passgsfelel1' ] = gg.l1_ele.passGsfEleL1Filter

    if hasattr(gg, 'l1_ele_filter') and gg.l1_ele_filter:
        tofill_gen['l1_ele_filter_pt'       ] = gg.l1_ele_filter.pt()
        tofill_gen['l1_ele_filter_eta'      ] = gg.l1_ele_filter.eta()
        tofill_gen['l1_ele_filter_phi'      ] = gg.l1_ele_filter.phi()
        tofill_gen['l1_ele_filter_charge'   ] = gg.l1_ele_filter.charge()
        
    if hasattr(gg, 'hlt_mom_to_gsfele') and gg.hlt_mom_to_gsfele:
        tofill_gen['hlt_gsfele_matched_mom' ] = 1
    if hasattr(gg, 'hlt_mom_to_ele') and gg.hlt_mom_to_ele:
        tofill_gen['hlt_ele_matched_mom' ] = 1
        
        

    for ifilt in range(1,10):
      if hasattr(gg, 'egamma%s'%ifilt) and getattr(gg, 'egamma%s'%ifilt):
        this_eg = getattr(gg, 'egamma%s'%ifilt)
        tofill_gen['filter%s_ele_pt'%ifilt     ] = this_eg.pt()
        tofill_gen['filter%s_ele_eta'%ifilt    ] = this_eg.eta()
        tofill_gen['filter%s_ele_phi'%ifilt    ] = this_eg.phi()
#         tofill_gen['filter%s_ele_charge'%ifilt ] = this_eg.charge()

    for ifilt in range(1,10):
      if hasattr(gg, 'tightele%s'%ifilt) and getattr(gg, 'tightele%s'%ifilt):
        this_eg = getattr(gg, 'tightele%s'%ifilt)
        tofill_gen['filter%s_tightele_pt'%ifilt     ] = this_eg.pt()
        tofill_gen['filter%s_tightele_eta'%ifilt    ] = this_eg.eta()
        tofill_gen['filter%s_tightele_phi'%ifilt    ] = this_eg.phi()


    tt = gen_had_taus[0]
    if tt.bestmom == None or tt.dau == None:  continue

    tofill_gen['gen_tau_ndau'      ] = tt.numberOfDaughters()
    tofill_gen['gen_tau_pt'        ] = tt.pt()
    tofill_gen['gen_tau_eta'       ] = tt.eta()
    tofill_gen['gen_tau_phi'       ] = tt.phi()
    tofill_gen['gen_tau_charge'    ] = tt.charge()
    tofill_gen['gen_tau_decaymode' ] = tt.decayMode
    tofill_gen['gen_tau_lxy'       ] = tt.lxy
    tofill_gen['gen_tau_dxy'       ] = tt.dxy
    tofill_gen['gen_tau_vx'        ] = tt.bestmom.vx()  ## user defined ones (mother prod. vertex)
    tofill_gen['gen_tau_vy'        ] = tt.bestmom.vy()
    tofill_gen['gen_tau_vz'        ] = tt.bestmom.vz()
    tofill_gen['gen_tau_cosxy'     ] = tt.cosxy
    tofill_gen['gen_tau_mom_mass'  ] = tt.bestmom.mass()
    tofill_gen['gen_tau_momct2d'   ] = tt.momct2d
    tofill_gen['gen_tau_vis_mass'  ] = tt.vismass()
    tofill_gen['gen_tau_vis_pt'    ] = tt.vispt()
    tofill_gen['gen_tau_vis_eta'   ] = tt.viseta()
    tofill_gen['gen_tau_vis_phi'   ] = tt.visphi()

    if hasattr(tt, 'l1_tau') and tt.l1_tau:
        tofill_gen['l1_tau_pt'       ] = tt.l1_tau.pt()
        tofill_gen['l1_tau_eta'      ] = tt.l1_tau.eta()
        tofill_gen['l1_tau_phi'      ] = tt.l1_tau.phi()
        tofill_gen['l1_tau_iso'      ] = tt.l1_tau.hwIso()

    if hasattr(tt, 'hlt_tau') and tt.hlt_tau:
        tofill_gen['hlt_tau_mass'       ] = tt.hlt_tau.mass()
        tofill_gen['hlt_tau_pt'         ] = tt.hlt_tau.pt()
        tofill_gen['hlt_tau_eta'        ] = tt.hlt_tau.eta()
        tofill_gen['hlt_tau_phi'        ] = tt.hlt_tau.phi()
        tofill_gen['hlt_tau_charge'     ] = tt.hlt_tau.charge()
        tofill_gen['hlt_tau_decaymode'  ] = tt.hlt_tau.decayMode()

        tofill_gen['hlt_tau_dxy'        ] = tt.hlt_tau.dxy 
        tofill_gen['hlt_tau_dxyerr'     ] = tt.hlt_tau.dxyerr 
        tofill_gen['hlt_tau_ip3d'       ] = tt.hlt_tau.ip3d 
        tofill_gen['hlt_tau_ip3derr'    ] = tt.hlt_tau.ip3derr
        tofill_gen['hlt_tau_passFilters'] = tt.hlt_tau.passFilters 
            

    tofill_gen['gen_tau_ele_dr'   ] = deltaR(gg.eta(),gg.phi(),tt.viseta(), tt.visphi())

    tofill_gen['pass_met_hlt'            ] = pass_met_hlt
    tofill_gen['pass_lep_hlt'            ] = pass_lep_hlt
    tofill_gen['pass_EGTau_hlt'          ] = pass_EGTau_hlt
    tofill_gen['pass_EGTau_obj'          ] = pass_emu_hlt_overlap
    tofill_gen['pass_EGOnly_obj'         ] = pass_eleonly_hlt
    
    for ipath in met_paths:
      tofill_gen['pass_%s'%ipath.strip('HLT_')] = pass_hlts[ipath]
    for ipath in lep_paths:
        tofill_gen['pass_%s'%ipath.strip('HLT_')] = pass_lep_hlts[ipath]
    ntuple_gen.Fill(array('f',tofill_gen.values()))


##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()