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
from treeVariablesHLT_mu import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer , genDecayModeGEANT, isAncestor, isGenLepTau# utility functions
from math import sqrt, pow
from copy import deepcopy as dc
from itertools import product as itertools_product


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
feat_list = ['lxy', 'dxy', 'visdxy', 'cosxy', 'momct', 'momct2d', 'mom_mass', 'pi_lxy', 'pi_cosxy' ]
pftau_feat_list = ['leadChargedCandPt', 'leadChargedCandPdgId', 'leadCandPt', 'leadCandPdgId', 'maxHCALPFClusterEt', 'nChargedHad', 'nGamma', 'sum_pt_charged', 'sum_pt_neutral', 'dxy', 'dxyerr', 'ip3d', 'ip3derr', 'isoVal', 'passChargedIso', 'passRelChargedIso', 'passAbsChargedIso', 'passFilters']
displmu_feat_list = ['dxy']
c_const = 299.792458
dR_cone = 0.1
dPt_cone = 0.2
good_gen_status = 2
mom_pdgId = [1000015]
if 'HNL' in sample:  mom_pdgId = [9900012]
if 'HNL' in sample and 'Dirac' in sample:  mom_pdgId = [9990012]
if 'gmsb' in sample:  
    mom_pdgId = [2000015, 1000015]
    good_gen_status = 8
if 'DY' in sample:  mom_pdgId = [23]

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_mu_tuple_{}_{}.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(list(zip(branches, [-9999.]*len(branches)))) # initialise all branches to unphysical -99       
print(outfile_gen.GetName()) 

##########################################################################################
# Get ahold of the events
f = open('%s_filelist.txt'%args.sample)
infile = f.readlines()[ifile]

os.environ['X509_USER_PROXY'] = '/afs/cern.ch/user/f/fiorendi/x509up_u58808' 
redirector = 'root://cms-xrd-global.cern.ch//'
if 'fnal' in sample:
    redirector = 'root://cmseos.fnal.gov//'
if 'eos' in sample:
    redirector = ''

events = Events(redirector+infile.strip())
# events = Events([
#                  '/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_3_0_pre6/src/HLTrigger/Configuration/test/outputHLT.root'
#                  '/eos/cms/store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC/crab_ntuples_gmsb_12_3_0pre6_fromMenu/220311_150825/0000_sara/outputHLT_24.root',
#                 ])

# print 'using a local file since fnal is down!!!!!!!!!!!!!!!!!!!!!'

print(infile)
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files


##########################################################################################
def findMatchToGen(gen_taus, hlt_taus, hlt_tau): 
 
    gen_taus_match = gen_taus 
    for gg in gen_taus_match: 
        setattr(gg, hlt_tau, None)
        bestcom = bestMatch(gg.visp4, hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
            setattr(gg,hlt_tau,bestcom[0])
    return gen_taus_match    


def findMatchToGenTAU(gen_taus, hlt_taus, hlt_tau): 
 
    gen_taus_match = gen_taus 
    for gg in gen_taus_match: 
        setattr(gg, hlt_tau, None)
        
        printGEN(gg)
        for ihlt in hlt_taus:
          printForDebug(ihlt, gg)
        
        bestcom = bestMatch(gg.visp4, hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
            setattr(gg,hlt_tau,bestcom[0])
            print ('matched')
    return gen_taus_match    

def findMatchToGenMom(gen_taus, hlt_taus, attr_name): 
 
    gen_taus_match = gen_taus 
#     gen_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
#     mom_collection = [itau.bestmom for itau in gen_taus]
    for gg in gen_taus_match: 
        setattr(gg, attr_name, None)
        bestcom = bestMatch(gg.bestmom.p4(), hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
            setattr(gg,attr_name,bestcom[0])
            print('found one match: ', gg.bestmom.p4())
            
    return gen_taus_match    



# def findMatchToTau(hlt_taus, hlt_tau): 
#  
#     print('tau match')
#     hlt_taus_match = hlt_taus 
#     for gg in hlt_taus_match: 
#         setattr(gg, hlt_tau, None)
#         setattr(gg, hlt_tau,0.2)
#     return hlt_taus_match    

def printForDebug(cc, gg):
    print('    pt: %.2f'%( cc.pt()), \
          '\t eta: %.2f'%( cc.eta()), \
          '\t phi: %.2f'%( cc.phi()), \
          '\t pdgID: %d'%( cc.pdgId()), \
          '\t dpT/pT = %.2f,  dR = %.2f '%(abs(cc.pt() - gg.vispt())/gg.vispt(), deltaR(cc.eta(),cc.phi(),gg.viseta(), gg.visphi())))

def printGEN(gg):
    print('GEN pt: %.2f'%(gg.vispt(),), \
          '\t eta: %.2f'%(gg.viseta()), \
          '\t phi: %.2f'%(gg.visphi()), \
          '\t dm: %.2f'%(gg.decayMode), '\n')

   
##########################################################################################
# instantiate the handles to the relevant collections.

process_name = 'RECO'
if 'UL' or 'HNL' in sample:  process_name = 'PAT'

# gen particles
handle_gen              = Handle('std::vector<reco::GenParticle>')
handle_reco_taus        = Handle('std::vector<pat::Tau>')
handle_L3mus            = Handle('std::vector<reco::RecoChargedCandidate>')
# L2 taus
handle_l2_taus          = Handle('std::vector<reco::CaloJet>')
handle_l2_isoTaus       = Handle('std::vector<reco::CaloJet>')
handle_hlt_pftaus       = Handle('std::vector<reco::PFTau>')
handle_hlt_pftaus_displ = Handle('std::vector<reco::PFTau>')

# PF
handle_hltPFs           = Handle('std::vector<reco::PFCandidate>')
## tracks 
handle_hltTracks        = Handle('std::vector<reco::Track>')
handle_hltIter4Tracks   = Handle('std::vector<reco::Track>')
handle_hltIter04Tracks  = Handle('std::vector<reco::Track>')
# vertices
handle_vtx              = Handle('std::vector<reco::Vertex>')
# L1 taus
handle_l1               = Handle('BXVector<l1t::Tau>')
handle_l1_mu            = Handle('BXVector<l1t::Muon>')
handle_my_l1_mu         = Handle('BXVector<l1t::Muon>')
#  offline tracks 
handle_lost_tracks      = Handle('std::vector<pat::PackedCandidate>')
handle_packed           = Handle('std::vector<pat::PackedCandidate>')

handle_IP_displ     = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_iso_displ    = Handle('reco::PFTauDiscriminator')
handle_absiso_displ = Handle('reco::PFTauDiscriminator')
handle_reliso_displ = Handle('reco::PFTauDiscriminator')
handle_isoval_displ = Handle('reco::PFTauDiscriminator')
ip_handle           = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_hlt_filter    = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l1mu_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3mu_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_isomu_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 

handle_l2disp_filter = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_casc_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_glbd_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 

handle_beamspot  = Handle('reco::BeamSpot') 


handles = {}
handles['gen_particles']     = [('prunedGenParticles', '', process_name) , handle_gen, False]
if 'gmsb' in sample:
    handles['gen_particles'] = [('genParticlePlusGeant', '', 'SIM') , handle_gen, False]

handles['beamspot']        = [('hltOnlineBeamSpot', '', 'MYHLT'), handle_beamspot, False]

handles['l1_muons']          = [('gmtStage2Digis','Muon','RECO'),         handle_l1_mu, False]
handles['my_l1_muons']       = [('hltGtStage2Digis','Muon','MYHLT'),      handle_my_l1_mu, False]

handles['l1_taus']           = [('caloStage2Digis','Tau','RECO'),            handle_l1, False]
handles['l2_taus']           = [('hltL2TauJetsL1TauSeeded', '', 'MYHLT'),    handle_l2_taus, False]
handles['l2_isoTaus']        = [('hltL2TauJetsIsoL1TauSeeded', '', 'MYHLT'), handle_l2_isoTaus, False]

## for release 12_3_0_pre3
if 'DY' not in sample:
    handles['l2_isoTaus']        = [('hltL2TauJetsIsoL1TauSeededGlob', '', 'MYHLT'), handle_l2_isoTaus, False]


handles['hlt_pftaus']        = [('hltHpsPFTauProducer', '', 'MYHLT'), handle_hlt_pftaus, False]
handles['hlt_pftaus_displ']  = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]

handles['hlt_IP_displ']   = [('hltHpsPFTauTransverseImpactParameters', '', 'MYHLT'), handle_IP_displ, False]
handles['hlt_iso_displ']  = [('hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator', '', 'MYHLT'), handle_iso_displ, False]
handles['hlt_absiso_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationDiscriminator', '', 'MYHLT'), handle_absiso_displ, False]
handles['hlt_reliso_displ']  = [('hltHpsDisplPFTauMediumRelativeChargedIsolationDiscriminator', '', 'MYHLT'), handle_reliso_displ, False]
handles['hlt_isoval_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationValue', '', 'MYHLT'), handle_isoval_displ, False]

handles['hlt_filter']  = [('hltHpsDisplacedMuMediumChargedIsoDisplPFTau20TrackPt1L1HLTMatchedGlob', '', 'MYHLT'), handle_hlt_filter, False]
# process.HLT_DisplacedMu24_glbDispl_displacedTau_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sBigORMu18erTauXXer2p1 + 
# hltPreDisplacedMu24glbDispldisplacedTau + hltL1fL1sBigORMu18erTauXXer2p1L1Filtered0 + 
#   HLTL2muonrecoSequenceNoVtx + hltL2fL1SingleMuf0L2NoVtxFiltered7DisplTau + 
#   HLTIterGlbDisplacedOneSequence + hltL3fSingleMuL1f0L2NVf7L3GlbDispl10 + 
#   HLTL2TauJetsL1TauSeededSequence + hltDisplMuL2Tau20eta2p2 + 
#   HLTL2p5IsoTauL1TauSeededGlobalSequence + hltDisplMuL2GlobIsoTau20eta2p2 + 
#   HLTGlobalPFTauDisplHPSSequence + HLTHPSSingleDisplPFTauPt20Eta2p1Trk1Glob + 
#   HLTHPSMediumChargedIsoDisplPFTauSequence + hltHpsSelectedPFTausTrackPt1MediumChargedIsolationGlobDispl + 
#   hltHpsL1JetsHLTDisplacedMuDisplPFTauTrackPt1MatchGlob + hltHpsDisplacedMuMediumChargedIsoDisplPFTau20TrackPt1L1HLTMatchedGlob + 
#   HLTDisplPFTauDxyProducer + hltHpsSingleMediumChargedIsoDisplPFTau20Dxy0p005 + HLTEndSequence )

handles['hlt_l1mu_filter']  = [('hltL1fL1sMu22L1Filtered0', '', 'MYHLT'), handle_l1mu_filter, False]
handles['hlt_l1mu_filter_18']  = [('hltL1sBigORMu18erTauXXer2p1', '', 'MYHLT'), handle_l1mu_filter, False]
handles['hlt_l3mu_filter']  = [('hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q', '', 'MYHLT'), handle_l3mu_filter, False]
handles['hlt_isomu_filter'] = [('hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07', '', 'MYHLT'), handle_isomu_filter, False]

handles['hlt_l2disp_filter']  = [('hltL2fL1SingleMuf0L2NoVtxFiltered7', '', 'MYHLT'), handle_l2disp_filter, False]
handles['hlt_cascade_filter'] = [('hltL3fSingleMuL1f0L2NVf7L3NoFiltersNoVtxFiltered10', '', 'MYHLT'), handle_casc_filter, False]
handles['hlt_glbdisp_filter'] = [('hltL3fSingleMuL1f0L2NVf7L3GlbDispl10', '', 'MYHLT'), handle_glbd_filter, False]

'''
trigger::TriggerFilterObjectWithRefs   "hltL1fForIterL3L1fL1sMu22L1Filtered0"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL1sSingleMu22"          ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL1fL1sMu22L1Filtered0"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL2fL1sSingleMu22L1f0L2Filtered10Q"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfecalIsoRhoFilteredEB0p14EE0p10"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3pfhcalIsoRhoFilteredHB0p16HE0p20"   ""                "MYHLT"
trigger::TriggerFilterObjectWithRefs   "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"   ""                "MYHLT"
'''

handles['hltTks']            = [('hltMergedTracks',  '', 'MYHLT'),                                        handle_hltTracks, False]
handles['hltIter04Tks']      = [('hltIter4MergedWithIter0ForTau', '', 'MYHLT'),                           handle_hltIter04Tracks, False]
handles['hltIter4Tks']       = [('hltDisplacedhltIter4PFlowTrackSelectionHighPurityForTau', '', 'MYHLT'), handle_hltIter4Tracks, False]
handles['hltMuMergedTks']    = [('hltPFMuonMerging', '', 'MYHLT'),     handle_hltTracks, False]

handles['reco_taus']         = [('slimmedTaus', '', process_name),        handle_reco_taus, False]
handles['lost_tracks']       = [('lostTracks', '',  process_name),        handle_lost_tracks, False]
handles['packed']            = [('packedPFCandidates', '', process_name), handle_packed, False]
handles['vtx']               = [('offlineSlimmedPrimaryVertices','','PAT'), handle_vtx, False]


if 'raw' in sample:
    handles['gen_particles'] = [('genParticles', '', 'SIM') , handle_gen, False]
    handles['l1_taus']       = [('hltGtStage2Digis','Tau','MYHLT'), handle_l1, False]




##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print('===> processing %d / %d event' %(i, totevents))
    
    for k, v in handles.items():
        setattr(ev, k, None)
        v[2] = False
        try:
            ev.getByLabel(v[0], v[1])
            setattr(ev, k, v[1].product())
            v[2] = True
        except:    
            v[2] = False
    
    # select only hadronically decaying taus
    gen_had_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                    pp.status()==good_gen_status and isGenHadTau(pp)]

    gen_mu_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                    pp.status()==good_gen_status and isGenLepTau(pp,13)]

    if len(gen_had_taus) < 1 or  len(gen_mu_taus) < 1 :  continue

    for gg in gen_had_taus:
        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)

    for gg in gen_mu_taus:
        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)


#     for gg in gen_mu_taus:   print('had mom: ', gg.bestmom.pdgId())
#     for gg in gen_had_taus:  print('mu mom: ', gg.bestmom.pdgId())

    if len(gen_mu_taus) == 0 or  len(gen_had_taus) == 0 : 
      continue
    
    ## temporary!
    if len(gen_mu_taus) > 1 : 
      print('more than one gen mu tau')
      continue
    if len(gen_had_taus) > 1 : 
      print('more than one gen had tau')
      continue
      
    
    if gen_mu_taus[0].bestmom.pdgId() != -gen_had_taus[0].bestmom.pdgId() :
      print('tau_mu and tau_had coming from different moms')

    ######################################################################
    ## save variables for tau->mu
    for gg in gen_mu_taus:
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)

        ### find first dau to be used for vertices 
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print ('is none')

        if gg.bestmom == None or gg.dau == None:
            continue
        
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

        if 'taugun' in sample:
            gg.lxy = sqrt(pow(gg.dau.vx(),2)+pow(gg.dau.vy(),2))
#             gg.visdxy = (gg.dau.vy()*gg.vispx() - gg.dau.vx()*gg.vispy())/gg.vispt()
        else:    
            gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
#             gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

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
        if 'gmsb' in sample:
            dm_string = genDecayModeGEANT( [d for d in finalDaughters(gg) if abs(d.pdgId()) not in [12, 14, 16]])
            gg.decayMode = tauDecayModes.nameToInt(dm_string)
        ### find first dau to be used for vertices 
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print ('is none')
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()
        gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))


    ######################################################################


    if handles['hlt_l1mu_filter'][2]:  
        mu_filter_product = ev.hlt_l1mu_filter.l1tmuonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'l1_mu')

    if handles['hlt_l1mu_filter_18'][2]:  
        mu18_filter_product = ev.hlt_l1mu_filter_18.l1tmuonRefs()
        gen_mu_taus_18 = findMatchToGen(gen_mu_taus, mu18_filter_product, 'l1_mu18')

    if handles['hlt_l3mu_filter'][2]:  
        mu_filter_product = ev.hlt_l3mu_filter.muonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'l3_mu')

    if handles['hlt_isomu_filter'][2]:  
        mu_filter_product = ev.hlt_isomu_filter.muonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'iso_mu')


    if handles['hlt_l2disp_filter'][2]:  
        mu_filter_product = ev.hlt_l2disp_filter.muonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'l2disp_mu')

    if handles['hlt_cascade_filter'][2]:  
        mu_filter_product = ev.hlt_cascade_filter.muonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'casc_mu')
        for gg in gen_mu_taus:
          if hasattr(gg, 'casc_mu') and gg.casc_mu:
            gg.casc_mu.dxy = gg.casc_mu.track().dxy(ev.beamspot)

    if handles['hlt_glbdisp_filter'][2]:          
        mu_filter_product = ev.hlt_glbdisp_filter.muonRefs()
        gen_mu_taus = findMatchToGen(gen_mu_taus, mu_filter_product, 'glb_mu')

        for gg in gen_mu_taus:
          if hasattr(gg, 'glb_mu') and gg.glb_mu:
            gg.glb_mu.dxy = gg.glb_mu.track().dxy(ev.beamspot)

            
            ### match taus only if muon is there
            if handles['hlt_pftaus_displ'][2]:  
              all_taus = [pp for pp in ev.hlt_pftaus_displ ]
              for thetau,ifeat in itertools_product(all_taus, pftau_feat_list):  
                setattr(thetau, ifeat, -9999.)

              ## first add info on IP and iso
              if handles['hlt_iso_displ'][2]:
                iso_product = ev.hlt_iso_displ
                pt_iso_dict = {}
                for k in range(len(iso_product)):
                  pt_iso_dict[iso_product.key(k).pt()] = iso_product.value(k)
      
              pt_filter_list = []
              if handles['hlt_filter'][2]:
                filter_product = ev.hlt_filter.pftauRefs()
                pt_filter_list = []
                for k in range(len(filter_product)):
                  pt_filter_list.append(filter_product[k].pt()) 
#                   print ('pt list: ', pt_filter_list)
              
              if handles['hlt_IP_displ'][2]:
                  ip_product = ev.hlt_IP_displ
                  pt_ip_dict = {}
                  for k in range(len(ip_product)):
                      pt_ip_dict[ip_product.key(k).pt()] = [ip_product.value(k).dxy(), ip_product.value(k).dxy_error(), 
                                                            ip_product.value(k).ip3d(), ip_product.value(k).ip3d_error()] 

              for itau,thetau in enumerate(all_taus):
                try:    
                  thetau.passChargedIso    = pt_iso_dict[thetau.pt()]
                except: thetau.passChargedIso    =  -9999.
          
                try:    thetau.passRelChargedIso = pt_reliso_dict[thetau.pt()]
                except: thetau.passRelChargedIso =  -9999.
          
                try:    thetau.passAbsChargedIso = pt_absiso_dict[thetau.pt()]
                except: thetau.passAbsChargedIso =  -9999.
          
                try:    thetau.isoVal            = pt_isoval_dict[thetau.pt()]
                except: thetau.isoVal            =  -9999.
          
                try:
                  if thetau.pt() in pt_filter_list:  thetau.passFilters = 1
                except: thetau.passFilters = 0

                try:  
                    thetau.dxy     = pt_ip_dict[thetau.pt()][0]
                    thetau.dxyerr  = pt_ip_dict[thetau.pt()][1]
                    thetau.ip3d    = pt_ip_dict[thetau.pt()][2]
                    thetau.ip3derr = pt_ip_dict[thetau.pt()][3]
                except:  pass
      
              gen_had_taus = findMatchToGen(gen_had_taus, all_taus, 'hlt_pftau_displ')
              
            ## access the l1 taus
            l1_taus = [tau for tau in ev.l1_taus if ev.l1_taus.getFirstBX()==0]
            l1_taus = [tau for tau in l1_taus if tau.pt()>15.]
            gen_had_taus = findMatchToGen(gen_had_taus, l1_taus, 'l1_tau')
  

        

#     # match L2 taus to gen taus
#     if handles['l2_taus'][2]:  
#         gen_taus = findMatchToGen(gen_taus, ev.l2_taus, 'l2_tau')
#  
#     # match L2 iso taus to gen taus
#     if handles['l2_isoTaus'][2]:  
#         gen_taus = findMatchToGen(gen_taus, ev.l2_isoTaus, 'l2_isotau')


#     # match hlt displced PF taus to gen taus
#     if handles['hlt_pftaus_displ'][2]:  
#         all_taus = [pp for pp in ev.hlt_pftaus_displ ]
#         for thetau,ifeat in itertools_product(all_taus, pftau_feat_list):  
#             setattr(thetau, ifeat, -9999.)
# 

#         if handles['hlt_reliso_displ'][2]:
#             reliso_product = ev.hlt_reliso_displ
#             pt_reliso_dict = {}
#             for k in range(len(reliso_product)):
#                 pt_reliso_dict[reliso_product.key(k).pt()] = reliso_product.value(k)
# 
#         if handles['hlt_absiso_displ'][2]:
#             absiso_product = ev.hlt_absiso_displ
#             pt_absiso_dict = {}
#             for k in range(len(absiso_product)):
#                 pt_absiso_dict[absiso_product.key(k).pt()] = absiso_product.value(k)
# 
#         if handles['hlt_isoval_displ'][2]:
#             isoval_product = ev.hlt_isoval_displ
#             pt_isoval_dict = {}
#             for k in range(len(isoval_product)):
#                 pt_isoval_dict[isoval_product.key(k).pt()] = isoval_product.value(k)
# 
#         if handles['hlt_filter'][2]:
#             filter_product = ev.hlt_filter.pftauRefs()
#             pt_filter_list = []
#             for k in range(len(filter_product)):
#                 pt_filter_list.append(filter_product[k].pt()) 
# 
#             for itau,thetau in enumerate(all_taus):
#                 try:    thetau.passChargedIso    = pt_iso_dict[thetau.pt()]
#                 except: thetau.passChargedIso    =  -9999.
#     
#                 try:    thetau.passRelChargedIso = pt_reliso_dict[thetau.pt()]
#                 except: thetau.passRelChargedIso =  -9999.
#     
#                 try:    thetau.passAbsChargedIso = pt_absiso_dict[thetau.pt()]
#                 except: thetau.passAbsChargedIso =  -9999.
#     
#                 try:    thetau.isoVal            = pt_isoval_dict[thetau.pt()]
#                 except: thetau.isoVal            =  -9999.
#     
#                 if thetau.pt() in pt_filter_list: thetau.passFilters = 1
#                 else:                             thetau.passFilters = 0
# 
#         if handles['hlt_IP_displ'][2]:
#             ip_product = ev.hlt_IP_displ
#             pt_ip_dict = {}
#             for k in range(len(ip_product)):
#                 pt_ip_dict[ip_product.key(k).pt()] = [ip_product.value(k).dxy(), ip_product.value(k).dxy_error(), 
#                                                       ip_product.value(k).ip3d(), ip_product.value(k).ip3d_error()] 
# 
#             for itau,thetau in enumerate(all_taus):
#                 try:  
#                     thetau.dxy     = pt_ip_dict[thetau.pt()][0]
#                     thetau.dxyerr  = pt_ip_dict[thetau.pt()][1]
#                     thetau.ip3d    = pt_ip_dict[thetau.pt()][2]
#                     thetau.ip3derr = pt_ip_dict[thetau.pt()][3]
# #                     print thetau.dxy 
#                 except:  pass
#         
# 
#         gen_taus = findMatchToGen(gen_taus, all_taus, 'hlt_pftau_displ')
#         for gg in gen_taus:
#             if hasattr(gg, 'hlt_pftau_displ') and gg.hlt_pftau_displ:
# 
#                 theLeadChargedCand = gg.hlt_pftau_displ.leadChargedHadrCand()
#                 theLeadPFCand      = gg.hlt_pftau_displ.leadCand()
#                 gg.hlt_pftau_displ.leadChargedCandPt     = theLeadChargedCand.pt()
#                 gg.hlt_pftau_displ.leadChargedCandPdgId  = theLeadChargedCand.pdgId()
#                 gg.hlt_pftau_displ.leadCandPt            = theLeadPFCand.pt()
#                 gg.hlt_pftau_displ.leadCandPdgId         = theLeadPFCand.pdgId()
# 
#                 gg.hlt_pftau_displ.maxHCALPFClusterEt    = gg.hlt_pftau_displ.maximumHCALPFClusterEt()
#                 gg.hlt_pftau_displ.nChargedHad           = gg.hlt_pftau_displ.signalChargedHadrCands().size() 
#                 gg.hlt_pftau_displ.nGamma                = gg.hlt_pftau_displ.signalGammaCands().size() 
# 
#                 sum_pt_charged = 0
#                 for ich in range(gg.hlt_pftau_displ.nChargedHad  ) :
#                     sum_pt_charged += gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()
# #                     print '\t ch: ', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pdgId(), '\t', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()
#     
#                 sum_pt_neutral = 0
#                 for ineu in range(gg.hlt_pftau_displ.nGamma  ) :
#                     sum_pt_neutral += gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()
# #                     print '\t neu: ', gg.hlt_pftau_displ.signalGammaCands()[ineu].pdgId(), '\t', gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()
# 
#                 gg.hlt_pftau_displ.sum_pt_charged =  sum_pt_charged  
#                 gg.hlt_pftau_displ.sum_pt_neutral =  sum_pt_neutral  
# 
#         gen_taus = findMatchToGenMom(gen_taus, all_taus, 'mom_matched_to_hlt')
#         
# 
#     ######################################################################################
# 
#     # match reco taus to gen taus
# #     for l1 in l1_taus : l1.gen_tau  = None # first initialise the matching to None
# #     for gg in gen_taus: gg.l1_tau   = None # first initialise the matching to None
# #     
# #     gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched
# #     
# #     for l1 in l1_taus:
# #         matches = [gg for gg in gen_taus_copy if deltaR(l1.p4(), gg.visp4)<dR_cone]
# #         ave_matches = len(matches)
# # #         if ave_matches > 1:  print ave_matches
# #         if not len(matches):
# #             continue
# #         matches.sort(key = lambda gg : deltaR(l1.p4(), gg.visp4))
# #         bestmatch = matches[0]
# #         l1.gen_tau = bestmatch
# #         bestmatch.l1_tau = l1
# #         gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]
# 
# #     ######################################################################################
#     # Save pileup information in MC #
    this_pu = -99.
# #     bx_vector = []
# # #     for ipuinfo in ev.pu:
# # #         if ipuinfo.getBunchCrossing() == 0:
# # #             this_pu = ipuinfo.getTrueNumInteractions()
# # #             break            
# # 
# # 
# #     ######################################################################################
#     
#     # fill the ntuple: each gen tau makes an entry
    for gg in gen_mu_taus:
        if gg.bestmom == None or gg.dau == None:  continue

        for k, v in tofill_gen.items(): tofill_gen[k] = -9999. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()
        tofill_gen['trueNI'            ] = this_pu

        tofill_gen['mu_gen_ndau'      ] = gg.numberOfDaughters()
        tofill_gen['mu_gen_pt'        ] = gg.pt()
        tofill_gen['mu_gen_eta'       ] = gg.eta()
        tofill_gen['mu_gen_phi'       ] = gg.phi()
        tofill_gen['mu_gen_charge'    ] = gg.charge()
#         tofill_gen['tau_gen_decaymode' ] = gg.decayMode
        tofill_gen['mu_gen_lxy'       ] = gg.lxy
        tofill_gen['mu_gen_dxy'       ] = gg.dxy
        tofill_gen['mu_gen_vx'        ] = gg.bestmom.vx()  ## user defined ones (mother prod. vertex)
        tofill_gen['mu_gen_vy'        ] = gg.bestmom.vy()
        tofill_gen['mu_gen_vz'        ] = gg.bestmom.vz()
        tofill_gen['mu_gen_cosxy'     ] = gg.cosxy
#         tofill_gen['tau_gen_momct'     ] = gg.momct
        tofill_gen['mu_gen_mom_mass'  ] = gg.bestmom.mass()
        tofill_gen['mu_gen_mom_pt'    ] = gg.bestmom.pt()
        tofill_gen['mu_gen_mom_eta'   ] = gg.bestmom.eta()
        tofill_gen['mu_gen_mom_phi'   ] = gg.bestmom.phi()
        tofill_gen['mu_gen_momct2d'   ] = gg.momct2d
        tofill_gen['mu_gen_vis_mass'  ] = gg.vismass()
        tofill_gen['mu_gen_vis_pt'    ] = gg.vispt()
        tofill_gen['mu_gen_vis_eta'   ] = gg.viseta()
        tofill_gen['mu_gen_vis_phi'   ] = gg.visphi()

        if hasattr(gg, 'l1_mu') and gg.l1_mu:
            tofill_gen['tau_l1mu_pt'       ] = gg.l1_mu.pt()
            tofill_gen['tau_l1mu_eta'      ] = gg.l1_mu.eta()
            tofill_gen['tau_l1mu_phi'      ] = gg.l1_mu.phi()
            tofill_gen['tau_l1mu_charge'   ] = gg.l1_mu.charge()

        if hasattr(gg, 'l1_mu18') and gg.l1_mu18:
            tofill_gen['tau_l1mu18_pt'       ] = gg.l1_mu18.pt()
            tofill_gen['tau_l1mu18_eta'      ] = gg.l1_mu18.eta()
            tofill_gen['tau_l1mu18_phi'      ] = gg.l1_mu18.phi()
            tofill_gen['tau_l1mu18_charge'   ] = gg.l1_mu18.charge()

        if hasattr(gg, 'l3_mu') and gg.l3_mu:
            tofill_gen['tau_l3mu_pt'       ] = gg.l3_mu.pt()
            tofill_gen['tau_l3mu_eta'      ] = gg.l3_mu.eta()
            tofill_gen['tau_l3mu_phi'      ] = gg.l3_mu.phi()
            tofill_gen['tau_l3mu_charge'   ] = gg.l3_mu.charge()

        if hasattr(gg, 'iso_mu') and gg.iso_mu:
            tofill_gen['tau_isomu_pt'       ] = gg.iso_mu.pt()
            tofill_gen['tau_isomu_eta'      ] = gg.iso_mu.eta()
            tofill_gen['tau_isomu_phi'      ] = gg.iso_mu.phi()
            tofill_gen['tau_isomu_charge'   ] = gg.iso_mu.charge()

        if hasattr(gg, 'l2disp_mu') and gg.l2disp_mu:
            tofill_gen['tau_l2disp_pt'       ] = gg.l2disp_mu.pt()
            tofill_gen['tau_l2disp_eta'      ] = gg.l2disp_mu.eta()
            tofill_gen['tau_l2disp_phi'      ] = gg.l2disp_mu.phi()
            tofill_gen['tau_l2disp_charge'   ] = gg.l2disp_mu.charge()
        if hasattr(gg, 'casc_mu') and gg.casc_mu:
            tofill_gen['tau_cascade_pt'       ] = gg.casc_mu.pt()
            tofill_gen['tau_cascade_eta'      ] = gg.casc_mu.eta()
            tofill_gen['tau_cascade_phi'      ] = gg.casc_mu.phi()
            tofill_gen['tau_cascade_charge'   ] = gg.casc_mu.charge()
            tofill_gen['tau_cascade_dxy'      ] = gg.casc_mu.dxy
        if hasattr(gg, 'glb_mu') and gg.glb_mu:
            tofill_gen['tau_glb_pt'       ] = gg.glb_mu.pt()
            tofill_gen['tau_glb_eta'      ] = gg.glb_mu.eta()
            tofill_gen['tau_glb_phi'      ] = gg.glb_mu.phi()
            tofill_gen['tau_glb_charge'   ] = gg.glb_mu.charge()
            tofill_gen['tau_glb_dxy'      ] = gg.glb_mu.dxy
        
        tau_gg =  gen_had_taus[0]   
        tofill_gen['tau_gen_ndau'      ] = tau_gg.numberOfDaughters()
        tofill_gen['tau_gen_pt'        ] = tau_gg.pt()
        tofill_gen['tau_gen_eta'       ] = tau_gg.eta()
        tofill_gen['tau_gen_phi'       ] = tau_gg.phi()
        tofill_gen['tau_gen_charge'    ] = tau_gg.charge()
        tofill_gen['tau_gen_decaymode' ] = tau_gg.decayMode
        tofill_gen['tau_gen_lxy'       ] = tau_gg.lxy
        tofill_gen['tau_gen_dxy'       ] = tau_gg.dxy
        tofill_gen['tau_gen_vx'        ] = tau_gg.bestmom.vx()  ## user defined ones (mother prod. vertex)
        tofill_gen['tau_gen_vy'        ] = tau_gg.bestmom.vy()
        tofill_gen['tau_gen_vz'        ] = tau_gg.bestmom.vz()
        tofill_gen['tau_gen_cosxy'     ] = tau_gg.cosxy

        tofill_gen['tau_gen_mom_mass'  ] = tau_gg.bestmom.mass()
        tofill_gen['tau_gen_momct2d'   ] = tau_gg.momct2d
        tofill_gen['tau_gen_vis_mass'  ] = tau_gg.vismass()
        tofill_gen['tau_gen_vis_pt'    ] = tau_gg.vispt()
        tofill_gen['tau_gen_vis_eta'   ] = tau_gg.viseta()
        tofill_gen['tau_gen_vis_phi'   ] = tau_gg.visphi()


        if hasattr(tau_gg, 'hlt_pftau_displ') and tau_gg.hlt_pftau_displ:
            tofill_gen['tau_hltPFdispltau_mass'     ] = tau_gg.hlt_pftau_displ.mass()
            tofill_gen['tau_hltPFdispltau_pt'       ] = tau_gg.hlt_pftau_displ.pt()
            tofill_gen['tau_hltPFdispltau_eta'      ] = tau_gg.hlt_pftau_displ.eta()
            tofill_gen['tau_hltPFdispltau_phi'      ] = tau_gg.hlt_pftau_displ.phi()
            tofill_gen['tau_hltPFdispltau_charge'   ] = tau_gg.hlt_pftau_displ.charge()
            tofill_gen['tau_hltPFdispltau_decaymode'] = tau_gg.hlt_pftau_displ.decayMode()
#             ## new vars
#             tofill_gen['tau_hltPFdispltau_leadChargedCandPdgId'] = tau_gg.hlt_pftau_displ.leadChargedCandPdgId
#             tofill_gen['tau_hltPFdispltau_leadChargedCandPt'   ] = tau_gg.hlt_pftau_displ.leadChargedCandPt
#             tofill_gen['tau_hltPFdispltau_leadPFCandPdgId'     ] = tau_gg.hlt_pftau_displ.leadCandPdgId
#             tofill_gen['tau_hltPFdispltau_leadPFCandPt'        ] = tau_gg.hlt_pftau_displ.leadCandPt
#             tofill_gen['tau_hltPFdispltau_maxHCALPFClusterEt'  ] = tau_gg.hlt_pftau_displ.maxHCALPFClusterEt
#             tofill_gen['tau_hltPFdispltau_nChargedHad'         ] = tau_gg.hlt_pftau_displ.nChargedHad      
#             tofill_gen['tau_hltPFdispltau_nGamma'              ] = tau_gg.hlt_pftau_displ.nGamma 
#             tofill_gen['tau_hltPFdispltau_sumPtCharged'        ] = tau_gg.hlt_pftau_displ.sum_pt_charged 
#             tofill_gen['tau_hltPFdispltau_sumPtNeutral'        ] = tau_gg.hlt_pftau_displ.sum_pt_neutral 
# 
            tofill_gen['tau_hltPFdispltau_dxy'           ] = tau_gg.hlt_pftau_displ.dxy 
            tofill_gen['tau_hltPFdispltau_dxyerr'        ] = tau_gg.hlt_pftau_displ.dxyerr 
            tofill_gen['tau_hltPFdispltau_ip3d'          ] = tau_gg.hlt_pftau_displ.ip3d 
            tofill_gen['tau_hltPFdispltau_ip3derr'       ] = tau_gg.hlt_pftau_displ.ip3derr
            tofill_gen['tau_hltPFdispltau_passChargedIso'] = tau_gg.hlt_pftau_displ.passChargedIso 
#             tofill_gen['tau_hltPFdispltau_passRelChargedIso'] = tau_gg.hlt_pftau_displ.passRelChargedIso 
#             tofill_gen['tau_hltPFdispltau_passAbsChargedIso'] = tau_gg.hlt_pftau_displ.passAbsChargedIso 
#             tofill_gen['tau_hltPFdispltau_isoval']         = tau_gg.hlt_pftau_displ.isoVal 
            tofill_gen['tau_hltPFdispltau_passFilters']         = tau_gg.hlt_pftau_displ.passFilters 

        if hasattr(tau_gg, 'l1_tau') and tau_gg.l1_tau:
            tofill_gen['tau_l1_pt'       ] = tau_gg.l1_tau.pt()
            tofill_gen['tau_l1_eta'      ] = tau_gg.l1_tau.eta()
            tofill_gen['tau_l1_phi'      ] = tau_gg.l1_tau.phi()
            tofill_gen['tau_l1_charge'   ] = tau_gg.l1_tau.charge()
            tofill_gen['tau_l1_iso'      ] = tau_gg.l1_tau.hwIso()

        ntuple_gen.Fill(array('f',tofill_gen.values()))

    ######################################################################################
    # printout some info, if you want
    # printer(taus, gen_taus)

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()