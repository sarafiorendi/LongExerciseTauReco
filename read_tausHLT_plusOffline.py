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
from treeVariablesHLT import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer , genDecayModeGEANT, isAncestor# utility functions
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
outfile_gen = ROOT.TFile('tau_gentau_tuple_{}_{}_withDispl.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(zip(branches, [-9999.]*len(branches))) # initialise all branches to unphysical -99       
print outfile_gen.GetName() 

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
# events = Events('/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_1_0/src/HLTrigger/Configuration/test/outputHLT_2.root')
# events = Events('/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/iter4/CMSSW_12_0_0_pre6/src/outputHLT_2.root')
# events = Events('/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/iter4/CMSSW_12_0_0_pre6/src/outputHLT_2.root')
# events = Events('/eos/cms/store/user/fiorendi/displacedTaus/Staus_M_500_10mm_14TeV_Run3MC/MINIAODSIM/210616_092721/0000/TSG-Run3Winter21DRMiniAOD-staus_M_500_10mm_inMINIAODSIM_467.root') # make sure this corresponds to your file name!
# print 'using a local file since fnal is down!!!!!!!!!!!!!!!!!!!!!'

print infile
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
#             if hlt_tau=='hlt_pftau_displ':  pdb.set_trace()
#             print 'found one', hlt_tau
    return gen_taus_match    

def findMatchToTau(hlt_taus, hlt_tau): 
 
    print 'tau match'
    hlt_taus_match = hlt_taus 
    for gg in hlt_taus_match: 
        setattr(gg, hlt_tau, None)
        setattr(gg, hlt_tau,0.2)
    return hlt_taus_match    

def printForDebug(cc, gg):
    print '    pt: %.2f'%( cc.pt()), \
          '\t eta: %.2f'%( cc.eta()), \
          '\t phi: %.2f'%( cc.phi()), \
          '\t pdgID: %d'%( cc.pdgId()), \
          '\t dpT/pT = %.2f,  dR = %.2f '%(abs(cc.pt() - gg.vispt())/gg.vispt(), deltaR(cc.eta(),cc.phi(),gg.viseta(), gg.visphi()))

def printGEN(gg):
    print 'GEN pt: %.2f'%(gg.vispt(),), \
          '\t eta: %.2f'%(gg.viseta()), \
          '\t phi: %.2f'%(gg.visphi()), '\n'

   
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
handle_hlt_pftaus_step0 = Handle('std::vector<reco::PFTau>')
handle_hlt_pftaus_step1 = Handle('std::vector<reco::PFTau>')
handle_hlt_pftaus_step2 = Handle('std::vector<reco::PFTau>')

# PF
handle_hltPFjets        = Handle('std::vector<reco::PFJet>')
handle_hltPFs           = Handle('std::vector<reco::PFCandidate>')
## tracks 
handle_hltTracks        = Handle('std::vector<reco::Track>')
handle_hltIter4Tracks   = Handle('std::vector<reco::Track>')
handle_hltIter04Tracks  = Handle('std::vector<reco::Track>')
# vertices
handle_vtx              = Handle('std::vector<reco::Vertex>')
# L1 taus
handle_l1               = Handle('BXVector<l1t::Tau>')
#  offline tracks 
handle_lost_tracks      = Handle('std::vector<pat::PackedCandidate>')
handle_packed           = Handle('std::vector<pat::PackedCandidate>')

handle_IP_displ     = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_iso_displ    = Handle('reco::PFTauDiscriminator')
handle_absiso_displ = Handle('reco::PFTauDiscriminator')
handle_reliso_displ = Handle('reco::PFTauDiscriminator')
handle_isoval_displ = Handle('reco::PFTauDiscriminator')
ip_handle           = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_hlt_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 


handles = {}
handles['gen_particles']     = [('prunedGenParticles', '', process_name) , handle_gen, False]
if 'gmsb' in sample:
    handles['gen_particles'] = [('genParticlePlusGeant', '', 'SIM') , handle_gen, False]

handles['l1_taus']           = [('caloStage2Digis','Tau','RECO'),            handle_l1, False]
handles['l2_taus']           = [('hltL2TauJetsL1TauSeeded', '', 'MYHLT'),    handle_l2_taus, False]
handles['l2_isoTaus']        = [('hltL2TauJetsIsoL1TauSeeded', '', 'MYHLT'), handle_l2_isoTaus, False]

handles['hlt_pfs']           = [('hltParticleFlowForTaus', '', 'MYHLT'),                           handle_hltPFs, False]
handles['hlt_pftaus']        = [('hltHpsL1JetsHLTSinglePFTauTrackPt1MatchGlob', '', 'MYHLT'),      handle_hlt_pftaus, False]
# handles['hlt_pftaus_displ']  = [('hltHpsSelectedPFTausTrackPt1GlobDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]
handles['hlt_pftaus_displ']  = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]

handles['hlt_pftaus_step0']  = [('hltHpsCombinatoricRecoTausDispl', '', 'MYHLT'), handle_hlt_pftaus_step0, False]
handles['hlt_pftaus_step1']  = [('hltHpsPFTauProducerSansRefsDispl', '', 'MYHLT'), handle_hlt_pftaus_step1, False]
handles['hlt_pftaus_step2']  = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_pftaus_step2, False]

handles['hlt_IP_displ']   = [('hltHpsPFTauTransverseImpactParameters', '', 'MYHLT'), handle_IP_displ, False]
handles['hlt_iso_displ']  = [('hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator', '', 'MYHLT'), handle_iso_displ, False]
handles['hlt_absiso_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationDiscriminator', '', 'MYHLT'), handle_absiso_displ, False]
handles['hlt_reliso_displ']  = [('hltHpsDisplPFTauMediumRelativeChargedIsolationDiscriminator', '', 'MYHLT'), handle_reliso_displ, False]
handles['hlt_isoval_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationValue', '', 'MYHLT'), handle_isoval_displ, False]

handles['hlt_filter']  = [('hltHpsSinglePFTau32TrackPt1L1HLTMatchedGlobDispl', '', 'MYHLT'), handle_hlt_filter, False]


## mycheck
# print 'change!!!!!!!!!!!!!!!!!!!!!!!!!!'
# handles['hlt_pftaus']        = [('hltHpsSelectedPFTausTrackPt1NoIsolationGlobDispl', '', 'MYHLT'),      handle_hlt_pftaus, False]
# handles['hlt_pftaus_displ']  = [('hltHpsL1JetsHLTDoublePFTauTrackPt1MatchGlobDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]


## test tau steps
# handle_hlt_pftaus_displ = Handle('std::vector<reco::PFTau>')
# handles['hlt_pftaus_displ']  = [('hltHpsL1JetsHLTDoublePFTauTrackPt1MatchGlobDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]
# 
 
handles['hltPFjets']         = [('hltAK4PFJetsReg', '', 'MYHLT'),                 handle_hltPFjets, False]

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

# handle_pu       = Handle('std::vector< PileupSummaryInfo>')
# handles['pu']   = [('slimmedAddPileupInfo', '', 'RECO'), handle_pu, False]




##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print '===> processing %d / %d event' %(i, totevents)
    
    for k, v in handles.iteritems():
        setattr(ev, k, None)
        v[2] = False
        try:
            ev.getByLabel(v[0], v[1])
            setattr(ev, k, v[1].product())
            v[2] = True
        except:    
            v[2] = False
            
    # select only hadronically decaying taus
    gen_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                pp.status()==good_gen_status and isGenHadTau(pp)]
    if len(gen_taus) == 0:  continue

    # determine gen decaymode and find mother
    for gg in gen_taus:

        ## reset other attributes
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])
        if 'gmsb' in sample:
            dm_string = genDecayModeGEANT( [d for d in finalDaughters(gg) if abs(d.pdgId()) not in [12, 14, 16]])
            gg.decayMode = tauDecayModes.nameToInt(dm_string)
        
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
        if gg.dau == None:  print 'is none'
#         else: pdb.set_trace()

        if gg.bestmom == None or gg.dau == None:
            continue
        

#     for gg in gen_taus:
#         print 
#         pdb.set_trace()
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

        ## gg = tau   bestmom = stau    vx = production point 
        ## now find tau decay point to stau production point distance 
        gg.pi_lxy = sqrt(pow(gg.dau.vx()-gg.bestmom.vx(),2) + pow(gg.dau.vy()-gg.bestmom.vy(),2) )

#             stau_vectorP = np.array([bestmom.px(), bestmom.py(), 0])
#             stau_vectorL = np.array([gg.vx()-bestmom.vx(), gg.vy()-bestmom.vy(), 0])
#             gg.tau_cosxy  = tau_vectorL.dot(tau_vectorP)/((np.linalg.norm(tau_vectorL) * np.linalg.norm(tau_vectorP)))
#             l3d = sqrt(pow(gg.vx()-bestmom.vx(),2) + pow(gg.vy()-bestmom.vy(),2) + pow(gg.vz()-bestmom.vz(),2))
#             if bestmom.p() != 0: gg.momct    = c_const*l3d*bestmom.mass()/bestmom.p()
        if gg.bestmom.pt() != 0:  gg.momct2d = c_const*gg.lxy*gg.bestmom.mass()/gg.bestmom.pt()
        
        ## calc dR with all other taus/dau leptons in the event
        gen_other_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                          pp.status()==good_gen_status]
#         for iothertau in gen_other_taus:
#           daus = finalDaughters(iothertau)
#           for idaus in daus:
#             print idaus.pdgId(), deltaR(idaus.p4(), gg.visp4)
#             deltaR(idaus.p4(), gg.visp4)
#           print '\n'                            
#         pdb.set_trace()                  
        
        
        
    ######################################################################################

    if handles['lost_tracks'][2] and  handles['packed'][2]:  
    	lt = [hh for hh in ev.lost_tracks if hh.pt() > 10]
    	pc = [hh for hh in ev.packed if hh.pt() > 10]
        alltks = lt + pc

        for gg in gen_taus :
            gg.offtk = None
            offtks_skim   = [cc for cc in alltks if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.p4(), gg.visp4)<0.4]
            bestcom = bestMatch(gg, offtks_skim )
            if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
                gg.offtk = bestcom[0]

#     if handles['reco_taus'][2]:  
#         gen_taus = findMatchToGen(gen_taus, ev.reco_taus, 'reco_tau')
# #         for recotau in ev.reco_taus:
# #             if recotau.leadChargedHadrCand().isNonnull():
# # #                 pdb.set_trace()
# #                 mytrk=recotau.leadChargedHadrCand().get()
# #                 print mytrk.pseudoTrack().originalAlgo(), mytrk.pseudoTrack().algo()
# # #             if recotau.leadTauChargedHadronCandidate().isNonnull():
# # #                 pdb.set_trace()
# 
# 
    # match L2 taus to gen taus
    if handles['l2_taus'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.l2_taus, 'l2_tau')
 
    # match L2 iso taus to gen taus
    if handles['l2_isoTaus'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.l2_isoTaus, 'l2_isotau')

# # #     print '---------------- matching pf jets -------------------'
# #     if handles['hltPFjets'][2]:  
# #         gen_taus, hltPFjets = findMatchToGen(gen_taus, ev.hltPFjets, 'hltPFjet')
# # 
    # match hlt PF taus to gen taus
    if handles['hlt_pftaus'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.hlt_pftaus, 'hlt_pftau')

    if handles['hlt_pftaus_step0'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.hlt_pftaus_step0, 'hlt_pftau_step0')
        
    if handles['hlt_pftaus_step1'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.hlt_pftaus_step1, 'hlt_pftau_step1')

    if handles['hlt_pftaus_step2'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.hlt_pftaus_step2, 'hlt_pftau_step2')

        ## printing stuff
#         for gg in gen_taus:
#             if hasattr(gg, 'hlt_pftau_step2') and gg.hlt_pftau_step2:
#                 print '---------------------- step 2 ---------------------------------'
#                 print  gg.vispt()
# #                 pdb.set_trace()
#                 for ich in range(gg.hlt_pftau_step2.signalChargedHadrCands().size()  ) :
#                     print '\t ch: ', gg.hlt_pftau_step2.signalChargedHadrCands()[ich].pdgId(), '\t', gg.hlt_pftau_step2.signalChargedHadrCands()[ich].pt()
#                 for ineu in range(gg.hlt_pftau_step2.signalGammaCands().size()   ) :
#                     print '\t neu: ', gg.hlt_pftau_step2.signalGammaCands()[ineu].pdgId(), '\t', gg.hlt_pftau_step2.signalGammaCands()[ineu].pt()


    # match hlt displced PF taus to gen taus
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

        if handles['hlt_reliso_displ'][2]:
            reliso_product = ev.hlt_reliso_displ
            pt_reliso_dict = {}
            for k in range(len(reliso_product)):
                pt_reliso_dict[reliso_product.key(k).pt()] = reliso_product.value(k)

        if handles['hlt_absiso_displ'][2]:
            absiso_product = ev.hlt_absiso_displ
            pt_absiso_dict = {}
            for k in range(len(absiso_product)):
                pt_absiso_dict[absiso_product.key(k).pt()] = absiso_product.value(k)

        if handles['hlt_isoval_displ'][2]:
            isoval_product = ev.hlt_isoval_displ
            pt_isoval_dict = {}
            for k in range(len(isoval_product)):
                pt_isoval_dict[isoval_product.key(k).pt()] = isoval_product.value(k)

        if handles['hlt_filter'][2]:
            filter_product = ev.hlt_filter.pftauRefs()
            pt_filter_list = []
            for k in range(len(filter_product)):
                pt_filter_list.append(filter_product[k].pt()) 

            for itau,thetau in enumerate(all_taus):
                try:    thetau.passChargedIso    = pt_iso_dict[thetau.pt()]
                except: thetau.passChargedIso    =  -9999.
    
                try:    thetau.passRelChargedIso = pt_reliso_dict[thetau.pt()]
                except: thetau.passRelChargedIso =  -9999.
    
                try:    thetau.passAbsChargedIso = pt_absiso_dict[thetau.pt()]
                except: thetau.passAbsChargedIso =  -9999.
    
                try:    thetau.isoVal            = pt_isoval_dict[thetau.pt()]
                except: thetau.isoVal            =  -9999.
    
                if thetau.pt() in pt_filter_list: thetau.passFilters = 1
                else:                             thetau.passFilters = 0

        if handles['hlt_IP_displ'][2]:
            ip_product = ev.hlt_IP_displ
            pt_ip_dict = {}
            for k in range(len(ip_product)):
                pt_ip_dict[ip_product.key(k).pt()] = [ip_product.value(k).dxy(), ip_product.value(k).dxy_error(), 
                                                      ip_product.value(k).ip3d(), ip_product.value(k).ip3d_error()] 

            for itau,thetau in enumerate(all_taus):
                try:  
                    thetau.dxy     = pt_ip_dict[thetau.pt()][0]
                    thetau.dxyerr  = pt_ip_dict[thetau.pt()][1]
                    thetau.ip3d    = pt_ip_dict[thetau.pt()][2]
                    thetau.ip3derr = pt_ip_dict[thetau.pt()][3]
#                     print thetau.dxy 
                except:  pass
        

        gen_taus = findMatchToGen(gen_taus, all_taus, 'hlt_pftau_displ')

        for gg in gen_taus:
            if hasattr(gg, 'hlt_pftau_displ') and gg.hlt_pftau_displ:

                theLeadChargedCand = gg.hlt_pftau_displ.leadChargedHadrCand()
                theLeadPFCand      = gg.hlt_pftau_displ.leadCand()
                gg.hlt_pftau_displ.leadChargedCandPt     = theLeadChargedCand.pt()
                gg.hlt_pftau_displ.leadChargedCandPdgId  = theLeadChargedCand.pdgId()
                gg.hlt_pftau_displ.leadCandPt            = theLeadPFCand.pt()
                gg.hlt_pftau_displ.leadCandPdgId         = theLeadPFCand.pdgId()

                gg.hlt_pftau_displ.maxHCALPFClusterEt    = gg.hlt_pftau_displ.maximumHCALPFClusterEt()
                gg.hlt_pftau_displ.nChargedHad           = gg.hlt_pftau_displ.signalChargedHadrCands().size() 
                gg.hlt_pftau_displ.nGamma                = gg.hlt_pftau_displ.signalGammaCands().size() 

                sum_pt_charged = 0
                for ich in range(gg.hlt_pftau_displ.nChargedHad  ) :
                    sum_pt_charged += gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()
#                     print '\t ch: ', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pdgId(), '\t', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()
    
                sum_pt_neutral = 0
                for ineu in range(gg.hlt_pftau_displ.nGamma  ) :
                    sum_pt_neutral += gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()
#                     print '\t neu: ', gg.hlt_pftau_displ.signalGammaCands()[ineu].pdgId(), '\t', gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()

                gg.hlt_pftau_displ.sum_pt_charged =  sum_pt_charged  
                gg.hlt_pftau_displ.sum_pt_neutral =  sum_pt_neutral  




    if handles['hlt_pfs'][2]:  
    	hltPFs_all   = [hh for hh in ev.hlt_pfs if hh.pt() > 10]
    	hltPFs_pions = [hh for hh in ev.hlt_pfs if abs(hh.pdgId())==211 and hh.pt() > 10]

        for gg in gen_taus :
            gg.hlt_pf = None; gg.hlt_pf_pi = None
            hltPFs_all_skim   = [cc for cc in hltPFs_all   if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.p4(), gg.visp4)<0.4]
            hltPFs_pions_skim = [cc for cc in hltPFs_pions if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.p4(), gg.visp4)<0.4]

            bestPF = bestMatch(gg, hltPFs_all_skim )
            if bestPF[0] != None and sqrt(bestPF[1]) < dR_cone :
                gg.hlt_pf = bestPF[0]
            bestPFpi = bestMatch(gg, hltPFs_pions_skim )
            if bestPFpi[0] != None and sqrt(bestPFpi[1]) < dR_cone :
                gg.hlt_pf_pi = bestPFpi[0]


    ######################### tracks
    # match hlt tracks to gen taus
    if handles['hltTks'][2]:  
        for gg in gen_taus :
            gg.hltTk = None 
            hltTks_skim  = [cc for cc in ev.hltTks if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.eta(),cc.phi(),gg.viseta(),gg.visphi())<0.4]
            bestcom = bestMatch(gg, hltTks_skim )
            if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
                gg.hltTk = bestcom[0]

    if handles['hltIter4Tks'][2]:  
        for gg in gen_taus :
            gg.hltIt4Tk = None 
            hltIt4Tks_skim  = [cc for cc in ev.hltIter4Tks if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.eta(),cc.phi(),gg.viseta(),gg.visphi())<0.4]
            bestcom = bestMatch(gg, hltIt4Tks_skim )
            if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
                gg.hltIt4Tk = bestcom[0]

    if handles['hltIter04Tks'][2]:  
        for gg in gen_taus :
            gg.hltIt04Tk = None 
            hltIt04Tks_skim  = [cc for cc in ev.hltIter04Tks if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.eta(),cc.phi(),gg.viseta(),gg.visphi())<0.4]
            bestcom = bestMatch(gg, hltIt04Tks_skim )
            if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
                gg.hltIt04Tk = bestcom[0]

    if handles['hltMuMergedTks'][2]:  
        for gg in gen_taus :
            gg.hltMuMTk = None 
            hltMuMTks_skim  = [cc for cc in ev.hltMuMergedTks if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and deltaR(cc.eta(),cc.phi(),gg.viseta(),gg.visphi())<0.4]
            bestcom = bestMatch(gg, hltMuMTks_skim )
            if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
                gg.hltMuMTk = bestcom[0]


    ######################################################################################
    # access the l1 taus
    l1_taus = [tau for tau in ev.l1_taus if ev.l1_taus.getFirstBX()==0]
    l1_taus = [tau for tau in l1_taus if tau.pt()>15.]
#     pdb.set_trace()

    gen_taus = findMatchToGen(gen_taus, l1_taus, 'l1_tau')

    # match reco taus to gen taus
#     for l1 in l1_taus : l1.gen_tau  = None # first initialise the matching to None
#     for gg in gen_taus: gg.l1_tau   = None # first initialise the matching to None
#     
#     gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched
#     
#     for l1 in l1_taus:
#         matches = [gg for gg in gen_taus_copy if deltaR(l1.p4(), gg.visp4)<dR_cone]
#         ave_matches = len(matches)
# #         if ave_matches > 1:  print ave_matches
#         if not len(matches):
#             continue
#         matches.sort(key = lambda gg : deltaR(l1.p4(), gg.visp4))
#         bestmatch = matches[0]
#         l1.gen_tau = bestmatch
#         bestmatch.l1_tau = l1
#         gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

#     ######################################################################################
    # Save pileup information in MC #
    this_pu = -99.
#     bx_vector = []
# #     for ipuinfo in ev.pu:
# #         if ipuinfo.getBunchCrossing() == 0:
# #             this_pu = ipuinfo.getTrueNumInteractions()
# #             break            
# 
# 
#     ######################################################################################
    
    # fill the ntuple: each gen tau makes an entry
    for gg in gen_taus:
        if gg.bestmom == None or gg.dau == None:  continue

        for k, v in tofill_gen.iteritems(): tofill_gen[k] = -9999. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()
        tofill_gen['trueNI'            ] = this_pu
# #         tofill_gen['nvtx'              ] = vertices.size()
# #         tofill_gen['PV_x'              ] = vertices[0].x()
# #         tofill_gen['PV_y'              ] = vertices[0].y()
# #         tofill_gen['PV_z'              ] = vertices[0].z()
        if hasattr(gg, 'l2_tau') and gg.l2_tau:
            tofill_gen['tau_L2_mass'     ] = gg.l2_tau.mass()
            tofill_gen['tau_L2_pt'       ] = gg.l2_tau.pt()
            tofill_gen['tau_L2_eta'      ] = gg.l2_tau.eta()
            tofill_gen['tau_L2_phi'      ] = gg.l2_tau.phi()
            tofill_gen['tau_L2_charge'   ] = gg.l2_tau.charge()
        if hasattr(gg, 'l2_isotau') and gg.l2_isotau:
            tofill_gen['tau_isoL2_mass'     ] = gg.l2_isotau.mass()
            tofill_gen['tau_isoL2_pt'       ] = gg.l2_isotau.pt()
            tofill_gen['tau_isoL2_eta'      ] = gg.l2_isotau.eta()
            tofill_gen['tau_isoL2_phi'      ] = gg.l2_isotau.phi()
            tofill_gen['tau_isoL2_charge'   ] = gg.l2_isotau.charge()
#         if hasattr(gg, 'hltPFjet') and gg.hltPFjet:
#             tofill_gen['tau_hltPFjet_mass'     ] = gg.hltPFjet.mass()
#             tofill_gen['tau_hltPFjet_pt'       ] = gg.hltPFjet.pt()
#             tofill_gen['tau_hltPFjet_eta'      ] = gg.hltPFjet.eta()
#             tofill_gen['tau_hltPFjet_phi'      ] = gg.hltPFjet.phi()
#             tofill_gen['tau_hltPFjet_charge'   ] = gg.hltPFjet.charge()
# 
        if hasattr(gg, 'hlt_pftau') and gg.hlt_pftau:
            tofill_gen['tau_hltPFtau_mass'     ] = gg.hlt_pftau.mass()
            tofill_gen['tau_hltPFtau_pt'       ] = gg.hlt_pftau.pt()
            tofill_gen['tau_hltPFtau_eta'      ] = gg.hlt_pftau.eta()
            tofill_gen['tau_hltPFtau_phi'      ] = gg.hlt_pftau.phi()
            tofill_gen['tau_hltPFtau_charge'   ] = gg.hlt_pftau.charge()
            tofill_gen['tau_hltPFtau_decaymode'] = gg.hlt_pftau.decayMode()
        if hasattr(gg, 'hlt_pftau_displ') and gg.hlt_pftau_displ:
            tofill_gen['tau_hltPFdispltau_mass'     ] = gg.hlt_pftau_displ.mass()
            tofill_gen['tau_hltPFdispltau_pt'       ] = gg.hlt_pftau_displ.pt()
            tofill_gen['tau_hltPFdispltau_eta'      ] = gg.hlt_pftau_displ.eta()
            tofill_gen['tau_hltPFdispltau_phi'      ] = gg.hlt_pftau_displ.phi()
            tofill_gen['tau_hltPFdispltau_charge'   ] = gg.hlt_pftau_displ.charge()
            tofill_gen['tau_hltPFdispltau_decaymode'] = gg.hlt_pftau_displ.decayMode()
            ## new vars
            tofill_gen['tau_hltPFdispltau_leadChargedCandPdgId'] = gg.hlt_pftau_displ.leadChargedCandPdgId
            tofill_gen['tau_hltPFdispltau_leadChargedCandPt'   ] = gg.hlt_pftau_displ.leadChargedCandPt
            tofill_gen['tau_hltPFdispltau_leadPFCandPdgId'     ] = gg.hlt_pftau_displ.leadCandPdgId
            tofill_gen['tau_hltPFdispltau_leadPFCandPt'        ] = gg.hlt_pftau_displ.leadCandPt
            tofill_gen['tau_hltPFdispltau_maxHCALPFClusterEt'  ] = gg.hlt_pftau_displ.maxHCALPFClusterEt
            tofill_gen['tau_hltPFdispltau_nChargedHad'         ] = gg.hlt_pftau_displ.nChargedHad      
            tofill_gen['tau_hltPFdispltau_nGamma'              ] = gg.hlt_pftau_displ.nGamma 
            tofill_gen['tau_hltPFdispltau_sumPtCharged'        ] = gg.hlt_pftau_displ.sum_pt_charged 
            tofill_gen['tau_hltPFdispltau_sumPtNeutral'        ] = gg.hlt_pftau_displ.sum_pt_neutral 

            tofill_gen['tau_hltPFdispltau_dxy'           ] = gg.hlt_pftau_displ.dxy 
            tofill_gen['tau_hltPFdispltau_dxyerr'        ] = gg.hlt_pftau_displ.dxyerr 
            tofill_gen['tau_hltPFdispltau_ip3d'          ] = gg.hlt_pftau_displ.ip3d 
            tofill_gen['tau_hltPFdispltau_ip3derr'       ] = gg.hlt_pftau_displ.ip3derr
            tofill_gen['tau_hltPFdispltau_passChargedIso'] = gg.hlt_pftau_displ.passChargedIso 
            tofill_gen['tau_hltPFdispltau_passRelChargedIso'] = gg.hlt_pftau_displ.passRelChargedIso 
            tofill_gen['tau_hltPFdispltau_passAbsChargedIso'] = gg.hlt_pftau_displ.passAbsChargedIso 
            tofill_gen['tau_hltPFdispltau_isoval']         = gg.hlt_pftau_displ.isoVal 
            tofill_gen['tau_hltPFdispltau_passFilters']         = gg.hlt_pftau_displ.passFilters 
            

        if hasattr(gg, 'hlt_pftau_step0') and gg.hlt_pftau_step0:
            tofill_gen['tau_hltPFdispltau_0_pt'       ] = gg.hlt_pftau_step0.pt()
            tofill_gen['tau_hltPFdispltau_0_eta'      ] = gg.hlt_pftau_step0.eta()
            tofill_gen['tau_hltPFdispltau_0_phi'      ] = gg.hlt_pftau_step0.phi()
            tofill_gen['tau_hltPFdispltau_0_charge'   ] = gg.hlt_pftau_step0.charge()
            tofill_gen['tau_hltPFdispltau_0_decaymode'] = gg.hlt_pftau_step0.decayMode()
        if hasattr(gg, 'hlt_pftau_step1') and gg.hlt_pftau_step1:
            tofill_gen['tau_hltPFdispltau_1_pt'       ] = gg.hlt_pftau_step1.pt()
            tofill_gen['tau_hltPFdispltau_1_eta'      ] = gg.hlt_pftau_step1.eta()
            tofill_gen['tau_hltPFdispltau_1_phi'      ] = gg.hlt_pftau_step1.phi()
            tofill_gen['tau_hltPFdispltau_1_charge'   ] = gg.hlt_pftau_step1.charge()
            tofill_gen['tau_hltPFdispltau_1_decaymode'] = gg.hlt_pftau_step1.decayMode()
        if hasattr(gg, 'hlt_pftau_step2') and gg.hlt_pftau_step2:
            tofill_gen['tau_hltPFdispltau_2_pt'       ] = gg.hlt_pftau_step2.pt()
            tofill_gen['tau_hltPFdispltau_2_eta'      ] = gg.hlt_pftau_step2.eta()
            tofill_gen['tau_hltPFdispltau_2_phi'      ] = gg.hlt_pftau_step2.phi()
            tofill_gen['tau_hltPFdispltau_2_charge'   ] = gg.hlt_pftau_step2.charge()
            tofill_gen['tau_hltPFdispltau_2_decaymode'] = gg.hlt_pftau_step2.decayMode()
             
            

        if hasattr(gg, 'hlt_pf') and gg.hlt_pf:
            tofill_gen['tau_hltPF_pt'       ] = gg.hlt_pf.pt()
            tofill_gen['tau_hltPF_eta'      ] = gg.hlt_pf.eta()
            tofill_gen['tau_hltPF_phi'      ] = gg.hlt_pf.phi()
            tofill_gen['tau_hltPF_charge'   ] = gg.hlt_pf.charge()
            tofill_gen['tau_hltPF_pdgid'    ] = gg.hlt_pf.pdgId()
        if hasattr(gg, 'hlt_pf_pi') and gg.hlt_pf_pi:
            tofill_gen['tau_hltPFpi_pt'       ] = gg.hlt_pf_pi.pt()
            tofill_gen['tau_hltPFpi_eta'      ] = gg.hlt_pf_pi.eta()
            tofill_gen['tau_hltPFpi_phi'      ] = gg.hlt_pf_pi.phi()
            tofill_gen['tau_hltPFpi_charge'   ] = gg.hlt_pf_pi.charge()
            
        if hasattr(gg, 'hltTk') and gg.hltTk:
            tofill_gen['tau_hltTk_pt'       ] = gg.hltTk.pt()
            tofill_gen['tau_hltTk_eta'      ] = gg.hltTk.eta()
            tofill_gen['tau_hltTk_phi'      ] = gg.hltTk.phi()
            tofill_gen['tau_hltTk_charge'   ] = gg.hltTk.charge()
        if hasattr(gg, 'hltIt4Tk') and gg.hltIt4Tk:
            tofill_gen['tau_hltIt4Tk_pt'       ] = gg.hltIt4Tk.pt()
            tofill_gen['tau_hltIt4Tk_eta'      ] = gg.hltIt4Tk.eta()
            tofill_gen['tau_hltIt4Tk_phi'      ] = gg.hltIt4Tk.phi()
            tofill_gen['tau_hltIt4Tk_charge'   ] = gg.hltIt4Tk.charge()
        if hasattr(gg, 'hltIt04Tk') and gg.hltIt04Tk:
            tofill_gen['tau_hltIt04Tk_pt'       ] = gg.hltIt04Tk.pt()
            tofill_gen['tau_hltIt04Tk_eta'      ] = gg.hltIt04Tk.eta()
            tofill_gen['tau_hltIt04Tk_phi'      ] = gg.hltIt04Tk.phi()
            tofill_gen['tau_hltIt04Tk_charge'   ] = gg.hltIt04Tk.charge()

        if hasattr(gg, 'l1_tau') and gg.l1_tau:
            tofill_gen['tau_l1_pt'       ] = gg.l1_tau.pt()
            tofill_gen['tau_l1_eta'      ] = gg.l1_tau.eta()
            tofill_gen['tau_l1_phi'      ] = gg.l1_tau.phi()
            tofill_gen['tau_l1_charge'   ] = gg.l1_tau.charge()
            tofill_gen['tau_l1_iso'      ] = gg.l1_tau.hwIso()

            
        if hasattr(gg, 'hltMuMTk') and gg.hltMuMTk:
            tofill_gen['tau_hltMuMTk_pt'       ] = gg.hltMuMTk.pt()
            tofill_gen['tau_hltMuMTk_eta'      ] = gg.hltMuMTk.eta()
            tofill_gen['tau_hltMuMTk_phi'      ] = gg.hltMuMTk.phi()
            tofill_gen['tau_hltMuMTk_charge'   ] = gg.hltMuMTk.charge()
# 
        if hasattr(gg, 'offtk') and gg.offtk:
            tofill_gen['tau_offtk_pt'       ] = gg.offtk.pt()
            tofill_gen['tau_offtk_eta'      ] = gg.offtk.eta()
            tofill_gen['tau_offtk_phi'      ] = gg.offtk.phi()
            tofill_gen['tau_offtk_charge'   ] = gg.offtk.charge()
#         if hasattr(gg, 'reco_tau') and gg.reco_tau:
#             tofill_gen['tau_reco_tau_mass'     ] = gg.reco_tau.mass()
#             tofill_gen['tau_reco_tau_pt'       ] = gg.reco_tau.pt()
#             tofill_gen['tau_reco_tau_eta'      ] = gg.reco_tau.eta()
#             tofill_gen['tau_reco_tau_phi'      ] = gg.reco_tau.phi()
#             tofill_gen['tau_reco_tau_charge'   ] = gg.reco_tau.charge()
#             tofill_gen['tau_reco_tau_decaymode'] = gg.reco_tau.decayMode()
#             
#         if hasattr(gg, 'l1_jet') and gg.l1_jet:
# #             tofill_gen['tau_l1_mass'     ] = gg.reco_tau.mass()
#             tofill_gen['jet_l1_pt'       ] = gg.l1_jet.pt()
#             tofill_gen['jet_l1_eta'      ] = gg.l1_jet.eta()
#             tofill_gen['jet_l1_phi'      ] = gg.l1_jet.phi()
#             tofill_gen['jet_l1_charge'   ] = gg.l1_jet.charge()
# 
        tofill_gen['tau_gen_ndau'      ] = gg.numberOfDaughters()
        tofill_gen['tau_gen_pt'        ] = gg.pt()
        tofill_gen['tau_gen_eta'       ] = gg.eta()
        tofill_gen['tau_gen_phi'       ] = gg.phi()
        tofill_gen['tau_gen_charge'    ] = gg.charge()
        tofill_gen['tau_gen_decaymode' ] = gg.decayMode
        tofill_gen['tau_gen_lxy'       ] = gg.lxy
        tofill_gen['tau_gen_dxy'       ] = gg.dxy
        tofill_gen['tau_gen_vx'        ] = gg.bestmom.vx()  ## user defined ones (mother prod. vertex)
        tofill_gen['tau_gen_vy'        ] = gg.bestmom.vy()
        tofill_gen['tau_gen_vz'        ] = gg.bestmom.vz()
        tofill_gen['tau_gen_cosxy'     ] = gg.cosxy
#         tofill_gen['tau_gen_momct'     ] = gg.momct
        tofill_gen['tau_gen_mom_mass'  ] = gg.bestmom.mass()
        tofill_gen['tau_gen_momct2d'   ] = gg.momct2d
        tofill_gen['tau_gen_vis_mass'  ] = gg.vismass()
        tofill_gen['tau_gen_vis_pt'    ] = gg.vispt()
        tofill_gen['tau_gen_vis_eta'   ] = gg.viseta()
        tofill_gen['tau_gen_vis_phi'   ] = gg.visphi()

        tofill_gen['tau_gen_pi_lxy'  ] = gg.pi_lxy  ## ok
        tofill_gen['tau_gen_pi_cosxy'] = gg.pi_cosxy

        ntuple_gen.Fill(array('f',tofill_gen.values()))

    ######################################################################################
    # printout some info, if you want
    # printer(taus, gen_taus)

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()








'''
#         for gg in gen_taus :
#             hltPFs_all_skim = []
#             hltPFs_pions_skim = []
# #             if gg.lxy > 10:
# #                 print '--------------' 
# #                 printGEN(gg) 
#             
#             for cc in hltPFs_pions:
# #                 if gg.lxy > 10:  printForDebug(cc, gg) 
#                 if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and ( deltaR(cc.p4(), gg.visp4)<0.4) :
#                     hltPFs_pions_skim.append(cc)
#             bestcom = bestMatch(gg, hltPFs_pions_skim )
#             if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
#                 gg.hlt_pf_pi = bestcom[0]
# #                 if gg.lxy > 10:
# #                     print 'good pfpion has pt: %.2f and pdgID %f' %(gg.hlt_pf_pi.pt(), abs(gg.hlt_pf_pi.pdgId()))
# 
# #             if hasattr(gg, 'hlt_pftau') and gg.hlt_pftau:
# #                 print 'hps pt: %.2f'%( gg.hlt_pftau.pt()), \
# #                       '\t eta: %.2f'%( gg.hlt_pftau.eta()), \
# #                       '\t phi: %.2f'%( gg.hlt_pftau.phi())
#                     
#             for cc in hltPFs_all:
# #                 if gg.lxy > 10:  printForDebug(cc, gg) 
#                 if (abs(cc.pt() - gg.vispt())<dPt_cone*gg.vispt()) and ( deltaR(cc.p4(), gg.visp4)<0.4) :
#                     hltPFs_all_skim.append(cc)
#             bestcom = bestMatch(gg, hltPFs_all_skim )
#             if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
#                 gg.hlt_pf = bestcom[0]
# #                 if gg.lxy > 10:
# #                     print 'good pfcand has pt: %.2f and pdgID %f' %(gg.hlt_pf.pt(), abs(gg.hlt_pf.pdgId()))
# 
# 
'''
