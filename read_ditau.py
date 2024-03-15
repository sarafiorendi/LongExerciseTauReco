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
from treeVariablesHLT_ditau import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer , genDecayModeGEANT, isAncestor# utility functions
from math import sqrt, pow
from copy import deepcopy as dc
from itertools import product as itertools_product

ROOT.gSystem.Load('libRecoTauTagRecoTau')
from ROOT import HLTJetTag
# from ROOT import JetTag

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
    good_gen_status = 2
if 'DY' in sample:  mom_pdgId = [23]

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_ditau_tuple_{}_{}_withDispl_compareToPromptTrigger.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(list(zip(branches, [-9999.]*len(branches)))) # initialise all branches to unphysical -99       
print(outfile_gen.GetName()) 

##########################################################################################
# Get ahold of the events
file_min = 0
file_max = 200
# file_max = 200
f = open('%s_filelist.txt'%args.sample)
infile = f.readlines()[file_min:file_max]

# infile = f.readlines()[ifile]

os.environ['X509_USER_PROXY'] = '/afs/cern.ch/user/f/fiorendi/x509up_u58808' 
redirector = 'root://cms-xrd-global.cern.ch//'
if 'fnal' in sample:  redirector = 'root://cmseos.fnal.gov//'
if 'eos' in sample:   redirector = ''

## to submit on condor
# events = Events(redirector+infile.strip())


files = []
for i in infile:
  files.append(i.strip())
#   files.append(base_filename.replace('_0.', '_%s.'%i))
# filenames = []
events = Events(files)
print ('running locally on n files!!!!!!!!!!!!!!!!!!!!')

# print(infile)
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


def findMatchToTau(hlt_taus, hlt_tau): 
 
    print('tau match')
    hlt_taus_match = hlt_taus 
    for gg in hlt_taus_match: 
        setattr(gg, hlt_tau, None)
        setattr(gg, hlt_tau,0.2)
    return hlt_taus_match    

def printForDebug(cc, gg):
    print('    pt: %.2f'%( cc.pt()), \
          '\t eta: %.2f'%( cc.eta()), \
          '\t phi: %.2f'%( cc.phi()), \
          '\t pdgID: %d'%( cc.pdgId()), \
          '\t dpT/pT = %.2f,  dR = %.2f '%(abs(cc.pt() - gg.vispt())/gg.vispt(), deltaR(cc.eta(),cc.phi(),gg.viseta(), gg.visphi())))

def printGEN(gg):
    print('GEN pt: %.2f'%(gg.vispt(),), \
          '\t eta: %.2f'%(gg.viseta()), \
          '\t phi: %.2f'%(gg.visphi()), '\n')

   
##########################################################################################
# instantiate the handles to the relevant collections.

process_name = 'RECO'
if 'UL' or 'HNL' in sample:  process_name = 'PAT'

# gen particles
handle_gen              = Handle('std::vector<reco::GenParticle>')
# L2 taus
handle_l2_taus          = Handle('std::vector<reco::CaloJet>')
handle_l2_isoTaus       = Handle('std::vector<reco::CaloJet>')
handle_hlt_pftaus       = Handle('std::vector<reco::PFTau>')
handle_hlt_pftaus_displ = Handle('std::vector<reco::PFTau>')

# vertices
handle_vtx              = Handle('std::vector<reco::Vertex>')
# L1 taus
handle_l1               = Handle('BXVector<l1t::Tau>')

handle_IP_displ     = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_iso_displ    = Handle('reco::PFTauDiscriminator')
handle_absiso_displ = Handle('reco::PFTauDiscriminator')
handle_reliso_displ = Handle('reco::PFTauDiscriminator')
handle_isoval_displ = Handle('reco::PFTauDiscriminator')
ip_handle           = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_hlt_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
# handle_hlt_res   = Handle('edm::TriggerResults') 
triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")

handle_isol2        = Handle('reco::JetTagCollection') 
# edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>
handle_tau_for_iso  = Handle('std::vector<reco::CaloJet>') 
# edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>

handles = {}
handles['gen_particles']     = [('prunedGenParticles', '', process_name) , handle_gen, False]
if 'gmsb' in sample:
    handles['gen_particles'] = [('genParticlePlusGeant', '', 'SIM') , handle_gen, False]

# handles['l1_taus']           = [('caloStage2Digis','Tau','RECO'),            handle_l1, False]
handles['l1_taus']           = [('hltGtStage2Digis','Tau','MYHLT'), handle_l1, False]
handles['l2_taus']           = [('hltL2TauJetsL1TauSeeded', '', 'MYHLT'),    handle_l2_taus, False]
handles['l2_isoTaus']        = [('hltL2TauJetsIsoL1TauSeededGlob', '', 'MYHLT'), handle_l2_isoTaus, False]

handles['hlt_l2_tau_for_iso']  = [('hltL2TausForPixelIsolationL1TauSeeded', '', 'MYHLT'), handle_tau_for_iso, False]
handles['hlt_l2iso']           = [('hltL2TauPixelIsoTagProducerL1TauSeededGlob', '', 'MYHLT'), handle_isol2, False]


handles['hlt_pftaus_displ']  = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_pftaus_displ, False]
handles['hlt_IP_displ']   = [('hltHpsPFTauTransverseImpactParameters', '', 'MYHLT'), handle_IP_displ, False]
handles['hlt_iso_displ']  = [('hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator', '', 'MYHLT'), handle_iso_displ, False]
handles['hlt_absiso_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationDiscriminator', '', 'MYHLT'), handle_absiso_displ, False]
handles['hlt_reliso_displ']  = [('hltHpsDisplPFTauMediumRelativeChargedIsolationDiscriminator', '', 'MYHLT'), handle_reliso_displ, False]
handles['hlt_isoval_displ']  = [('hltHpsDisplPFTauMediumAbsoluteChargedIsolationValue', '', 'MYHLT'), handle_isoval_displ, False]

handles['hlt_filter']  = [('hltHpsDoubleMediumChargedIsoDisplPFTau32TrackPt1L1HLTMatchedGlob', '', 'MYHLT'), handle_hlt_filter, False]
 
handles['vtx']               = [('offlineSlimmedPrimaryVertices','','PAT'), handle_vtx, False]

# handles['trig_results'] = [('TriggerResults','','HLT'), handle_hlt_res, False]


if 'raw' in sample:
    handles['gen_particles'] = [('genParticles', '', 'SIM') , handle_gen, False]
    handles['l1_taus']       = [('hltGtStage2Digis','Tau','MYHLT'), handle_l1, False]

# handle_pu       = Handle('std::vector< PileupSummaryInfo>')
# handles['pu']   = [('slimmedAddPileupInfo', '', 'RECO'), handle_pu, False]


met_paths = [
  "HLT_PFMET120_PFMHT120_IDTight", 
  "HLT_PFMET130_PFMHT130_IDTight", 
#   "HLT_PFMET140_PFMHT140_IDTight", 
#   "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", 
  "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", 
  "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", 
#   "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", 
  "HLT_MET105_IsoTrk50", 
  "HLT_MET120_IsoTrk50"
]

lep_paths = [
  "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1", 
  "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130", 
#   "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140", 
#   "HLT_MET105_IsoTrk50", 
#   "HLT_MET120_IsoTrk50"
]


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
#             if k =='hlt_l2iso':  pdb.set_trace()
            ev.getByLabel(v[0], v[1])
            setattr(ev, k, v[1].product())
            v[2] = True
        except:    
            v[2] = False
    

    # select only hadronically decaying taus
    gen_all_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and pp.status()==good_gen_status]
    # select only hadronically decaying taus
    gen_taus = [pp for pp in gen_all_taus if isGenHadTau(pp)]
    if len(gen_taus) < 2:  continue

    # select accepted tau mothers
    gen_moms = [imom for imom in ev.gen_particles if abs(imom.pdgId()) in mom_pdgId]

    ## calculate MET 
    gen_neutrinos = [pp for pp in ev.gen_particles if (abs(pp.pdgId())==12 or abs(pp.pdgId())==16)]
    good_gen_neutrinos = []
    for ineu in gen_neutrinos:
      if abs(ineu.mother(0).pdgId()) == 15:
        good_gen_neutrinos.append(ineu)
                    
    good_gen_lsps = []
    gen_lsps = [pp for pp in ev.gen_particles if abs(pp.pdgId())==1000022 and pp.status()==1]
    for ilsp in gen_lsps:
      if abs(ilsp.mother(0).pdgId()) in mom_pdgId:
        good_gen_lsps.append(ilsp)
    
    gen_met = 0
    for pp in good_gen_lsps+good_gen_neutrinos:
      gen_met = gen_met + pp.pt()

    # info about tau_h
    for gg in gen_taus:
        gg.bestmom = None
        tau_moms = [imom for imom in gen_moms if isAncestor(imom, gg) ]
        if len(tau_moms) > 0:  gg.bestmom = tau_moms[0]

    gen_taus = [gg for gg in gen_taus if gg.bestmom !=None]

    # determine gen decaymode and find mother
    for gg in gen_taus:

        ## reset other attributes
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])
        ### find first dau to be used for vertices 
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print('is none')
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()
        if 'taugun' in sample:
            gg.lxy = sqrt(pow(gg.dau.vx(),2)+pow(gg.dau.vy(),2))
        else:    
            gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))

        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))
        ## now find tau decay point to stau production point distance 
        gg.pi_lxy = sqrt(pow(gg.dau.vx()-gg.bestmom.vx(),2) + pow(gg.dau.vy()-gg.bestmom.vy(),2) )

        if gg.bestmom.pt() != 0:  gg.momct2d = c_const*gg.lxy*gg.bestmom.mass()/gg.bestmom.pt()

        ## calc dR with all other taus/dau leptons in the event
#         gen_other_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
#                           pp.status()==good_gen_status]
    ######################################################################################

    # match L2 taus to gen taus
    if handles['l2_taus'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.l2_taus, 'l2_tau')
 
    # match L2 iso taus to gen taus
    if handles['l2_isoTaus'][2]:  
        gen_taus = findMatchToGen(gen_taus, ev.l2_isoTaus, 'l2_isotau')
        
    if handles['hlt_l2_tau_for_iso'][2] and handles['hlt_l2iso'][2]:
        l2isoval_product = ev.hlt_l2iso
        pt_l2isoval_dict = {}
        for k in range(len(l2isoval_product)):
            pt_l2isoval_dict[l2isoval_product.key(k).pt()] = l2isoval_product.value(k)

        all_l2taus = [pp for pp in ev.hlt_l2_tau_for_iso ]
        for thetau in all_l2taus:  
            setattr(thetau, 'l2isoval', -9999.)

        for itau,thetau in enumerate(all_l2taus):
            try:    thetau.l2isoval            = pt_l2isoval_dict[thetau.pt()]
            except: thetau.l2isoval            =  -9999.

        gen_taus = findMatchToGen(gen_taus, all_l2taus, 'l2_tau_for_iso')

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
                sum_pt_neutral = 0
                for ineu in range(gg.hlt_pftau_displ.nGamma  ) :
                    sum_pt_neutral += gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()

                gg.hlt_pftau_displ.sum_pt_charged =  sum_pt_charged  
                gg.hlt_pftau_displ.sum_pt_neutral =  sum_pt_neutral  

    ######################################################################################
    # access the l1 taus
    if handles['l1_taus'][2]:  
        all_l1taus = []
        for iii in range(ev.l1_taus.size(0)):
            all_l1taus.append(ev.l1_taus.at(0, iii))
        l1_taus = [tau for tau in all_l1taus if tau.pt()>15.]
        gen_taus = findMatchToGen(gen_taus, l1_taus, 'l1_tau',  0.3)

#     ######################################################################################
    pass_met_hlt = 0
    pass_hlts = {}
    pass_lep_hlt = 0
    pass_lep_hlts = {}

    ev.getByLabel(triggerBitLabel, triggerBits)
    names = ev.object().triggerNames(triggerBits.product())
    for i in range(triggerBits.product().size()):   
#         if 'Tau' in names.triggerName(i).split('_v')[0]:  
#           print (names.triggerName(i).split('_v')[0])
        if names.triggerName(i).split('_v')[0] in met_paths:
#           print ("Trigger ", names.triggerName(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)"))
          if triggerBits.product().accept(i) :  
#             print ("Trigger ", names.triggerName(i), ": ", ("PASS" if triggerBits.product().accept(i) else "fail (or not run)"))
            pass_met_hlt = 1
            pass_hlts[names.triggerName(i).split('_v')[0]] = 1
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

        for k, v in tofill_gen.items(): tofill_gen[k] = -9999. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()
        tofill_gen['trueNI'            ] = this_pu
        tofill_gen['gen_met'           ] = gen_met
        tofill_gen['pass_met_hlt'       ] = pass_met_hlt
        tofill_gen['pass_lep_hlt'     ] = pass_lep_hlt
        
        for ipath in met_paths:
            tofill_gen['pass_%s'%ipath.strip('HLT_')] = pass_hlts[ipath]
        for ipath in lep_paths:
            tofill_gen['pass_%s'%ipath.strip('HLT_')] = pass_lep_hlts[ipath]

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
        if hasattr(gg, 'l2_tau_for_iso') and gg.l2_tau_for_iso:
            tofill_gen['tau_l2_for_iso_pt'       ] = gg.l2_tau_for_iso.pt()
            tofill_gen['tau_l2_for_iso_eta'      ] = gg.l2_tau_for_iso.eta()
            tofill_gen['tau_l2_for_iso_phi'      ] = gg.l2_tau_for_iso.phi()
            tofill_gen['tau_l2_for_iso_l2iso'    ] = gg.l2_tau_for_iso.l2isoval
            tofill_gen['tau_l2_for_iso_charge'   ] = gg.l2_tau_for_iso.charge()
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
            
        if hasattr(gg, 'l1_tau') and gg.l1_tau:
            tofill_gen['tau_l1_pt'       ] = gg.l1_tau.pt()
            tofill_gen['tau_l1_eta'      ] = gg.l1_tau.eta()
            tofill_gen['tau_l1_phi'      ] = gg.l1_tau.phi()
            tofill_gen['tau_l1_charge'   ] = gg.l1_tau.charge()
            tofill_gen['tau_l1_iso'      ] = gg.l1_tau.hwIso()

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

        ntuple_gen.Fill(array('f',list(tofill_gen.values())))

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()
