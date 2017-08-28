import ROOT
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

# Get ahold of the events
events = Events('outputFULL.root')
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files

# PAT taus
label_taus  = ('patTaus', '', 'TAURECO')
handle_taus = Handle('std::vector<pat::Tau>')

# gen particles
label_gen  = ('prunedGenParticles', '', 'PAT')
handle_gen = Handle('std::vector<reco::GenParticle>')

# identify hadronically decaying generator level taus
def isGenHadTau(gen_tau):
    ndau = gen_tau.numberOfDaughters()
    isgentau = sum([abs(gen_tau.daughter(dd).pdgId()) in (11,13) for dd in range(ndau)]) == 0
    gen_tau.visp4 = gen_tau.p4()
    for neu in [gen_tau.daughter(dd) for dd in range(ndau) if abs(gen_tau.daughter(dd).pdgId()) in (12,14,16)]:
        gen_tau.visp4 -= neu.p4()
    gen_tau.vispt  = gen_tau.visp4.pt
    gen_tau.viseta = gen_tau.visp4.eta
    gen_tau.visphi = gen_tau.visp4.phi
    return isgentau

# grab the final daughters (after any possible intermediate decay and radiative process)
def finalDaughters(gen, daughters=None):
    if daughters is None:
        daughters = []
    for i in range(gen.numberOfDaughters()):
        daughter = gen.daughter(i)
        if daughter.numberOfDaughters() == 0:
            daughters.append(daughter)
        else:
            finalDaughters(daughter, daughters)
    return daughters

# example histogram
histos = OrderedDict()
histos['pull_pt'] = ROOT.TH1F('pull_pt', 'pull_pt', 50, -2, 2)
histos['pull_pt'].GetXaxis().SetTitle('(p_{T}^{off} - p_{T}^{gen})/p_{T}^{gen}')
histos['pull_pt'].GetYaxis().SetTitle('counts')

# start looping on the events
for i, ev in enumerate(events):
    
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print '===> processing %d / %d event' %(i, totevents)
    
    # access the taus
    ev.getByLabel(label_taus, handle_taus)
    taus = handle_taus.product()
    
    # loosely filter the taus 
    taus = [tau for tau in taus if tau.pt()>18.]

    # get the gen taus
    ev.getByLabel (label_gen, handle_gen)
    gen_particles = handle_gen.product()
    
    # select only hadronically decaying taus
    gen_taus = [pp for pp in gen_particles if abs(pp.pdgId())==15 and pp.status()==2 and isGenHadTau(pp)]

    # determine gen decaymode
    for gg in gen_taus:
        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])

    # match reco taus to gen taus
    for tt in taus    : tt.gen_tau  = None
    for gg in gen_taus: gg.reco_tau = None
    
    gen_taus_copy = gen_taus
    
    for tt in taus:
        matches = [gg for gg in gen_taus_copy if deltaR(tt.p4(), gg.visp4)<0.3]
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(tt.p4(), gg.visp4))
        bestmatch = matches[0]
        tt.gen_tau = bestmatch
        bestmatch.reco_tau = tt
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

    # fill histograms
    for gg in gen_taus:
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            pull = (gg.reco_tau.pt() - gg.vispt())/gg.vispt()
            histos['pull_pt'].Fill(pull)

    

#     if len(gen_taus)>0 and len(taus)>0:
#         import pdb ; pdb.set_trace()
#         print '\n\n===========> RECO TAUS'
#         for tt in taus: 
#             print 'reco tau         pt %.3f eta %.3f phi %.3f' %(tt.pt(), tt.eta(), tt.phi())
#             if hasattr(tt, 'gen_tau') and tt.gen_tau:
#                 print 'matched gen tau  pt %.3f eta %.3f phi %.3f' %(tt.gen_tau.vispt(), tt.gen_tau.viseta(), tt.gen_tau.visphi())
#             else:
#                 print 'matched gen tau  MISSING'
#         print '===========> GEN TAUS'
#         for tt in gen_taus: 
#             print 'gen tau          pt %.3f eta %.3f phi %.3f' %(tt.vispt(), tt.viseta(), tt.visphi())
#             if hasattr(tt, 'reco_tau') and tt.reco_tau:
#                 print 'matched reco tau pt %.3f eta %.3f phi %.3f' %(tt.reco_tau.pt(), tt.reco_tau.eta(), tt.reco_tau.phi())
#             else:
#                 print 'matched reco tau MISSING'
