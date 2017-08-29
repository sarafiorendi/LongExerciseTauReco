'''
Define here utility functions to be imported in the main FWLite analyser.
'''

# identify hadronically decaying generator level taus
def isGenHadTau(gen_tau):
    ndau = gen_tau.numberOfDaughters()
    # identify as had decaying if no electron/muon appears as daughter
    isgentau = sum([abs(gen_tau.daughter(dd).pdgId()) in (11,13) for dd in range(ndau)]) == 0
    # assign the gen tau a new attribut, visp4, 
    # which is the visible 4-momentum of the generated tau.
    # This is done by subtracting the neutrino p4's from the total gen tau p4.
    gen_tau.visp4 = gen_tau.p4()
    for neu in [gen_tau.daughter(dd) for dd in range(ndau) if abs(gen_tau.daughter(dd).pdgId()) in (12,14,16)]:
        gen_tau.visp4 -= neu.p4()
    gen_tau.vispt   = gen_tau.visp4.pt
    gen_tau.viseta  = gen_tau.visp4.eta
    gen_tau.visphi  = gen_tau.visp4.phi
    gen_tau.vismass = gen_tau.visp4.mass
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

# print tau reco/gen information in a nicely laid out format
def printer(taus, gen_taus):
    print '\n\n===========> RECO TAUS'
    for tt in taus: 
        print 'reco tau         pt %.3f eta %.3f phi %.3f' %(tt.pt(), tt.eta(), tt.phi())
        if hasattr(tt, 'gen_tau') and tt.gen_tau:
            print 'matched gen tau  pt %.3f eta %.3f phi %.3f' %(tt.gen_tau.vispt(), tt.gen_tau.viseta(), tt.gen_tau.visphi())
        else:
            print 'matched gen tau  MISSING'
    print '===========> GEN TAUS'
    for gg in gen_taus: 
        print 'gen tau          pt %.3f eta %.3f phi %.3f' %(gg.vispt(), gg.viseta(), gg.visphi())
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            print 'matched reco tau pt %.3f eta %.3f phi %.3f' %(gg.reco_tau.pt(), gg.reco_tau.eta(), gg.reco_tau.phi())
        else:
            print 'matched reco tau MISSING'
