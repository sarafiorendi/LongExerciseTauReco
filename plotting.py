'''
Example to reade the ntuples and produce plots
'''

import ROOT

# open the input file and get the tree 
infile = ROOT.TFile.Open('tau_gen_tuple.root', 'read')
infile.cd()
tree = infile.Get('tree')

# define the histogram(s)
histo_pull = ROOT.TH1F('pull', 'pull', 40, -1, 1)
histo_pull.GetXaxis().SetTitle('(p_{T}^{off} - p_{T}^{gen})/p_{T}^{gen}')
histo_pull.GetYaxis().SetTitle('counts')

# draw the histogram
# require that the generator level pt to be larger than 15 GeV and that 
# a good reconstructed tau is properly matched to the generator level tau
tree.Draw('(tau_reco_pt - tau_gen_vis_pt)/tau_gen_vis_pt >> pull', 'tau_gen_vis_pt>15 & tau_reco_pt>0')

# save the current pad in pdf format
ROOT.gPad.SaveAs('pull.pdf')
