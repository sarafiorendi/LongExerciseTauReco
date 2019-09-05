'''
Example to read the ntuples and produce plots
'''

import ROOT
import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser(description="Plot quantities from ntuples")
parser.add_argument("--file",
		choices=['ZTT','QCD'],
		required=True,
		help='Specify the sample you want to use for plotting')
args = parser.parse_args()
sample = args.file


# open the input file and get the tree 
infile = ROOT.TFile.Open('tau_{}_tuple_{}.root'.format("jet" if sample=="QCD" else "gentau", sample), 'read')
infile.cd()
tree = infile.Get('tree')

# define the histogram(s)
histo_pull = ROOT.TH1F('pull', 'pull', 40, -1, 1)
histo_pull.GetXaxis().SetTitle('(p_{T}^{off} - p_{T}^{gen})/p_{T}^{gen}')
histo_pull.GetYaxis().SetTitle('counts')

# draw the histogram
# require that the generator level pt to be larger than 18 GeV and that 
# a good reconstructed tau is properly matched to the generator level tau
tree.Draw('(tau_reco_pt - tau_gen_vis_pt)/tau_gen_vis_pt >> pull', 'tau_gen_vis_pt>18 & abs(tau_gen_vis_eta)<2.4 & tau_reco_pt>0')
c = ROOT.TCanvas("c","c")
c.cd()
histo_pull.Draw()
# save the current canvas
c.SaveAs('pull_{}.pdf'.format(sample))
c.SaveAs('pull_{}.png'.format(sample))
c.Clear()


# efficiency example using TEfficiency https://root.cern.ch/doc/master/classTEfficiency.html
# reconstruction efficiency as a funciton of gen vis tau pt
bins = np.array([18., 20., 22., 24., 28., 32., 36., 40., 45., 50., 60., 80., 100.]) # variable binning
histo_den = ROOT.TH1F('den', 'den', len(bins)-1, bins)
histo_num = ROOT.TH1F('num', 'num', len(bins)-1, bins)

# apply minimal cuts on the denominator in order to select only taus that 
# fall in the detector acceptance (pt > 18 GeV and abs(eta)<2.4, tracker coverage)
tree.Draw('tau_gen_vis_pt >> den', 'tau_gen_vis_pt>18 & abs(tau_gen_vis_eta)<2.4')
tree.Draw('tau_gen_vis_pt >> num', 'tau_gen_vis_pt>18 & abs(tau_gen_vis_eta)<2.4 & tau_reco_pt>0')
c = ROOT.TCanvas("c","c")
c.cd()
eff = ROOT.TEfficiency(histo_num, histo_den)
eff.SetTitle(';generator-level #tau p_{T}^{vis} [GeV]; reconstruction efficiency')
eff.SetMarkerStyle(8)
eff.Draw('APL')

# save the current canvas
c.SaveAs('reco_efficiency_{}.pdf'.format(sample))
c.SaveAs('reco_efficiency_{}.png'.format(sample))
c.Clear()