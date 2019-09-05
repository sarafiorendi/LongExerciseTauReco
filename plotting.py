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

## SOLUTION
# Exercise 3: Calculate the overall efficiency of the tau selection.
# Simple loop - This could also be done in less lines with tree.Draw()
if sample=="ZTT":
	total = 0
	passed = 0
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		if (tree.tau_gen_vis_pt > 18.0 and abs(tree.tau_gen_eta)<2.4):  # we are only interested in taus with sufficient pT in the central detector region
			total += 1
			if tree.tau_reco_pt>0.0: # > 0 means that is exists i.e. has been reconstructed as we initialized the array with -99
				passed += 1
	print("{} taus in total.".format(total))
	print("{} were reconstructed as taus.".format(passed))
	print("Efficieny = {}".format(round(float(passed)/float(total),4)))

## SOLUTION
# Exercise 4: Calculate the overall misidentification rate.
# Simple loop - This could also be done in less lines with tree.Draw()
if sample=="QCD":
	total = 0
	passed = 0
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		if (tree.jet_pt > 18.0 and abs(tree.jet_eta)<2.4):  # we are only interested in taus with sufficient pT in the central detector region
			total += 1
			if tree.tau_reco_pt>0.0: # > 0 means that is exists i.e. has been reconstructed as we initialized the array with -99
				passed += 1
	print("{} jets in total.".format(total))
	print("{} were reconstructed as taus.".format(passed))
	print("Misidentification rate = {}".format(round(float(passed)/float(total),4)))


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

## SOLUTION
# Exercise 6: Do the same for jets misidentified as tau leptons
# efficiency example using TEfficiency https://root.cern.ch/doc/master/classTEfficiency.html
# reconstruction efficiency as a funciton of jet vis tau pt
bins = np.array([18., 20., 22., 24., 28., 32., 36., 40., 45., 50., 60., 80., 100.]) # variable binning
histo_den = ROOT.TH1F('den2', 'den2', len(bins)-1, bins)
histo_num = ROOT.TH1F('num2', 'num2', len(bins)-1, bins)

# apply minimal cuts on the denominator in order to select only taus that 
# fall in the detector acceptance (pt > 18 GeV and abs(eta)<2.4, tracker coverage)
tree.Draw('jet_pt >> den2', 'jet_pt>18 & abs(jet_eta)<2.4')
tree.Draw('jet_pt >> num2', 'jet_pt>18 & abs(jet_eta)<2.4 & tau_reco_pt>0')
c = ROOT.TCanvas("c","c")
c.cd()
eff = ROOT.TEfficiency(histo_num, histo_den)
eff.SetTitle(';Jet p_{T} [GeV]; reconstruction efficiency')
eff.SetMarkerStyle(8)
eff.Draw('APL')

# save the current canvas
c.SaveAs('jet_efficiency_{}.pdf'.format(sample))
c.SaveAs('jet_efficiency_{}.png'.format(sample))
c.Clear()

## SOLUTION
# Exercise 7: Plot true vs. reconstructed decay mode
# define the histogram(s)
dm_gen = ROOT.TH1F('dm_gen', 'dm_gen', 11, -0.5, 10.5)
dm_reco = ROOT.TH1F('dm_reco', 'dm_reco', 11, -0.5, 10.5)

# Same as before: Fill the histograms from the TTree
tree.Draw('tau_gen_decaymode >> dm_gen', 'tau_gen_vis_pt>18 & abs(tau_gen_vis_eta)<2.4 & tau_reco_pt>0')
tree.Draw('tau_reco_decaymode >> dm_reco', 'tau_gen_vis_pt>18 & abs(tau_gen_vis_eta)<2.4 & tau_reco_pt>0')

# Create canvas
c = ROOT.TCanvas("c","c")
c.cd()

# Set different colors for the two histograms
dm_gen.SetLineColor(1)
dm_reco.SetLineColor(2)
dm_gen.GetXaxis().SetTitle('Decay Mode')
dm_gen.GetYaxis().SetTitle('counts')
dm_gen.Draw()
dm_reco.Draw("same")

# Ensure that all data points fit on the plot.
dm_gen.SetAxisRange(0, 1.1*max(dm_gen.GetMaximum(),dm_reco.GetMaximum()), "Y")
dm_reco.SetAxisRange(0, 1.1*max(dm_gen.GetMaximum(),dm_reco.GetMaximum()), "Y")

# Create legend and place on coordinates (x1,y1,x2,y2)
leg = ROOT.TLegend(0.4,0.6,0.8,0.8)
# No box around legend
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.AddEntry(dm_gen,"gen level","L")
leg.AddEntry(dm_reco,"reco level","L")
leg.Draw()

# save the current canvas
c.SaveAs('dm_{}.pdf'.format(sample))
c.SaveAs('dm_{}.png'.format(sample))
c.Clear()

# SOLUTION: Exercise 9
# Add '& (tree.)tau_reco_tight_id' to lines 34, 50, 84, 108
