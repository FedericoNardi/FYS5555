import ROOT
from ROOT import *
import os, sys
import teststat as stat
import numpy as np

from infofile import infos

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetGridStyle(2)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.TH1.AddDirectory(kFALSE)

channel = sys.argv[1]

# Variables

#variables = ['pt1', 'pt2', 'eta1', 'eta2', 'phi1', 'phi2', 'mll', 'met']
variables = ['mll', 'mlljj', 'pt']

#xtitles = {'pt1':'Leading lepton p_{T} (GeV)', 'pt2':'Subleading lepton p_{T} (GeV)', 'eta1':'Leading lepton #eta', 'eta2':'Subleading lepton #eta', 'phi1':'Leading lepton #phi', 'phi2':'Subleading lepton #phi', 'mll':'m_{ll} (GeV)', 'met':'E_{T}^{miss} (GeV)'}
xtitles = {'mll':'m_{ll} (GeV)', 'mlljj':'m_{lljj} (GeV)','pt':'Total transverse momentum p_{T}^{tot}'}


# Backgrounds

backgrounds = ['Zjets', 'Top', 'Diboson', 'Wjets']

signals = ['M_W3000N1500', 'M_W2400N1800']
#signals = ['M_W5000N2000', 'M_W4200N2100', 'M_W3600N1800', 'M_W3000N1500', 'M_W2400N1800']

Zjets = [364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108, 364109, 364110, 364111, 364112, 364113, 364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122, 364123, 364124,
         364125, 364126, 364127, 364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136, 364137, 364138, 364139, 364140, 364141]

Wjets = [364156, 364157, 364158, 364159, 364160, 364161, 364162, 364163, 364164, 364165, 364166, 364167, 364168, 364169, 364170, 364171, 364172, 364173, 364174, 364175, 364176, 364177, 364178, 364179, 364180,
         364181, 364182, 364183, 364184, 364185, 364186, 364187, 364188, 364189, 364190, 364191, 364192, 364193, 364194, 364195, 364196, 364197]

Diboson = [361600, 361601, 361602, 361603, 361604, 361606, 361607, 361609, 361610]

Top = [410000, 410011, 410012, 4100013, 410014, 410025, 410026]

fileIDs = {'Diboson':Diboson, 'Zjets':Zjets, 'Wjets':Wjets, 'Top':Top}
fileIDs_signal = {'M_W5000N2000':[309070], 'M_W4200N2100':[302708], 'M_W3600N1800':[302701], 'M_W3000N1500':[302687], 'M_W2400N1800':[302681]}

hist_bkg = {}
for var in variables:
  hist_bkg[var] = {}
  for bkg in backgrounds:
    hist_bkg[var][bkg] = TH1F()

hist_sign = {}
for var in variables:
    hist_sign[var]={}
    for sig in signals:
      hist_sign[var][sig]=TH1F()

colours = dict(Diboson=kAzure+1, Top=kRed+1, Zjets=kOrange-2, Wjets=kGray)
colours_sign = dict(M_W5000N2000=kBlue, M_W4200N2100=kBlue+1, M_W3600N1800=kBlue+2, M_W3000N1500=kBlue+3, M_W2400N1800=kBlue+4)


# Extract info about cross section and sum of weights from infofile

info = {}
for key in infos.keys():
  ID = infos[key]['DSID']
  info[ID] = {}
  info[ID]['xsec'] = infos[key]['xsec']
  info[ID]['sumw'] = infos[key]['sumw']
  info[ID]['events'] = infos[key]['events']


# Function for making histograms

L = 10.6 # integrated luminosity = 10.6 fb^-1

def fill_hist(h, h_name, key, ID):

  h_midl = infile.Get(h_name).Clone("h_midl")

  xsec = 1000*info[ID]['xsec']
  nev = info[ID]['sumw']

  N_mc = xsec*L

  sf = N_mc/nev

  if not h.GetName():
    h=infile.Get(h_name)
    h.Scale(sf)
    n = h.GetNbinsX()
    for i in range(n):
      bc = h.GetBinContent(i)
      if bc < 0:
        h.SetBinContent(i,0)
    if 'M_' in key:
      h.SetLineColor(colours_sign[key])
      h.SetLineWidth(2)
      h.SetLineStyle(1)
    else:
      h.SetFillColor(colours[key])
      h.SetLineColor(colours[key])
  else:
    h_midl.Scale(sf)
    n = h_midl.GetNbinsX()
    for i in range(n):
      bc = h_midl.GetBinContent(i)
      if bc < 0:
        h_midl.SetBinContent(i,0)
    h.Add(h_midl)


  return h


# Loop over files in MC directory

for filename in os.listdir('Histograms/MC/'):
  if '.root' in filename:
    filepath = 'Histograms/MC/'+filename
    infile = TFile(filepath)
    file_id = int(filename.split('.')[2])
    #print filename
    for var in variables:
      for bkg in backgrounds:
        if file_id in fileIDs[bkg]:
          hist_bkg[var][bkg] = fill_hist(hist_bkg[var][bkg], 'h_'+channel+'_'+var, bkg, file_id)
      for sig in signals:
        if file_id in fileIDs_signal[sig]:
          hist_sign[var][sig] = fill_hist(hist_sign[var][sig], 'h_'+channel+'_'+var, sig, file_id)

# Get data

data = TFile('Histograms/Data/hist.Data.2016.root')
hist_d ={}

for var in variables:
  hist_d[var] = data.Get('h_'+channel+'_'+var)
  hist_d[var].SetMarkerStyle(20)
  hist_d[var].SetMarkerSize(0.7)
  hist_d[var].SetLineColor(kBlack)
  hist_d[var].GetYaxis().SetTitle("Events")
  hist_d[var].GetXaxis().SetTitle(xtitles[var])
  hist_d[var].GetXaxis().SetTitleFont(43)
  hist_d[var].GetXaxis().SetTitleSize(16)
  hist_d[var].GetYaxis().SetTitleFont(43)
  hist_d[var].GetYaxis().SetTitleSize(16)
  hist_d[var].GetXaxis().SetLabelFont(43)
  hist_d[var].GetXaxis().SetLabelSize(16)
  hist_d[var].GetYaxis().SetLabelFont(43)
  hist_d[var].GetYaxis().SetLabelSize(16)
  hist_d[var].GetXaxis().SetTitleOffset(4)
  hist_d[var].GetYaxis().SetTitleOffset(1.5)


# Style histograms, make stack and histograms with full background

stack = {}
stack_sign = {}
hist_r = {}
hist_mc = {}

for var in variables:
  stack[var] = THStack(var, "")
  stack_sign[var] = THStack(var,"")
  hist_mc[var] = TH1F()
  hist_r[var] = TH1F()
  for sig in reversed(signals):
    hist_sign[var][sig].GetYaxis().SetTitle("Events")
    hist_sign[var][sig].GetXaxis().SetTitle(xtitles[var])
    hist_sign[var][sig].GetXaxis().SetTitleFont(43)
    hist_sign[var][sig].GetXaxis().SetTitleSize(16)
    hist_sign[var][sig].GetYaxis().SetTitleFont(43)
    hist_sign[var][sig].GetYaxis().SetTitleSize(16)
    hist_sign[var][sig].GetXaxis().SetLabelFont(43)
    hist_sign[var][sig].GetXaxis().SetLabelSize(16)
    hist_sign[var][sig].GetYaxis().SetLabelFont(43)
    hist_sign[var][sig].GetYaxis().SetLabelSize(16)
    hist_sign[var][sig].GetXaxis().SetTitleOffset(4)
    hist_sign[var][sig].GetYaxis().SetTitleOffset(1.5)
    #stack_sign[var].Add(hist_sign[var][sig])
  for bkg in reversed(backgrounds):
    hist_bkg[var][bkg].GetYaxis().SetTitle("Events")
    hist_bkg[var][bkg].GetXaxis().SetTitle(xtitles[var])
    hist_bkg[var][bkg].GetXaxis().SetTitleFont(43)
    hist_bkg[var][bkg].GetXaxis().SetTitleSize(16)
    hist_bkg[var][bkg].GetYaxis().SetTitleFont(43)
    hist_bkg[var][bkg].GetYaxis().SetTitleSize(16)
    hist_bkg[var][bkg].GetXaxis().SetLabelFont(43)
    hist_bkg[var][bkg].GetXaxis().SetLabelSize(16)
    hist_bkg[var][bkg].GetYaxis().SetLabelFont(43)
    hist_bkg[var][bkg].GetYaxis().SetLabelSize(16)
    hist_bkg[var][bkg].GetXaxis().SetTitleOffset(4)
    hist_bkg[var][bkg].GetYaxis().SetTitleOffset(1.5)
    stack[var].Add(hist_bkg[var][bkg])
    if not hist_mc[var].GetName():
      hist_mc[var] = hist_bkg[var][bkg].Clone()
    else:
      hist_mc[var].Add(hist_bkg[var][bkg])
    hist_r[var] = hist_d[var].Clone()
    hist_r[var].Divide(hist_mc[var])
    hist_r[var].SetTitle("")
    hist_r[var].GetXaxis().SetTitle(xtitles[var])
    hist_r[var].GetYaxis().SetTitle("Data/#SigmaMC")
    hist_r[var].GetYaxis().SetNdivisions(506)
    hist_r[var].SetMarkerStyle(20)
    hist_r[var].SetMarkerSize(0.7)


# Make plot legend

leg = TLegend(0.70,0.50,0.88,0.88)
leg.SetFillStyle(4000)
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetBorderSize(0)

bkg_labels = {'Zjets':'Z+jets', 'Top':'Top', 'Diboson':'Diboson', 'Wjets':'W+jets'}
signal_labels = {'M_W3000N1500':'W(3000)N(1500)', 'M_W2400N1800':'W(2400)N(1800)'}
#signal_labels = {'M_W5000N2000':'W(5000)N(2000)', 'M_W4200N2100':'W(4200)N(2100)', 'M_W3600N1800':'W(3600)N(1800)', 'M_W3000N1500':'W(3000)N(1500)', 'M_W2400N1800':'W(2400)N(1800)'}


for bkg in backgrounds:
  leg.AddEntry(hist_bkg['mll'][bkg], bkg_labels[bkg], "f")
for sig in signals:
  leg.AddEntry(hist_sign['mll'][sig], signal_labels[sig],'f')
leg.AddEntry(hist_d['mll'],"Data","ple")

selection = ""
if channel == "ee":
  selection = "ee"
if channel == "uu":
  selection = "#mu#mu"

# Make plots

for var in variables:

  cnv = TCanvas("cnv_"+var,"transparent pad", 500, 500)
  cnv.SetTicks(1,1)
  cnv.SetLeftMargin(0.13)
  #cnv.SetLogy()

  p1 = TPad("p1", "", 0, 0.35, 1, 1)
  p2 = TPad("p2", "", 0, 0.0, 1, 0.35)
  p2.SetFillStyle(4000)

  p1.SetLogy()
  p1.SetBottomMargin(0.0)
  p1.Draw()
  p1.cd()

  stack[var].Draw("hist")
  stack[var].SetMinimum(10E-2)
  stack[var].GetYaxis().SetTitle("Events")
  stack[var].GetYaxis().SetTitleFont(43)
  stack[var].GetYaxis().SetTitleSize(16)
  stack[var].GetYaxis().SetLabelFont(43)
  stack[var].GetYaxis().SetLabelSize(16)
  stack[var].GetYaxis().SetTitleOffset(1.5)
  if var in ['eta1', 'eta2', 'phi1', 'phi2']:
    maximum = stack[var].GetMaximum()
    stack[var].SetMaximum(maximum*10E4)

  for sig in signals:
    hist_sign[var][sig].Draw("same hist")

  hist_d[var].Draw("same e0")
  leg.Draw("same")

  s = TLatex()
  s.SetNDC(1);
  s.SetTextAlign(13);
  s.SetTextColor(kBlack);
  s.SetTextSize(0.044);
  s.DrawLatex(0.4,0.86,"#font[72]{ATLAS} Open Data");
  s.DrawLatex(0.4,0.81,"#bf{#sqrt{s} = 13 TeV,^{}%.1f^{}fb^{-1}}" % (L));
  s.DrawLatex(0.4,0.76,"#bf{"+selection+" selection}");


  p1.Update()
  p1.RedrawAxis()

  cnv.cd()

  p2.Draw()
  p2.cd()

  p2.SetGridy()

  hist_r[var].SetMaximum(1.99)
  hist_r[var].SetMinimum(0.01)
  hist_r[var].Draw("0PZ")

  p2.SetTopMargin(0)
  p2.SetBottomMargin(0.35)
  p2.Update()
  p2.RedrawAxis()

  cnv.cd()
  cnv.Update()
  cnv.Print('Plots/'+channel+'_'+var+'.png')
  cnv.Close()

  for sig in signals:
    stat.TestStat(hist_d[var],stack[var],hist_sign[var][sig],var,sig)
