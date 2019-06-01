import ROOT
from ROOT import *
import os, sys
import numpy as np


# Function for test statistic and significance

def TestStat(Data, MC, Signal, var, sig):

  cvs = TCanvas("teststatic","ts",500,500)
  Nbins = Data.GetNbinsX()
  Zg = TGraph(Nbins)

  for i in range(1,Nbins+1):
    x = Data.GetXaxis().GetBinCenter(i)
    n = Data.GetBinContent(i)
    b = MC.GetStack().Last().GetBinContent(i)
    s = Signal.GetBinContent(i)
    logL_b = TMath.Log(TMath.Poisson(n,b))
    logL_sb = TMath.Log(TMath.Poisson(n,s+b))
    # Profile likelihood:
    Q = -2*(logL_b-logL_sb)
    # fix negative values
    if Q<0 or Q==float("inf"):
      Q=0
    #print 'Q: '+ str(Q)
    Z = np.sqrt(Q)
    Zg.SetPoint(i,x,Z)

  Zg.Draw("AC")

  Zg.GetYaxis().SetTitle("significance Z")
  Zg.GetYaxis().SetTitleFont(43)
  Zg.GetYaxis().SetTitleSize(16)
  Zg.GetYaxis().SetLabelFont(43)
  Zg.GetYaxis().SetLabelSize(16)
  Zg.GetYaxis().SetTitleOffset(1.5)
  xtitles = {'mll':'m_{ll} (GeV)', 'mlljj':'m_{lljj} (GeV)','pt':'Total transverse momentum p_{T}^{tot}'}
  Zg.GetXaxis().SetTitle(xtitles[var])
  Zg.GetXaxis().SetTitleFont(43)
  Zg.GetXaxis().SetTitleSize(16)
  Zg.GetXaxis().SetLabelFont(43)
  Zg.GetXaxis().SetLabelSize(16)
  Zg.GetXaxis().SetTitleOffset(1.5)
  signal = {'M_W3000N1500':'m_{W}=3000 GeV, m_{N}=1500 GeV', 'M_W2400N1800':'m_{W}=2400, m_{N}=1800'}
  Zg.SetTitle(signal[sig])

  cvs.Print("Plots/teststat_"+var+"_"+sig+".png")
  cvs.Close()
