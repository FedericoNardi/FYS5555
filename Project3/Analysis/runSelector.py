from ROOT import *
import sys, os, time

t0 = time.time()

arg1 = sys.argv[1]

if arg1 == 'Data':
  input_dir = '../2lep/Data/'
elif arg1 == 'MC':
  input_dir = '../2lep/MC/BSM_Signal_Samples/'
  input_dir1 = '../2lep/MC/SM_Backgrounds/'

#Wmass = 4200
#Nmass = 2100

myChain = TChain('mini')


for filename in os.listdir(input_dir):
  if not '.root' in filename: continue
  if arg1 == 'MC':
    if not 'LRSM' in filename: continue
    #if not str(Wmass) in filename: continue
    #if not str(Nmass) in filename: continue
  print filename
  myChain.Add(input_dir+filename)
if arg1 == 'MC':
  for filename in os.listdir(input_dir1):
    if not '.root' in filename: continue
    print filename
    myChain.Add(input_dir1+filename)


entries = myChain.GetEntries()

print "-------------------------------------------"
if arg1 == 'Data':
  print "Running on real data!"
else:
  print "Running on Monte Carlo!"
print "Number of events to process: %d" %entries
print "-------------------------------------------"

if arg1 == 'Data':
  myChain.Process("MySelector.C", "Data")
else:
  myChain.Process("MySelector.C", "MC")

t = int( time.time()-t0 )/60

print "Time spent: %d min" %t 
