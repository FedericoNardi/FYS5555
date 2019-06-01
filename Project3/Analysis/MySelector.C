#define MySelector_cxx
// The class definition in MySelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("MySelector.C")
// root> T->Process("MySelector.C","some options")
// root> T->Process("MySelector.C+")
//


#include "MySelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include <map>
#include <math.h>

// Define som global variables

TString option;
TString channel;

vector<TString> channels;

int events_tot = 0;
int events_ee = 0;
int events_uu = 0;

Int_t DSID = 0;
Int_t prev_DSID = 0;

TLorentzVector l1;
TLorentzVector l2;
TLorentzVector jet1;
TLorentzVector jet2;
TLorentzVector finalstates;
TLorentzVector dileptons;

bool TwoJets = true;
bool SameSign = true;

Float_t wgt;

void MySelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  option = GetOption();

  channels = {"ee", "uu"}; // Channels to be considered

  // Initialize all histograms
  for( const auto & chn:channels ){
    // Normal histograms
    //h_pt1[chn] = new TH1D("h_"+chn+"_pt1", chn+"_pt1", 30, 0, 1500);
    //h_pt2[chn] = new TH1D("h_"+chn+"_pt2", chn+"_pt2", 30, 0, 1500);
    //h_eta1[chn] = new TH1D("h_"+chn+"_eta1", chn+"_eta1", 40, -3, 3);
    //h_eta2[chn] = new TH1D("h_"+chn+"_eta2", chn+"_eta2", 40, -3, 3);
    //h_phi1[chn] = new TH1D("h_"+chn+"_phi1", chn+"_phi1", 40, -M_PI, M_PI);
    //h_phi2[chn] = new TH1D("h_"+chn+"_phi2", chn+"_phi2", 40, -M_PI, M_PI);
    h_pt[chn] = new TH1D("h_"+chn+"_pt", chn+"_pt", 10, 0 , 1500);
    h_mlljj[chn] = new TH1D("h_"+chn+"_mlljj", chn+"_mlljj", 10, 200, 1800);
    h_mll[chn] = new TH1D("h_"+chn+"_mll", chn+"_mll", 30, 0, 2000);
  }


}

void MySelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  option = GetOption();

}

Bool_t MySelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);

  // Count events
  events_tot++;
  if( events_tot % 1000000 == 0 ){ cout << events_tot/1E6 << " million events processed" << endl; }

  // Write histograms to file if we enter a new MC DSID (DSID = DataSet ID)
  if( option == "MC" ){
    DSID = *channelNumber;
    if(prev_DSID == 0){ prev_DSID = *channelNumber; }
    if( DSID != prev_DSID ){ WriteToFile(Form("%d", prev_DSID), option); }
    prev_DSID = DSID;
  }


  // Do event selection

  // Cut #1: Require (exactly) 2 leptons
  if(*lep_n != 2){ return kTRUE; }
  // Cut #2: Require opposite charge
  if(!SameSign){
    if(lep_charge[0] == lep_charge[1]){ return kTRUE; }
  }
  else{
    if(lep_charge[0] != lep_charge[1]){ return kTRUE; }
  }
  // Cut #3: Require same flavour (2 electrons or 2 muons)
  if(lep_type[0] != lep_type[1]){ return kTRUE; }

  // Cut #4: require only 2 jets
  if(TwoJets){
    if(*jet_n!=2){return kTRUE;}
  }

  // Cut #5: require high momentum on jets
  if(TwoJets){
    if(jet_pt[0]/1000. < 30.){return kTRUE;}
    if(jet_pt[1]/1000. < 30.){return kTRUE;}
  }


  // Identify the leptons
  if(fabs(lep_type[0])==11){ channel = "ee"; events_ee++; } // Electrons
  if(fabs(lep_type[0])==13){ channel = "uu"; events_uu++; } // Muons

  if(channel == "ee"){
    if(lep_pt[0]/1000. < 50){return kTRUE;}
    if(lep_pt[1]/1000. < 50){return kTRUE;}
    if(*trigE != 1){return kTRUE;}
  }

  if(channel == "uu"){
    if(lep_pt[0]/1000. < 50){return kTRUE;}
    if(lep_pt[1]/1000. < 50){return kTRUE;}
    if(*trigM != 1){return kTRUE;}
  }

  // Cut #6: Cut on rapidity considering detector geometry
  if(channel=="ee") {
    if (fabs(lep_eta[0]>2.47)||(fabs(lep_eta[0])>1.37 && fabs(lep_eta[0])<1.52)){return kTRUE;}
    if (fabs(lep_eta[0]>2.47)||(fabs(lep_eta[1])>1.37 && fabs(lep_eta[1])<1.52)){return kTRUE;}
  }
  if(channel=="uu"){
    if(fabs(lep_eta[0])>2.5){return kTRUE;}
    if(fabs(lep_eta[1])>2.5){return kTRUE;}
  }

  // Set Lorentz vectors:
  if(lep_pt[0]>lep_pt[1]){
    l1.SetPtEtaPhiE(lep_pt[0]/1000., lep_eta[0], lep_phi[0], lep_E[0]/1000.);
    l2.SetPtEtaPhiE(lep_pt[1]/1000., lep_eta[1], lep_phi[1], lep_E[1]/1000.);
  }
  if(lep_pt[1]>lep_pt[0]){
    l1.SetPtEtaPhiE(lep_pt[1]/1000., lep_eta[1], lep_phi[1], lep_E[1]/1000.);
    l2.SetPtEtaPhiE(lep_pt[0]/1000., lep_eta[0], lep_phi[0], lep_E[0]/1000.);
  }
  if(TwoJets) {
    if (jet_pt[0] > jet_pt[1]) {
      l1.SetPtEtaPhiE(jet_pt[0] / 1000., jet_eta[0], jet_phi[0], jet_E[0] / 1000.);
      l2.SetPtEtaPhiE(jet_pt[1] / 1000., jet_eta[1], jet_phi[1], jet_E[1] / 1000.);
    }
    if (lep_pt[1] > lep_pt[0]) {
      l1.SetPtEtaPhiE(jet_pt[1] / 1000., jet_eta[1], jet_phi[1], jet_E[1] / 1000.);
      l2.SetPtEtaPhiE(jet_pt[0] / 1000., jet_eta[0], jet_phi[0], jet_E[0] / 1000.);
    }
  }
  // Variables are stored in the TTree with unit MeV, so we need to divide by 1000
  // to get GeV, which is a more practical and commonly used unit.

  // Cut the Z peak

  dileptons = l1 + l2;
  if(fabs(dileptons.M() - 91.0) < 20){return kTRUE;}

  finalstates = dileptons + jet1 + jet2;
  if(fabs(finalstates.M() - 200)<20){return kTRUE;}



  // Calculate event weight
  wgt = 1.0;
  if( option == "MC" ){
    wgt = (*mcWeight)*(*scaleFactor_PILEUP)*(*scaleFactor_ELE)*(*scaleFactor_MUON)*(*scaleFactor_BTAG)*(*scaleFactor_LepTRIGGER);
  }

  // Fill histograms

  //h_pt1[channel]->Fill(l1.Pt(), wgt);
  //h_pt2[channel]->Fill(l2.Pt(), wgt);
  //h_eta1[channel]->Fill(l1.Eta(), wgt);
  //h_eta2[channel]->Fill(l2.Eta(), wgt);
  //h_phi1[channel]->Fill(l1.Phi(), wgt);
  //h_phi2[channel]->Fill(l2.Phi(), wgt);
  h_pt[channel]->Fill(finalstates.Pt(), wgt);
  //h_met[channel]->Fill(*met_et/1000.0, wgt);
  h_mlljj[channel]->Fill(finalstates.M(),wgt);
  h_mll[channel]->Fill(dileptons.M(), wgt);



  return kTRUE;
}

void MySelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void MySelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  if( option == "Data" ){ WriteToFile(Form("%d", 2016), option); }
  else{WriteToFile(Form("%d", DSID), option);}
  cout << "Total number of processed events: " << events_tot << endl;
  cout << "Number of events in ee channel: " << events_ee << endl;
  cout << "Number of events in uu channel: " << events_uu << endl;

}

void MySelector::WriteToFile(TString fileid, TString data_type)
{

  TString filename;

  filename = "hist."+data_type+"."+fileid+".root";

  TFile file("Histograms/"+data_type+"/"+filename, "RECREATE");

  for( const auto & chn:channels ){
    //h_pt1[chn]->Write();
    //h_pt2[chn]->Write();
    //h_eta1[chn]->Write();
    //h_eta2[chn]->Write();
    //h_phi1[chn]->Write();
    //h_phi2[chn]->Write();
    h_pt[chn]->Write();
    h_mlljj[chn]->Write();
    h_mll[chn]->Write();
  }

  file.Close();

  for( const auto & chn:channels ){
    //h_pt1[chn]->Reset();
    //h_pt2[chn]->Reset();
    //h_eta1[chn]->Reset();
    //h_eta2[chn]->Reset();
    //h_phi1[chn]->Reset();
    //h_phi2[chn]->Reset();
    h_pt[chn]->Reset();
    h_mlljj[chn]->Reset();
    h_mll[chn]->Reset();
  }
}