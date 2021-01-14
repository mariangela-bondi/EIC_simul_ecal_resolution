#define Cluster_Ecal_Selector_cxx
// The class definition in Cluster_Ecal_Selector.h has been generated automatically
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
// root> T->Process("Cluster_Ecal_Selector.C")
// root> T->Process("Cluster_Ecal_Selector.C","some options")
// root> T->Process("Cluster_Ecal_Selector.C+")
//


#include "Cluster_Ecal_Selector.h"
#include <TH2.h>
#include <TStyle.h>

void Cluster_Ecal_Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void Cluster_Ecal_Selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();
    float dE = 10;
  //  if(Ene>1999) dE = 50;
    float  min = Ene - 200;
    if(Ene>999) min =Ene-2000;
    float   max = Ene + 200;
    int N = (int) (max-min) / dE;
    cout<< "N "<< N<<"min "<< min<<"max "<<max<<"de "<<dE<<endl;
    hEnergy = new TH1D("hEnergy", "hEnergy", N, min, max);
    hNumber_crystals = new TH1D("hNumber_crystals","hNumber_crystals", 70, -0.5, 69.5 );
}

Bool_t Cluster_Ecal_Selector::Process(Long64_t entry)
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

    int nsize =ce_emcal_Etot_dep.GetSize();
    vector<Hit> hhit;
    for(int j=0; j<nsize; j++){
        
        Hit hit;
        hit.x_crs = ce_emcal_xcrs[j];
        hit.y_crs = ce_emcal_ycrs[j];
        hit.z_crs = ce_emcal_zcrs[j];
        hit.Et_dep = ce_emcal_Etot_dep[j];
        hit.E_digi = ce_emcal_ADC[j];
        hit.time = ce_emcal_TDC[j];
        hit.npe = ce_emcal_Npe[j];
        hhit.push_back(hit);
    
             }
    
  Cluster cluster;
  cluster =  ComputeCluster(hhit);
  //  cout << "calibrtion"<<cluster.C_energy<<" "<<Ene<< " "<<cluster.C_energy/Ene<<endl;
    hEnergy->Fill( cluster.C_energy);
    hNumber_crystals->Fill(cluster.C_size);
    
    
    
   return kTRUE;
}

void Cluster_Ecal_Selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Cluster_Ecal_Selector::Terminate()
{
//    hEnergy->Draw();

    
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
