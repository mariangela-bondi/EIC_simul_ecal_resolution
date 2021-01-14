//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 12 15:42:30 2021 by ROOT version 6.18/04
// from TTree events/Flattened root tree with event data
// found on file: out_500.root
//////////////////////////////////////////////////////////

#ifndef Cluster_Ecal_Selector_h
#define Cluster_Ecal_Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>

// Headers needed by this particular selector
#include <vector>

#include <string>
using namespace std;


class Cluster_Ecal_Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<ULong64_t> event_id = {fReader, "event_id"};
   TTreeReaderValue<Double_t> evt_true_q2 = {fReader, "evt_true_q2"};
   TTreeReaderValue<Double_t> evt_true_x = {fReader, "evt_true_x"};
   TTreeReaderValue<Double_t> evt_true_y = {fReader, "evt_true_y"};
   TTreeReaderValue<Double_t> evt_true_w2 = {fReader, "evt_true_w2"};
   TTreeReaderValue<Double_t> evt_true_nu = {fReader, "evt_true_nu"};
   TTreeReaderValue<Double_t> evt_true_t_hat = {fReader, "evt_true_t_hat"};
   TTreeReaderValue<Char_t> evt_has_dis_info = {fReader, "evt_has_dis_info"};
   TTreeReaderValue<Double_t> evt_weight = {fReader, "evt_weight"};
   TTreeReaderValue<ULong64_t> hit_count = {fReader, "hit_count"};
   TTreeReaderArray<unsigned long> hit_id = {fReader, "hit_id"};
   TTreeReaderArray<unsigned long> hit_trk_id = {fReader, "hit_trk_id"};
   TTreeReaderArray<unsigned long> hit_ptr_id = {fReader, "hit_ptr_id"};
   TTreeReaderArray<unsigned long> hit_parent_trk_id = {fReader, "hit_parent_trk_id"};
   TTreeReaderArray<string> hit_vol_name = {fReader, "hit_vol_name"};
   TTreeReaderArray<double> hit_x = {fReader, "hit_x"};
   TTreeReaderArray<double> hit_y = {fReader, "hit_y"};
   TTreeReaderArray<double> hit_z = {fReader, "hit_z"};
   TTreeReaderArray<unsigned long> hit_i_rep = {fReader, "hit_i_rep"};
   TTreeReaderArray<unsigned long> hit_j_rep = {fReader, "hit_j_rep"};
   TTreeReaderArray<double> hit_e_loss = {fReader, "hit_e_loss"};
   TTreeReaderValue<ULong64_t> trk_count = {fReader, "trk_count"};
   TTreeReaderArray<unsigned long> trk_id = {fReader, "trk_id"};
   TTreeReaderArray<long> trk_pdg = {fReader, "trk_pdg"};
   TTreeReaderArray<unsigned long> trk_parent_id = {fReader, "trk_parent_id"};
   TTreeReaderArray<long> trk_create_proc = {fReader, "trk_create_proc"};
   TTreeReaderArray<unsigned long> trk_level = {fReader, "trk_level"};
   TTreeReaderArray<double> trk_vtx_x = {fReader, "trk_vtx_x"};
   TTreeReaderArray<double> trk_vtx_y = {fReader, "trk_vtx_y"};
   TTreeReaderArray<double> trk_vtx_z = {fReader, "trk_vtx_z"};
   TTreeReaderArray<double> trk_vtx_dir_x = {fReader, "trk_vtx_dir_x"};
   TTreeReaderArray<double> trk_vtx_dir_y = {fReader, "trk_vtx_dir_y"};
   TTreeReaderArray<double> trk_vtx_dir_z = {fReader, "trk_vtx_dir_z"};
   TTreeReaderArray<double> trk_mom = {fReader, "trk_mom"};
   TTreeReaderValue<ULong64_t> gen_prt_count = {fReader, "gen_prt_count"};
   TTreeReaderArray<unsigned long> gen_prt_id = {fReader, "gen_prt_id"};
   TTreeReaderArray<unsigned long> gen_prt_vtx_id = {fReader, "gen_prt_vtx_id"};
   TTreeReaderArray<unsigned long> gen_prt_pdg = {fReader, "gen_prt_pdg"};
   TTreeReaderArray<unsigned long> gen_prt_trk_id = {fReader, "gen_prt_trk_id"};
   TTreeReaderArray<double> gen_prt_charge = {fReader, "gen_prt_charge"};
   TTreeReaderArray<double> gen_prt_dir_x = {fReader, "gen_prt_dir_x"};
   TTreeReaderArray<double> gen_prt_dir_y = {fReader, "gen_prt_dir_y"};
   TTreeReaderArray<double> gen_prt_dir_z = {fReader, "gen_prt_dir_z"};
   TTreeReaderArray<double> gen_prt_tot_mom = {fReader, "gen_prt_tot_mom"};
   TTreeReaderArray<double> gen_prt_tot_e = {fReader, "gen_prt_tot_e"};
   TTreeReaderArray<double> gen_prt_time = {fReader, "gen_prt_time"};
   TTreeReaderArray<double> gen_prt_polariz_x = {fReader, "gen_prt_polariz_x"};
   TTreeReaderArray<double> gen_prt_polariz_y = {fReader, "gen_prt_polariz_y"};
   TTreeReaderArray<double> gen_prt_polariz_z = {fReader, "gen_prt_polariz_z"};
   TTreeReaderValue<ULong64_t> gen_vtx_count = {fReader, "gen_vtx_count"};
   TTreeReaderArray<unsigned long> gen_vtx_id = {fReader, "gen_vtx_id"};
   TTreeReaderArray<unsigned long> gen_vtx_part_count = {fReader, "gen_vtx_part_count"};
   TTreeReaderArray<double> gen_vtx_x = {fReader, "gen_vtx_x"};
   TTreeReaderArray<double> gen_vtx_y = {fReader, "gen_vtx_y"};
   TTreeReaderArray<double> gen_vtx_z = {fReader, "gen_vtx_z"};
   TTreeReaderArray<double> gen_vtx_time = {fReader, "gen_vtx_time"};
   TTreeReaderArray<double> gen_vtx_weight = {fReader, "gen_vtx_weight"};
   TTreeReaderArray<string> ce_emcal_name = {fReader, "ce_emcal_name"};
   TTreeReaderArray<double> ce_emcal_Etot_dep = {fReader, "ce_emcal_Etot_dep"};
   TTreeReaderArray<int> ce_emcal_Npe = {fReader, "ce_emcal_Npe"};
   TTreeReaderArray<double> ce_emcal_ADC = {fReader, "ce_emcal_ADC"};
   TTreeReaderArray<double> ce_emcal_TDC = {fReader, "ce_emcal_TDC"};
   TTreeReaderArray<double> ce_emcal_xcrs = {fReader, "ce_emcal_xcrs"};
   TTreeReaderArray<double> ce_emcal_ycrs = {fReader, "ce_emcal_ycrs"};
   TTreeReaderArray<double> ce_emcal_zcrs = {fReader, "ce_emcal_zcrs"};


   Cluster_Ecal_Selector(TTree * /*tree*/ =0) {
       hEnergy =0;
       hNumber_crystals=0;
   }
   virtual ~Cluster_Ecal_Selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

    float Ene;
    TH1D *hEnergy;
    TH1D *hNumber_crystals;
    
    
    
    
    
    struct Hit {
        double  x_crs;
        double  y_crs;
        double  z_crs;
        double  Et_dep;
        double  E_digi;
        double  time;
        int  npe;

    };

    struct Cluster{
        
        double C_seed_energy;
        int C_seed_npe;
        double C_seed_x;
        double C_seed_y;
        double C_seed_z;
        double C_energy;
        double C_x;
        double C_y;
        double C_radius;
        double C_theta;
        double C_phi;
        double C_size;
        double C_Energy_tot_simul;
        double C_size_simul;
        
    };
      
        Cluster ComputeCluster(vector<Hit> hit){
            int Size = hit.size();
            double CRYS_ZPOS = 2110; //ECAL zpos from the center of EIC detector
            double Ethr =10.; //MeV
            double Rmoliere = 22.01; // in mm    ->Rmolier for PbWO is 20 mm
            int NR =3;
            double ClusSeed_Ene=0;
            int ClusSeed_npe=0;
            double ClusSeed_xcrs=0;
            double ClusSeed_ycrs=0;
            double ClusSeed_zcrs=0;
            double Clus_Etot = 0;
            double Clus_xx=0;
            double Clus_yy=0;
            double Clus_x=0;
            double Clus_y=0;
            int Clus_size=0;
            double Clus_sigmaX=0;
            double Clus_sigmaY=0;
            double Clus_Radius=0;
            double Clus_Theta=0; //in deg;
            double Clus_phi=0; //in deg;
            int Clus_size_simul=0;
            double Clus_Energy_tot_simul=0;
            TRandom *random_c = new TRandom();
            
            
            double npe_to_ene = 8; //gamma/MeV
            
            //  Cluster Seed;
            for(int i=0; i<Size; i++){
           //     if(hit.at(i).E_digi>Ethr&& hit.at(i).E_digi>ClusSeed_Ene) {
                if( hit.at(i).npe>0 && hit.at(i).npe>ClusSeed_npe){
                    ClusSeed_Ene = hit.at(i).npe /npe_to_ene;
                    //      ClusSeed_Ene = hit.at(i).E_digi;
                    ClusSeed_xcrs =hit.at(i).x_crs;
                    ClusSeed_ycrs =hit.at(i).y_crs;
                    ClusSeed_zcrs =hit.at(i).z_crs;
                    ClusSeed_npe = hit.at(i).npe;
              //      cout << "SEED "<<ClusSeed_xcrs<< " "<<ClusSeed_ycrs<<endl;
                    
            }
        }

            
            //ENERGY TOT simul starting fro Et_dep
            for(int i=0; i<Size; i++){
             if(hit.at(i).Et_dep>Ethr){
                 double Dx = fabs(hit.at(i).x_crs - ClusSeed_xcrs);
                 double Dy = fabs(hit.at(i).y_crs - ClusSeed_ycrs);
               //   if(sqrt(Dx*Dx+Dy*Dy)<=1*Rmoliere){
                 if(Dx<=NR*Rmoliere && Dy<=NR*Rmoliere){
                 Clus_Energy_tot_simul +=hit.at(i).Et_dep;
                 Clus_size_simul++;
                  }
             }
            }
            
            
            //Cluster Energy tot
              for(int i=0; i<Size; i++){
           //       npe_to_ene = 8*random_c->Gaus(1,0.03);
           //       cout<<"******"<< npe_to_ene<<endl;
                  if(hit.at(i).E_digi>Ethr){
                 // if(hit.at(i).npe/npe_to_ene>10){
                      double Dx = fabs(hit.at(i).x_crs - ClusSeed_xcrs);
                      double Dy = fabs(hit.at(i).y_crs - ClusSeed_ycrs);
                  
                  if(Dx<=NR*Rmoliere && Dy<=NR*Rmoliere){
                      
                   //   cout <<hit.at(i).E_digi<< " "<<hit.at(i).x_crs<< " "<<hit.at(i).y_crs<<" "<< Dx << " "<< Dy<< " "<<sqrt(Dx*Dx+Dy*Dy)<<endl;
                   //   Clus_Etot += hit.at(i).npe/npe_to_ene;
                      Clus_Etot += hit.at(i).E_digi;
                      Clus_size++;
                      
                  }
                  } //end if ethr
              }   // for energy tot
         
            // Cluster Center
            
            double w_tot = 0;
            double x,y;
            x=0;
            y=0;
            
            for(int i=0; i<Size; i++){
                  double w1 = std::max(0., (3.45 + std::log(hit.at(i).E_digi / Clus_Etot)));
                  x += w1 * hit.at(i).x_crs;
                  y += w1 * hit.at(i).y_crs;
                  Clus_xx += w1*hit.at(i).x_crs*hit.at(i).x_crs;
                  Clus_yy += w1*hit.at(i).y_crs*hit.at(i).y_crs;
                  w_tot += w1;
              }
            Clus_x = x/w_tot;
            Clus_y = y/w_tot;
            Clus_xx /=w_tot;
            Clus_yy /=w_tot;
            
             // Cluster sigma
            
            double sigmax2 = Clus_xx - std::pow(Clus_x, 2.);
            if (sigmax2 < 0) sigmax2 = 0;
            Clus_sigmaX = std::sqrt(sigmax2);

            double sigmay2 = Clus_yy - std::pow(Clus_y, 2.);
            if (sigmay2 < 0) sigmay2 = 0;
            Clus_sigmaY = std::sqrt(sigmay2);
            
            //Cluster radius
            double radius2 = (sigmax2 + sigmay2);
            if (radius2 < 0) radius2 = 0;
            Clus_Radius = std::sqrt(radius2);
            
            //Cluster theta
          Clus_Theta = (std::atan((std::sqrt(std::pow(Clus_x, 2.) + std::pow(Clus_y, 2.))) /(CRYS_ZPOS+ClusSeed_zcrs))) * (180. / M_PI);
            
            //Cluster phi
            Clus_phi = std::atan2(Clus_x, Clus_y) * (180. / M_PI); //

            Cluster cluster;
            
            cluster.C_seed_energy = ClusSeed_Ene;
            cluster.C_energy = Clus_Etot;
            cluster.C_seed_x = ClusSeed_xcrs;
            cluster.C_seed_y = ClusSeed_ycrs;
            cluster.C_seed_z = ClusSeed_zcrs;
            cluster.C_x = Clus_x;
            cluster.C_y = Clus_y;
            cluster.C_radius = Clus_Radius;
            cluster.C_theta = Clus_Theta;
            cluster.C_phi = Clus_phi;
            cluster.C_size = Clus_size;
            cluster.C_Energy_tot_simul = Clus_Energy_tot_simul;
            cluster.C_size_simul = Clus_size_simul;
            cluster.C_seed_npe = ClusSeed_npe;
            
            
            return cluster;
            }
     
            
    
    
   ClassDef(Cluster_Ecal_Selector,0);

};

#endif

#ifdef Cluster_Ecal_Selector_cxx
void Cluster_Ecal_Selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Cluster_Ecal_Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Cluster_Ecal_Selector_cxx

