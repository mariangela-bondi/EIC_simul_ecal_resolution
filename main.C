#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include <TLine.h>
#include <iostream>
#include <string>
#include <vector>
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TProof.h"
#include "TF1.h"
#include "TH1.h"
#include <TGraphErrors.h>
#include "Cluster_Ecal_Selector.h"
using namespace std;

int main(){
    string run[20] = {"500","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","6.5","7","7.5","8","8.5","9","9.5","10"};
    float Energy[20] = {500, 1000, 1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000};
    float Energy_gr[20] = {0.5, 1, 1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10};
    float EResolution[20];
    float EResolution_error[20];

    for(int ii=0; ii<20; ii++){
        Cluster_Ecal_Selector *mySelector = new Cluster_Ecal_Selector();
        cout <<"Beam energy: "<< run[ii]<<endl;
        TFile *f=new TFile(Form("/Users/Mariangela/work/simul_eic/rootfile/ecal_20cm_12x12sipmSat_pbwo_fiber_carbon_1mm_morepoint/out_%s.root",run[ii].c_str()));

        TTree *t=(TTree*)f->Get("events");
        
        mySelector->Ene =Energy[ii];
        t->Process(mySelector);
        mySelector->hEnergy->Draw();
        TF1  *fit_ene = new TF1("fit_ene","gaus",-19000,19000);
        mySelector->hEnergy->Fit("fit_ene","R");
        
        double par_ene[3];
        double ErrPar[3];
        fit_ene->GetParameters(par_ene);
        ErrPar[0]= fit_ene->GetParError(0);
        ErrPar[1]= fit_ene->GetParError(1);
        ErrPar[2]= fit_ene->GetParError(2);
        EResolution[ii] = (par_ene[2]/par_ene[1])*100;
     //   EResolution[ii] = (par_ene[2]/Energy[ii])*100;
        EResolution_error[ii] = (sqrt((1/pow(par_ene[1],2))*pow(ErrPar[2],2))+pow((-par_ene[2]/pow(par_ene[1],2)),2)*pow(ErrPar[1],2))*100;
        EResolution_error[ii] = (sqrt((1/pow(Energy[ii],2))*pow(ErrPar[2],2)))*100;
        cout <<EResolution[ii]<< " "<<EResolution_error[ii]<<endl;
        
        TFile *f_out=new TFile(Form("/Users/Mariangela/work/simul_eic/rootfile/cluster_20cm_12x12sipmSat_pbwo_fiber_carbon_1mm_morepoint/out_%s.root",run[ii].c_str()),"RECREATE");
        f_out->cd();
        mySelector->hEnergy->Write();
        mySelector->hNumber_crystals->Write();
    //    delete fit_npe;
        delete mySelector;
        f_out->Close();

    }
    TGraphErrors* gr_resolution = new TGraphErrors(20,Energy_gr, EResolution, 0,EResolution_error);
    TFile *f_out2=new TFile(Form("/Users/Mariangela/work/simul_eic/rootfile/cluster_20cm_12x12sipmSat_pbwo_fiber_carbon_1mm_morepoint/resolution_3x3.root"),"RECREATE");
    f_out2->cd();
    gr_resolution->Write();
    gr_resolution->Draw("");
    f_out2->Close();
}
