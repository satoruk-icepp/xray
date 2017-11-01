#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "function.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44
#define Nfitparam 5


Double_t SiPMZPosDesign[nMPPC];
Double_t SiPMZResult[nMPPC][Nfitparam];
Double_t SiPMZErr[nMPPC][Nfitparam];
Double_t SiPMZChiSq[nMPPC];
Bool_t SiPMZMeasured[nMPPC];

Double_t SiPMPhiPosDesign[nMPPC];
Double_t SiPMPhiResult[nMPPC][Nfitparam];
Double_t SiPMPhiErr[nMPPC][Nfitparam];
Double_t SiPMPhiChiSq[nMPPC];
Bool_t SiPMPhiMeasured[nMPPC];

void UCIRunAnalysis(int run, Bool_t PhiScan,Bool_t visualize=false) {

  /*-----Define rec tree to be read----*/
  TFile* frec;
  if(PhiScan == true ){
    frec = new TFile(Form("$(MEG2SYS)/analyzer/macros/xec/xray/ForUTokyo/corrected/PhiScan/xray%06d_corr.root", run),"READ");
  }else{	
    frec = new TFile(Form("$(MEG2SYS)/analyzer/macros/xec/xray/ForUTokyo/corrected/ZScan/xray%06d_corr.root", run),"READ");
  }

  std::cout<<"Read rec "<<run<<std::endl;

  TF1* FitFunc[32];
  TGraphErrors* grScaler[32];
  for(int i=0;i<32;i++){
    grScaler[i]=(TGraphErrors*)frec->Get(Form("mppc%d_corr",i));
    std::cout<<"debug"<<std::endl;
    if(grScaler[i]==0){
      break;
    }
    TString graphname=grScaler[i]->GetTitle();
    stringstream ss;
    Int_t ChNum;
    Double_t GraphZPos;
    Double_t GraphPhiPos;
    Double_t Beam;
    ss << graphname;
    ss.ignore(4);
    ss >> ChNum;
    ss.ignore(3);
    ss >> GraphZPos;
    ss.ignore(8);
    ss >> GraphPhiPos;
    //ss.ignore(12);
    //ss >> Beam;
    std::cout<<"debug"<<std::endl;
    SiPMPhiPosDesign[ChNum]=GraphPhiPos;
    SiPMZPosDesign[ChNum]=GraphZPos;
    TString fitfuncname;
    fitfuncname.Form("fit%d",i);
    FitFunc[i] = new TF1(fitfuncname,TwoGaus,-300,300,5);
    Double_t FitResult[5];
    Double_t FitErr[5];
    Double_t ChiSquare;
    if(PhiScan==true){
      FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
      //FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,60);
      FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,60);
      FitFunc[i]->SetParLimits(1,0,5);
      FitFunc[i]->SetParLimits(2,0,5);
      FitFunc[i]->SetParLimits(3,0,1000);
      FitFunc[i]->SetParLimits(4,0,200);//isosceles
      grScaler[i]->Fit(Form("fit%d",i),"MNQ");
      std::cout<<"debug"<<std::endl;
    }else{
      std::cout<<"debug"<<std::endl;
      FitFunc[i]->SetRange(GraphZPos-80,GraphZPos+80);
      //FitFunc[i]->SetParameters(GraphZPos,10,2,500,40);
      FitFunc[i]->SetParameters(GraphZPos,10,2,80,40);
      FitFunc[i]->SetParLimits(0,GraphZPos-80,GraphZPos+80);
      FitFunc[i]->SetParLimits(1,0,15);
      FitFunc[i]->SetParLimits(2,0,25);
      FitFunc[i]->SetParLimits(3,0,5000);
      FitFunc[i]->SetParLimits(4,0,200);//TwoGaus
      grScaler[i]->Fit(Form("fit%d",i),"MNQ");
    }

      std::cout<<"debug"<<std::endl;
    ChiSquare = FitFunc[i]->GetChisquare();
    for(int j=0;j<5;j++){
      FitResult[j] = FitFunc[i]->GetParameter(j);
      FitErr[j] = FitFunc[i]->GetParError(j);
    }

    if(PhiScan==true){
      SiPMPhiChiSq[ChNum]=ChiSquare;
      SiPMPhiMeasured[ChNum]=true;
      for(int k=0;k<Nfitparam;k++){
        SiPMPhiResult[ChNum][k]=FitResult[k];
        SiPMPhiErr[ChNum][k]=FitErr[k];
      }
    }else{
      SiPMZChiSq[ChNum]=ChiSquare;
      SiPMZMeasured[ChNum]=true;
      for(int k=0;k<Nfitparam;k++){
        SiPMZResult[ChNum][k]=FitResult[k];
        SiPMZErr[ChNum][k]=FitErr[k];
      }
    }	
    if(visualize==true){
      std::cout<<"mppc"<<i<<"Channel: "<<ChNum<<"row: "<<ChNum/44<<"line: "<<ChNum%44<<" Position: "<<FitResult[0]<<std::endl;	
    }
  }
  if(visualize==true){
    Int_t visch;
    std::cout<<"channel to visualize: ";
    std::cin>>visch;
    grScaler[visch]->Draw("ap");
    FitFunc[visch]->SetLineColor(kRed);
    FitFunc[visch]->Draw("same");
  }

  return;
}

