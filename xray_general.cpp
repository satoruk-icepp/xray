#include "xray_general.hpp"

void xray_general(){
  //TCanvas* canvas1=new TCanvas("canvas1","map",600,600);

  TFile *fall =new TFile("$(MEG2SYS)/analyzer/x-ray/xray_raw_dg.root","RECREATE");
  TTree *tall =new TTree("xrayac","xrayac");
  //TGraphErrors* grPhiPosAllch =new TGraphErrors();
  Int_t ChNum;

  Double_t XPos;
  Double_t YPos;

  Double_t PhiPos;
  Double_t PhiPosDesign;
  Double_t PhiFitErr[5];
  Double_t PhiChiSq;
  Bool_t PhiMeasured;

  Double_t ZPos;
  Double_t ZPosDesign;
  Double_t ZFitErr[5];
  Double_t ZChiSq;
  Bool_t ZMeasured;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("XPos",&XPos);
  tall->Branch("YPos",&YPos);

  tall->Branch("PhiPos",&PhiPos);
  tall->Branch("PhiPosDesign",&PhiPosDesign);
  tall->Branch("PhiChiSq",&PhiChiSq);
  tall->Branch("PhiMeasured",&PhiMeasured);

  tall->Branch("ZPos",&ZPos);
  tall->Branch("ZPosDesign",&ZPosDesign);
  tall->Branch("ZChiSq",&ZChiSq);
  tall->Branch("ZMeasured",&ZMeasured);

  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("PhiFitErr%d",i),&PhiFitErr[i]);
	tall->Branch(Form("ZFitErr%d",i),&ZFitErr[i]);
  }

  for(int i=0;i<PhiRunNum;i++){
    OneRunAnalysis(PhiRunList[i],true);
  }

  for(int i=0;i<ZRunNum;i++){
    OneRunAnalysis(ZRunList[i],false);
  }


  Int_t cnt=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;

	XPos=XPosAllch[iCh];
	YPos=YPosAllch[iCh];
	//Phi Position
    PhiPos=PhiPosAllch[iCh];
    PhiPosDesign=PhiPosDesignAllch[iCh];
    PhiChiSq=PhiChiSqAllch[iCh];
	PhiMeasured=PhiMeasuredAllch[iCh];
	//Z Position
	ZPos=ZPosAllch[iCh];
    ZPosDesign=ZPosDesignAllch[iCh];
    ZChiSq=ZChiSqAllch[iCh];
	ZMeasured=ZMeasuredAllch[iCh];
	for(int i=0;i<5;i++){
	  PhiFitErr[i]=PhiFitErrAllch[iCh][i];
	  ZFitErr[i]=ZFitErrAllch[iCh][i];
	}
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;
}



