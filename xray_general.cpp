#include "xray_general.hpp"

void xray_general(){
  //TCanvas* canvas1=new TCanvas("canvas1","map",600,600);

  TFile *fall =new TFile("$(MEG2SYS)/analyzer/x-ray/xray_allch_dg.root","RECREATE");
  TTree *tall =new TTree("xrayac","xrayac");
  //TGraphErrors* grPhiPosAllch =new TGraphErrors();
  Int_t ChNum;

  Double_t OneChXPos;
  Double_t OneChYPos;

  Double_t OneChPhiPos;
  Double_t OneChPhiPosDesign;
  Double_t OneChPhiPosErr;

  Double_t OneChZPos;
  Double_t OneChZPosDesign;
  Double_t OneChZPosErr;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("XPos",&OneChXPos);
  tall->Branch("YPos",&OneChYPos);

  tall->Branch("PhiPos",&OneChPhiPos);
  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiPosErr",& OneChPhiPosErr);

  tall->Branch("ZPos",&OneChZPos);
  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZPosErr",& OneChZPosErr);

  for(int i=0;i<PhiRunNum;i++){
    OneRunAnalysis(PhiRunList[i],true);
  }

  for(int i=0;i<ZRunNum;i++){
    OneRunAnalysis(ZRunList[i],false);
  }


  Int_t cnt=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;

	OneChXPos=XPosAllch[iCh];
	OneChYPos=YPosAllch[iCh];
	//Phi Position
    OneChPhiPos=PhiPosAllch[iCh];
    OneChPhiPosDesign=PhiPosDesignAllch[iCh];
	//Z Position
	OneChZPos=ZPosAllch[iCh];
    OneChZPosDesign=ZPosDesignAllch[iCh];
	for(i=0;i<5;i++){
    OneChPhiPosErr=PhiFitErrAllch[iCh][i];
    OneChZPosErr=ZFitErrAllch[iCh][i];
	}
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;
}



