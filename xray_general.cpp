#include "xray_general.hpp"
/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/


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
  Double_t OneChPhiPosGap;
  Bool_t OneChDataQualPhi;

  Double_t OneChZPos;
  Double_t OneChZPosDesign;
  Double_t OneChZPosErr;
  Double_t OneChZPosGap;
  Bool_t OneChDataQualZ;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("XPos",&OneChXPos);
  tall->Branch("YPos",&OneChYPos);

  tall->Branch("PhiPos",&OneChPhiPos);
  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiPosErr",& OneChPhiPosErr);
  tall->Branch("DataQualPhi",& OneChDataQualPhi);

  tall->Branch("ZPos",&OneChZPos);
  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZPosErr",& OneChZPosErr);
  tall->Branch("DataQualZ",& OneChDataQualZ);

  for(int i=0;i<PhiRunNum;i++){
    OneRunAnalysis(PhiRunList[i],0);
  }

  for(int i=0;i<ZRunNum;i++){
    OneRunAnalysis(ZRunList[i],1);
  }


  Int_t cnt=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;

	OneChXPos=SiPMAllXPos[iCh];
	OneChYPos=SiPMAllYPos[iCh];
	//Phi Position
    OneChPhiPos=SiPMAllPhiPos[iCh];
    OneChPhiPosErr=SiPMAllPhiPosErr[iCh];
    OneChPhiPosDesign=SiPMAllPhiPosDesign[iCh];
    OneChDataQualPhi=SiPMAllDataQualPhi[iCh];
	SiPMAllPhiPosGap[iCh]=OneChPhiPos-OneChPhiPosDesign;
	//Z Position
	OneChZPos=SiPMAllZPos[iCh];
    OneChZPosErr=SiPMAllZPosErr[iCh];
    OneChZPosDesign=SiPMAllZPosDesign[iCh];
	OneChZPosGap=OneChZPos-OneChZPosDesign;
    OneChDataQualZ=SiPMAllDataQualZ[iCh];
	SiPMAllZPosGap[iCh]=OneChZPos-OneChZPosDesign;
	/*if(OneChDataQualPhi==true){
	//grPhiPosAllch->SetPoint(cnt,iCh,OneChPhiPos-OneChPhiPosDesign);
	//grPhiPosAllch->SetPointError(cnt,0,OneChPhiPosErr);
	cnt++;
	}*/
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  //grPhiPosAllch->Draw("ap");
  //canvas1->cd();
  //InnerGeometry(SiPMAllPhiPosGap,5,-5);
  return;
}



