#include "xray_general.hpp"
/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

Double_t SiPMAllXPos[nMPPC];
Double_t SiPMAllYPos[nMPPC];
Double_t SiPMAllPhiPos[nMPPC];
Double_t SiPMAllPhiPosErr[nMPPC];
Double_t SiPMAllPhiPosDesign[nMPPC];
Double_t SiPMAllPhiPosGap[nMPPC];
Double_t SiPMAllZPos[nMPPC];
Double_t SiPMAllZPosErr[nMPPC];
Double_t SiPMAllZPosDesign[nMPPC];
Double_t SiPMAllZPosGap[nMPPC];
Bool_t SiPMAllDataQualPhi[nMPPC];
Bool_t SiPMAllDataQualZ[nMPPC];

Int_t PhiRunList[]={305000,305001,305002,305003,305004,305005,305006,//phi_configA
		    305031,305032,305033,305034,305035,305036,305037,//phi_configB
		    305051,305062,305063,305064,305069,305070,305071,//phi_configC
					305085,/*305133,*/305086,305087,305098,305108,//phi_configD
		    305080,305079,305248,305077,305074,//phi_configE
		    305120,305121,305125,305127,305128,//phi_configF
		    305145,305144,305143,305142,305141,//phi_configG
		    305132,305134,305135,305136,305137,//phi_configH
					305235,305234,305233,305232,305231,305243,305242,/*305244,*///phi_configI
		    305208,305209,305217,305223,305224,//phi_configJ
		    305159,305158,305150,305149,305148,//phi_configK
		    305179,305187,305191,305192,305193,//phi_configL
};

Int_t ZRunList[]={305311,305312,305313,305314,305315,305316,305317,305318,//z_configA
				  305322,305323,305324,305330,305331,305332,305333,305334,305341,//z_configB
				  305350,305351,305352,305353,305354,305355,305357,305358,//z_configC
				  305299,305300,305301,305302,305303,305304,305305,305306,//z_configD
				  305363,305365,305366,305367,305368,305370,305371,305372,//z_configE
				  305381,305382,305383,305384,305385,305389,305390,305391,//z_configF
				  305395,305396,305397,305398,305399,305400,305401,305402,//z_configG
				  305413,305414,305415,305416,305417,305418,305419,305420,//z_configH
				  305439,305440,305441,305442,305448,305451,305452,305453,//z_configI
				  305456,305457,305458,305459,305460,305461,305462,305466,//z_configJ
				  305472,305473,305474,305475,305476,305477,305478,//z_configK
				  305512,305513,305514,305515,305516,305517,305518//z_configL
};

Int_t PhiRunNum=sizeof(PhiRunList)/sizeof(PhiRunList[0]);
Int_t ZRunNum=sizeof(ZRunList)/sizeof(ZRunList[0]);

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



