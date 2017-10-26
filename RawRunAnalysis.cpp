// #include "TFile.h"
// #include "TTree.h"
// #include "TBranch.h"
// #include "TClonesArray.h"
//#include "Riostream.h"
//#ifndef __CINT__
//#include "ROMETreeInfo.h"
#include "function.h"
#include "RawRunList.h"
#include "transcoordinate.h"
//#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44
#define Nfitparam 5

Double_t XPosAllch[nMPPC];
Double_t YPosAllch[nMPPC];

Double_t PhiPosDesignAllch[nMPPC];
Double_t PhiFitResultAllch[nMPPC][Nfitparam];
Double_t PhiFitErrAllch[nMPPC][Nfitparam];
Double_t PhiChiSqAllch[nMPPC];
Bool_t PhiMeasuredAllch[nMPPC];

Double_t ZPosDesignAllch[nMPPC];
Double_t ZFitResultAllch[nMPPC][Nfitparam];
Double_t ZFitErrAllch[nMPPC][Nfitparam];
Double_t ZChiSqAllch[nMPPC];
Bool_t ZMeasuredAllch[nMPPC];

Int_t PhiRunNum=sizeof(PhiRunList)/sizeof(PhiRunList[0]);
Int_t ZRunNum=sizeof(ZRunList)/sizeof(ZRunList[0]);

void RawRunAnalysis(int run,Bool_t PhiScan, Bool_t visualize =false) {
  Int_t ChNum;

  Double_t PosDesign;

  /*-----Define rec tree to be read----*/
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/recfiles/rec%06d.root", run),"READ");
  TTree *rec = (TTree*)frec->Get("rec");
  Int_t nEvent = rec->GetEntries();
  std::cout<<"Read rec "<<run<<". "<<nEvent<<" scaler events."<<std::endl;

  /*-----Read PM RunHeader----*/
  TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");
  MEGXECPMRunHeader *pmrh = 0;

  Int_t valid[kNchforXray] ={};
  Int_t MPPCindex[kNchforXray] ={};
  Float_t MPPCXYZ[kNchforXray][3] ={};
  Float_t MPPCPhi[kNchforXray]={};
  Int_t ScalerCh =0;
  // Double_t theta= -TMath::Pi()*0.0/180;
  for (Int_t iPM = 0; iPM < nMPPC; iPM++) {
    pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
    if (pmrh->GetDRSAddress() >= 0) {
      ScalerCh = pmrh->GetDRSChipID() * 8 + pmrh->GetDRSChannelOnChip();
      //define whether measurement is valid or not
      if (ScalerCh < 0 || ScalerCh > kNchforXray) continue;
      valid[ScalerCh] = 1;
      MPPCindex[ScalerCh] = iPM;
      MPPCXYZ[ScalerCh][0] = pmrh->GetXYZAt(0) * 10;
	  MPPCXYZ[ScalerCh][1] = pmrh->GetXYZAt(1) * 10;
	  MPPCXYZ[ScalerCh][2] = pmrh->GetXYZAt(2) * 10;
      MPPCPhi[ScalerCh] = XYZ2Phi(MPPCXYZ[ScalerCh][0], MPPCXYZ[ScalerCh][1], MPPCXYZ[ScalerCh][2]);
    }
  }
  /*-----Define rec tree to be read----*/
  TClonesArray* recTRGScaler = new TClonesArray("MEGTRGScaler");
  TBranch *branchrectrgscaler = 0;
  branchrectrgscaler = rec->GetBranch("trgscaler.");
  if (branchrectrgscaler) branchrectrgscaler->SetAddress(&recTRGScaler);
  MEGXRAYData* recXRAYData = new MEGXRAYData();
  TBranch *branchrecxraydata = 0;
  branchrecxraydata = rec->GetBranch("xraydata.");
  if (branchrecxraydata) branchrecxraydata->SetAddress(&recXRAYData);

  gStyle->SetTitleW(2);
  TGraph* grBeamPosition= new TGraph();
  TGraph* grScaler[kNchforXray];
  TF1* FitFunc[kNchforXray];

  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    grScaler[iScalerCh] = new TGraph();
    TString fitfuncname;
    fitfuncname.Form("fit%d",iScalerCh);
    FitFunc[iScalerCh] = new TF1(fitfuncname,TwoGaus,-300,300,5);
  }

  Int_t Scaler_buf;
  Float_t BeamZ_buf, BeamPhi_buf;
  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {
    branchrectrgscaler->GetEntry(iEvent);
    branchrecxraydata->GetEntry(iEvent);
    // read the data quality
    if (recXRAYData->GetIsGood()) continue;
    for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      Scaler_buf = ((MEGTRGScaler*)(recTRGScaler->At(iScalerCh)))->GetScaler();
      BeamZ_buf = recXRAYData->GetBeamZ();
      BeamPhi_buf = recXRAYData->GetBeamPhi();
      if (PhiScan == true) {
		if (TMath::Abs(BeamZ_buf - MPPCXYZ[iScalerCh][2]) < 18. ){
		  grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamPhi_buf, Scaler_buf);
		}
      } 
      else if (PhiScan == false) {
		if (TMath::Abs(BeamPhi_buf - MPPCPhi[iScalerCh]) < 1.5 && TMath::Abs(MPPCXYZ[iScalerCh][2])< 157 ){
		  grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamZ_buf, Scaler_buf);
		}
      }
    }
  }

  /*Writing to File...*/
  //  fout->cd();
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    if (valid[iScalerCh]&&grScaler[iScalerCh]->GetN()>10){
	  Double_t FitResult[Nfitparam];
	  Double_t FitErr[Nfitparam];
      ChNum=MPPCindex[iScalerCh];
	  XPosAllch[ChNum]=MPPCXYZ[iScalerCh][0];
	  YPosAllch[ChNum]=MPPCXYZ[iScalerCh][1];
      if(PhiScan==true){
		FitFunc[iScalerCh]->SetRange(MPPCPhi[iScalerCh]-15,MPPCPhi[iScalerCh]+15);
		FitFunc[iScalerCh]->SetParLimits(0,MPPCPhi[iScalerCh]-5,MPPCPhi[iScalerCh]+5);
		FitFunc[iScalerCh]->SetParLimits(1,0,5);
		FitFunc[iScalerCh]->SetParLimits(2,0,5);
		FitFunc[iScalerCh]->SetParLimits(3,0,80);
		FitFunc[iScalerCh]->SetParameters(MPPCPhi[iScalerCh],0.7,0.2,10,60);
		grScaler[iScalerCh]->Fit(Form("fit%d",iScalerCh),"MNQ");
		PosDesign=MPPCPhi[iScalerCh];
		PhiPosDesignAllch[ChNum]=MPPCPhi[iScalerCh];
		PhiChiSqAllch[ChNum]=FitFunc[iScalerCh]->GetChisquare();
		PhiMeasuredAllch[ChNum]=true;


      }else if(PhiScan==false){
		FitFunc[iScalerCh]->SetRange(MPPCXYZ[iScalerCh][2]-80,MPPCXYZ[iScalerCh][2]+80);
		FitFunc[iScalerCh]->SetParLimits(0,MPPCXYZ[iScalerCh][2]-80,MPPCXYZ[iScalerCh][2]+80);	
		FitFunc[iScalerCh]->SetParLimits(1,0,15);
		FitFunc[iScalerCh]->SetParLimits(2,0,25);
		FitFunc[iScalerCh]->SetParLimits(3,0,5000);
		FitFunc[iScalerCh]->SetParLimits(4,0,200);
		FitFunc[iScalerCh]->SetParameters(MPPCXYZ[iScalerCh][2],10,2,80,40);
		grScaler[iScalerCh]->Fit(Form("fit%d",iScalerCh),"MNQ");
		PosDesign=MPPCXYZ[iScalerCh][2];
		ZPosDesignAllch[ChNum]=MPPCXYZ[iScalerCh][2];
		ZChiSqAllch[ChNum]=FitFunc[iScalerCh]->GetChisquare();
		ZMeasuredAllch[ChNum]=true;
      }

	  for(int i=0;i<Nfitparam;i++){
		FitResult[i]=FitFunc[iScalerCh]->GetParameter(i);
		FitErr[i]=FitFunc[iScalerCh]->GetParError(i);
	  }

	  if(PhiScan==true){
		for(int i=0;i<5;i++){
		  PhiFitResultAllch[ChNum][i]=FitResult[i];
		  PhiFitErrAllch[ChNum][i]=FitErr[i];
		}
		if(visualize==true){
		  std::cout<<"iScalerCh: "<<iScalerCh<<" Channel: "<<ChNum<<" Position: "<<PhiFitResultAllch[ChNum][0]<<std::endl;
		}
	  }else{
		for(int i=0;i<5;i++){
		  ZFitResultAllch[ChNum][i]=FitResult[i];
		  ZFitErrAllch[ChNum][i]=FitErr[i];
		}
		if(visualize==true){
		  std::cout<<"iScalerCh: "<<iScalerCh<<" Channel: "<<ChNum<<" Position: "<<ZFitResultAllch[ChNum][0]<<std::endl;
		}
	  }
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
