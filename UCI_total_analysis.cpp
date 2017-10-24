#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "function.h"
#include "RunList.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44
#define Nfitparam 5

void UCIRunAnalysis(int run,Bool_t PhiScan);
Bool_t DataQual(Double_t *FitErr, Bool_t PhiScan, Double_t Position, Double_t Design);
Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t TwoGaus(Double_t *x,Double_t *par);
Double_t SiPMZPos[nMPPC];
Double_t SiPMZPosDesign[nMPPC];
Bool_t SiPMZDataQual[nMPPC];
Double_t SiPMZErr[Nfitparam][nMPPC];

Double_t SiPMPhiPos[nMPPC];
Double_t SiPMPhiPosDesign[nMPPC];
Bool_t SiPMPhiDataQual[nMPPC];
Double_t SiPMPhiErr[Nfitparam][nMPPC];

/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

Int_t PhiRunNum=sizeof(PhiRunList)/sizeof(PhiRunList[0]);
Int_t ZRunNum=sizeof(ZRunList)/sizeof(ZRunList[0]);

void UCI_total_analysis(){
  TFile *fall =new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_DK.root","RECREATE");
  TTree *tall =new TTree("uci","uci");

  Int_t ChNum;

  Double_t OneChPhiPos;
  Double_t OneChPhiPosDesign;
  Bool_t OneChPhiDataQual;
  Double_t OneChPhiErr[Nfitparam];

  Double_t OneChZPos;
  Double_t OneChZPosDesign;
  Double_t OneChZErr[Nfitparam];
  Bool_t OneChZDataQual;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("PhiPos",&OneChPhiPos);
  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiDataQual",& OneChPhiDataQual);
  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("PhiFitErr%d",i),&OneChPhiErr[i]);
  }

  tall->Branch("ZPos",&OneChZPos);
  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZDataQual",&OneChZDataQual);  
  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("ZFitErr%d",i),&OneChZErr[i]);
  }

  for(int i=0;i<PhiRunNum;i++){
    UCIRunAnalysis(PhiRunList[i],true);
  }

  for(int i=0;i<ZRunNum;i++){
    UCIRunAnalysis(ZRunList[i],false);
  }

  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;
	//Phi Position
    OneChPhiPos=SiPMPhiPos[iCh];

    OneChPhiPosDesign=SiPMPhiPosDesign[iCh];
    OneChPhiDataQual=SiPMPhiDataQual[iCh];
	for(int i=0;i<Nfitparam;i++){
	  OneChPhiErr[i]=SiPMPhiErr[i][iCh];
	}

	//Z Position
	OneChZPos=SiPMZPos[iCh];
	OneChZPosDesign=SiPMZPosDesign[iCh];
    OneChZDataQual=SiPMZDataQual[iCh];
	for(int i=0;i<Nfitparam;i++){
	  OneChZErr[i]=SiPMZErr[i][iCh];
	}
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;
}

void UCIRunAnalysis(int run, Bool_t PhiScan) {

  /*-----Define rec tree to be read----*/
  TFile* frec;
  if(PhiScan == true ){
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/PhiScan/xray%06d_graphs.root", run),"READ");
  }else{	
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/ZScan/xray%06d_graphs.root", run),"READ");
  }

  std::cout<<"Read rec "<<run<<std::endl;

  TF1* FitFunc[32];
  TGraphErrors* grScaler[32];
  for(int i=0;i<32;i++){
	grScaler[i]=(TGraphErrors*)frec->Get(Form("mppc%d",i));
	TString graphname=grScaler[i]->GetTitle();
	stringstream ss;
	Int_t MPPCch;
	Double_t GraphZPos;
	Double_t GraphPhiPos;
	Double_t Beam;
	ss << graphname;
	ss.ignore(4);
	ss >> MPPCch;
	ss.ignore(3);
	ss >> GraphZPos;
	ss.ignore(8);
	ss >> GraphPhiPos;
	ss.ignore(12);
	ss >> Beam;

	SiPMPhiPosDesign[MPPCch]=GraphPhiPos;
	SiPMZPosDesign[MPPCch]=GraphZPos;
	TString fitfuncname;
	fitfuncname.Form("fit%d",i);
	FitFunc[i] = new TF1(fitfuncname,ArbFunc,-300,300,5);
	Double_t FitErr[5];
	Double_t MeasPos;
	Double_t TopWidth;
	Double_t sigma;
	Double_t Height;
	Double_t Baseline;
	if(PhiScan==true){
	  FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
	  //FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,60);
	  FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,60);
	  FitFunc[i]->SetParLimits(1,0,5);
	  FitFunc[i]->SetParLimits(2,0,5);
	  FitFunc[i]->SetParLimits(3,0,1000);
	  FitFunc[i]->SetParLimits(4,0,200);//isosceles
	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
	}else{
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
	
	MeasPos = FitFunc[i]->GetParameter(0);
	TopWidth = FitFunc[i]->GetParameter(1);
	sigma = FitFunc[i]->GetParameter(2);
	Height = FitFunc[i]->GetParameter(3);
	Baseline = FitFunc[i]->GetParameter(4);

	for(int j=0;j<5;j++){
	  FitErr[j] = FitFunc[i]->GetParError(j);
	}

	if(PhiScan==true){
	  SiPMPhiPos[MPPCch]=MeasPos;
	  SiPMPhiDataQual[MPPCch]=DataQual(FitErr,PhiScan,MeasPos,GraphPhiPos);
	  for(int i=0;i<Nfitparam;i++){
		SiPMPhiErr[i][MPPCch]=FitErr[i];
	  }
	}else{
	  SiPMZPos[MPPCch]=MeasPos;	 
	  SiPMZDataQual[MPPCch]=DataQual(FitErr,PhiScan,MeasPos,GraphZPos);
	  for(int i=0;i<Nfitparam;i++){
		SiPMZErr[i][MPPCch]=FitErr[i];
	  }
	}
  }
  return;
}

  //________________________________________________________________________________


  Int_t arraysearch(std::vector<Double_t> array, Double_t value){
	Bool_t found=false;
	Int_t sizearr=array.size();
	int i=-1;
	if(sizearr!=0){
	  for (i = 0; i < sizearr; i++) {
		if (array[i]==value) {
		  found=true;
		  break;
		}
	  }
	}
	if(found==false){
	  i=-1;  
	}
	return i;
  }

Bool_t DataQual(Double_t *FitErr, Bool_t PhiScan, Double_t Position, Double_t Design){
  Bool_t Quality=true;
  Double_t PhiErrMax[5]={0.1,0.5,0.2,8,0.5};
  Double_t ZErrMax[5]={0.5,2,1,100,0.5};
  Double_t Gap=std::abs(Position-Design);
  Double_t PhiGapMax=1;
  Double_t ZGapMax=20;
  if(PhiScan==true){
	for(int i=0;i<5;i++){
	  if(FitErr[i]>PhiErrMax[i]){
		Quality=false;
		break;
	  }
	  if(Gap>PhiGapMax){
		Quality=false;
	  }
	}
  }else{
	for(int i=0;i<5;i++){
	  if(FitErr[i]>ZErrMax[i]){
		Quality=false;
		break;
	  }
	}
	if(Gap>ZGapMax){
	  Quality=false;	  
	}
  }
  return Quality;
}
