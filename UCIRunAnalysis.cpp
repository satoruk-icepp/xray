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

Double_t SiPMZPos[nMPPC];
Double_t SiPMZPosDesign[nMPPC];
Double_t SiPMZErr[nMPPC][Nfitparam];
Double_t SiPMZChiSq[nMPPC];
Bool_t SiPMZMeasured[nMPPC];

Double_t SiPMPhiPos[nMPPC];
Double_t SiPMPhiPosDesign[nMPPC];
Double_t SiPMPhiErr[nMPPC][Nfitparam];
Double_t SiPMPhiChiSq[nMPPC];
Bool_t SiPMPhiMeasured[nMPPC];

void UCIRunAnalysis(int run, Bool_t PhiScan,Bool_t visualize=false) {

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
	ss.ignore(12);
	ss >> Beam;

	SiPMPhiPosDesign[ChNum]=GraphPhiPos;
	SiPMZPosDesign[ChNum]=GraphZPos;
	TString fitfuncname;
	fitfuncname.Form("fit%d",i);
	FitFunc[i] = new TF1(fitfuncname,ArbFunc,-300,300,5);
	Double_t FitErr[5];
	Double_t MeasPos;
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
	Baseline = FitFunc[i]->GetParameter(4);
	ChiSquare = FitFunc[i]->GetChisquare();
	for(int j=0;j<5;j++){
	  FitErr[j] = FitFunc[i]->GetParError(j);
	}

	if(PhiScan==true){
	  SiPMPhiPos[ChNum]=MeasPos;
	  SiPMPhiChiSq[ChNum]=ChiSquare;
	  SiPMPhiMeasured[ChNum]=true;
	  for(int i=0;i<Nfitparam;i++){
		SiPMPhiErr[ChNum][i]=FitErr[i];
	  }

	}else{
	  SiPMZPos[ChNum]=MeasPos;	 
	  SiPMZChiSq[ChNum]=ChiSquare;
	  SiPMZMeasured[ChNum]=true;
	  for(int i=0;i<Nfitparam;i++){
		SiPMZErr[ChNum][i]=FitErr[i];
	  }
	}	
	if(visualize==true){
	  std::cout<<"Channel: "<<ChNum<<" Position: "<<MeasPos<<std::endl;	
	}
  }

  return;
}




// void UCIRunAnalysis(int run, Bool_t PhiScan,int one) {
//   TCanvas* canvas1 =new TCanvas("canvas1","fitting",900,600);
//   TCanvas* canvas2 =new TCanvas("canvas2","one channel",600,600);
//   /*-----Define rec tree to be read----*/
//   TFile* frec;
//   if(PhiScan == true ){
// 	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/PhiScan/xray%06d_graphs.root", run),"READ");
//   }else{	
// 	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/ZScan/xray%06d_graphs.root", run),"READ");
//   }

//   std::cout<<"Read rec "<<run<<std::endl;

//   canvas1->Divide(6,6);
//   TF1* FitFunc[32];
//   TGraphErrors* grScaler[32];
//   for(int i=0;i<32;i++){
// 	grScaler[i]=(TGraphErrors*)frec->Get(Form("mppc%d",i));
// 	TString graphname=grScaler[i]->GetTitle();
// 	stringstream ss;
// 	std::cout<<graphname<<std::endl;
// 	Int_t ChNum;
// 	Double_t GraphZPos;
// 	Double_t GraphPhiPos;
// 	Double_t Beam;
// 	ss << graphname;
// 	ss.ignore(4);
// 	ss >> ChNum;
// 	ss.ignore(3);
// 	ss >> GraphZPos;
// 	ss.ignore(8);
// 	ss >> GraphPhiPos;
// 	ss.ignore(12);
// 	ss >> Beam;
// 	//std::cout<<ChNum<<std::endl;
// 	//std::cout<<GraphZPos<<std::endl;
// 	//std::cout<<GraphPhiPos<<std::endl;
// 	//std::cout<<Beam<<std::endl;
// 	//	sscanf(graphname,"MPPC %d: Z %f mm, Phi %f deg; Beam: %f",&ChNum,&GraphZPos,&GraphPhiPos,&Beam);
// 	TString fitfuncname;
// 	fitfuncname.Form("fit%d",i);
// 	FitFunc[i] = new TF1(fitfuncname,ArbFunc,-300,300,5);
// 	if(PhiScan==true){
// 	  FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
// 	  FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,40);
// 	  FitFunc[i]->SetParLimits(1,0,2);
// 	  FitFunc[i]->SetParLimits(2,0,2);
// 	  FitFunc[i]->SetParLimits(3,0,200);
// 	  FitFunc[i]->SetParLimits(4,0,200);
// 	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
// 	}else{
// 	  FitFunc[i]->SetRange(GraphZPos-80,GraphZPos+80);
// 	  FitFunc[i]->SetParameters(GraphZPos,10,2,500,40);
// 	  FitFunc[i]->SetParLimits(0,GraphZPos-80,GraphZPos+80);
// 	  FitFunc[i]->SetParLimits(1,0,15);
// 	  FitFunc[i]->SetParLimits(2,0,25);
// 	  FitFunc[i]->SetParLimits(3,0,5000);
// 	  FitFunc[i]->SetParLimits(4,0,200);
// 	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
// 	}
// 	Double_t MeasPos = FitFunc[i]->GetParameter(0);
// 	Double_t MeasPosErr = FitFunc[i]->GetParError(0);
// 	Double_t TopWidth = FitFunc[i]->GetParameter(1);
// 	Double_t sigma = FitFunc[i]->GetParameter(2);
// 	Double_t Baseline = FitFunc[i]->GetParameter(4);
// 	std::cout<< "Position:  "<<MeasPos<<"+-"<<MeasPosErr<<std::endl;
// 	std::cout<< "Top Width:  "<< TopWidth <<"  sigma:  "<<sigma<<"  Baseline:  "<<Baseline<<std::endl;
// 	canvas1->cd(i+1);
// 	grScaler[i]->GetXaxis()->SetLimits(MeasPos-30,MeasPos+30);
// 	grScaler[i]->Draw("ap");
// 	FitFunc[i]->SetLineColor(kRed);
// 	FitFunc[i]->Draw("same");
//   }
//   canvas2->cd();
//   grScaler[one]->Draw("ap");
//   FitFunc[one]->SetLineColor(kRed);
//   FitFunc[one]->Draw("same");

//   return;
// }

