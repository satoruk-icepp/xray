#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t TwoGaus(Double_t *x,Double_t *par);


Double_t SiPMZPos[nMPPC];
Double_t SiPMPhiPos[nMPPC];
Double_t SiPMZPosErr[nMPPC];
Double_t SiPMPhiPosErr[nMPPC];
Double_t BaseLine[nMPPC];

/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

void UCIRunAnalysis(int run, Bool_t PhiScan,int one) {
  TCanvas* canvas1 =new TCanvas("canvas1","fitting",900,600);
  TCanvas* canvas2 =new TCanvas("canvas2","one channel",600,600);
  /*-----Define rec tree to be read----*/
  TFile* frec;
  if(PhiScan == true ){
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/PhiScan/xray%06d_graphs.root", run),"READ");
  }else{	
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/UCIdata/ZScan/xray%06d_graphs.root", run),"READ");
  }

  std::cout<<"Read rec "<<run<<std::endl;

  canvas1->Divide(6,6);
  TF1* FitFunc[32];
  TGraphErrors* grScaler[32];
  for(int i=0;i<32;i++){
	grScaler[i]=(TGraphErrors*)frec->Get(Form("mppc%d",i));
	TString graphname=grScaler[i]->GetTitle();
	stringstream ss;
	//char graphname[64]=grScaler[i]->GetTitle();
	std::cout<<graphname<<std::endl;
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
	//std::cout<<MPPCch<<std::endl;
	//std::cout<<GraphZPos<<std::endl;
	//std::cout<<GraphPhiPos<<std::endl;
	//std::cout<<Beam<<std::endl;
	//	sscanf(graphname,"MPPC %d: Z %f mm, Phi %f deg; Beam: %f",&MPPCch,&GraphZPos,&GraphPhiPos,&Beam);
	TString fitfuncname;
	fitfuncname.Form("fit%d",i);
	FitFunc[i] = new TF1(fitfuncname,TwoGaus,-300,300,5);
	if(PhiScan==true){
	  FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
	  FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,40);
	  FitFunc[i]->SetParLimits(1,0,2);
	  FitFunc[i]->SetParLimits(2,0,2);
	  FitFunc[i]->SetParLimits(3,0,200);
	  FitFunc[i]->SetParLimits(4,0,200);
	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
	}else{
	  FitFunc[i]->SetRange(GraphZPos-80,GraphZPos+80);
	  FitFunc[i]->SetParameters(GraphZPos,10,2,500,40);
	  FitFunc[i]->SetParLimits(0,GraphZPos-80,GraphZPos+80);
	  FitFunc[i]->SetParLimits(1,0,15);
	  FitFunc[i]->SetParLimits(2,0,25);
	  FitFunc[i]->SetParLimits(3,0,5000);
	  FitFunc[i]->SetParLimits(4,0,200);
	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
	}
	Double_t MeasPos = FitFunc[i]->GetParameter(0);
	Double_t MeasPosErr = FitFunc[i]->GetParError(0);
	Double_t TopWidth = FitFunc[i]->GetParameter(1);
	Double_t sigma = FitFunc[i]->GetParameter(2);
	Double_t Baseline = FitFunc[i]->GetParameter(4);
	std::cout<< "Position:  "<<MeasPos<<"+-"<<MeasPosErr<<std::endl;
	std::cout<< "Top Width:  "<< TopWidth <<"  sigma:  "<<sigma<<"  Baseline:  "<<Baseline<<std::endl;
	canvas1->cd(i+1);
	//	grScaler[i]->GetXaxis()->SetLimits(MeasPos-30,MeasPos+30);
	grScaler[i]->Draw("ap");
	FitFunc[i]->SetLineColor(kRed);
	FitFunc[i]->Draw("same");
  }
  canvas2->cd();
  grScaler[one]->Draw("ap");
  FitFunc[one]->SetLineColor(kRed);
  FitFunc[one]->Draw("same");


  /*
	canvas3->cd();
	grPosition->SetTitle("Z Position Measurement;SiPM Channel;Z Position[mm]");
	grPosition->SetMarkerStyle(22);
	grPosition->  SetMarkerColor(2);
	grPosition->Draw("AP");
	grPositionDesign->SetMarkerStyle(22);
	grPositionDesign->Draw("p");
	canvas4->cd();
	grScalerABLS[OneCh]->SetTitle("X-ray Signal;Phi Position[deg];Trigger Rate(Back ground Subtracted)[Hz]");
	grScalerABLS[OneCh]->Draw("apl");
	FitFunc[OneCh]->SetLineColor(2);
	FitFunc[OneCh]->Draw("same");
	grScalerResidual[OneCh]->Draw("same");
	fout->cd();
	tout->Write();
	fout->Close();*/
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

Double_t TwoGaus(Double_t *x,Double_t *par){
  Double_t xx=x[0];
  Double_t mean=par[0];
  Double_t topwidth=par[1];
  Double_t meanfg=par[0]-par[1]/2;
  Double_t meansg=par[0]+par[1]/2;
  Double_t sigma=par[2];
  Double_t height=par[3];
  Double_t offset=par[4];
  Double_t f;
  if(xx<meanfg){
	f=offset+height*TMath::Gaus(xx,meanfg,sigma,true);
  }else if(xx>=meanfg&&xx<meansg){
	f=offset+ height*TMath::Gaus(meanfg,meanfg,sigma,true);
  }else{
	f=offset+ height*TMath::Gaus(xx,meansg,sigma,true);
  }
  return f;
}

  Double_t ArbFunc(Double_t *x,Double_t *par){
	Float_t xx=x[0];
	//Double_t f = par[0]/sqrt(2*TMath::Pi()*par[1])*TMath::Exp(-pow(xx-par[2],2)/(2*par[1]));
	//par[0]:mean
	//par[1]:topwidth
	//par[2]:zerowidth-topwidth
	//par[3]:height
	Double_t topstart=par[0]-par[1]/2;
	Double_t topend=par[0]+par[1]/2;
	Double_t redge=topstart-par[2]/2;
	Double_t fedge=topend+par[2]/2;
	Double_t height=par[3];
	Double_t f;
	if(xx<redge||xx>=fedge){
	  f=0;
	}else if(xx>=redge&&xx<topstart){
	  Double_t ratio=std::abs((xx-redge)/(topstart-redge));
	  f= height*ratio;
	}else if(xx>=topstart&&xx<topend){
	  f=height;  
	}else{
	  Double_t ratio=std::abs((xx-topend)/(fedge-topend));
	  f=height*(1-ratio);
	}
	return f;
  }
                                  
