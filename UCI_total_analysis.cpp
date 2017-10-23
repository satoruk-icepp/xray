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

void UCIRunAnalysis(int run,Bool_t PhiScan);
Bool_t DataQual(Double_t *FitErr, Bool_t PhiScan, Double_t Position, Double_t Design);
Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t TwoGaus(Double_t *x,Double_t *par);
Double_t SiPMZPos[nMPPC];
Double_t SiPMZPosDesign[nMPPC];
Bool_t SiPMZDataQual[nMPPC];
Double_t SiPMZPosErr[nMPPC];
Double_t SiPMZHeightErr[nMPPC];
Double_t SiPMZTopWidthErr[nMPPC];
Double_t SiPMZSigmaErr[nMPPC];
Double_t SiPMZBaseLineErr[nMPPC];

Double_t SiPMPhiPos[nMPPC];
Double_t SiPMPhiPosDesign[nMPPC];
Bool_t SiPMPhiDataQual[nMPPC];
Double_t SiPMPhiPosErr[nMPPC];
Double_t SiPMPhiTopWidthErr[nMPPC];
Double_t SiPMPhiSigmaErr[nMPPC];
Double_t SiPMPhiHeightErr[nMPPC];
Double_t SiPMPhiBaseLineErr[nMPPC];


/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

Int_t PhiRunList[]={
  305000,/*305001,*/305002,305003,305004,305005,305006,//phi_configA
  305031,305032,305033,305034,305035,305036,305037,//phi_configB
  305051,305062,305063,305064,305069,305070,305071,//phi_configC
  305085,/*305133,*/305086,305087,305098,305108,//phi_configD
  305080,305079,/*305248,*/305077,305074,//phi_configE
  305120,305121,305125,305127,305128,//phi_configF
  305145,305144,305143,305142,305141,//phi_configG
  305132,305134,305135,305136,305137,//phi_configH
  305235,305234,305233,305232,305231,/*305243,305242,305244,*///phi_configI
  305208,305209,305217,305223,305224,//phi_configJ
  305159,305158,305150,305149,305148,//phi_configK
  305179,305187,305191,305192,305193,//phi_configL
};


Int_t ZRunList[]={
  305311,305312,305313,305314,305315,305316,305317,305318,//z_configA
  305322,305323,305324,305330,305331,305332,305333,305334,305341,//z_configB
  305350,305351,305352,305353,305354,305355,305357,305358,//z_configC
  305299,305300,305301,305302,305303,305304,305305,305306,//z_configD
  305363,305365,305366,305367,305368,305370,305371,/*305372,*///z_configE
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

void UCI_total_analysis(){
  TFile *fall =new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch.root","RECREATE");
  TTree *tall =new TTree("uci","uci");

  Int_t ChNum;

  Double_t OneChPhiPos;
  Double_t OneChPhiPosDesign;
  Bool_t OneChPhiDataQual;
  Double_t OneChPhiPosErr;
  Double_t OneChPhiTopWidthErr;
  Double_t OneChPhiSigmaErr;
  Double_t OneChPhiHeightErr;
  Double_t OneChPhiBaseLineErr;


  Double_t OneChZPos;
  Double_t OneChZPosDesign;
  Double_t OneChZPosErr;
  Double_t OneChZTopWidthErr;
  Double_t OneChZSigmaErr;
  Double_t OneChZHeightErr;
  Double_t OneChZBaseLineErr;

  //Double_t OneChZPosGap;
  Bool_t OneChZDataQual;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("PhiPos",&OneChPhiPos);
  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiDataQual",& OneChPhiDataQual);
  tall->Branch("PhiPosErr",& OneChPhiPosErr);
  tall->Branch("PhiTopWidthErr",& OneChPhiTopWidthErr);
  tall->Branch("PhiSigmaErr",& OneChPhiSigmaErr);
  tall->Branch("PhiHeightErr",& OneChPhiHeightErr);
  tall->Branch("PhiBaseLineErr",& OneChPhiBaseLineErr);



  tall->Branch("ZPos",&OneChZPos);
  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZDataQual",& OneChZDataQual);  
  tall->Branch("ZPosErr",& OneChZPosErr);
  tall->Branch("ZTopWidthErr",& OneChZTopWidthErr);
  tall->Branch("ZSigmaErr",& OneChZSigmaErr);
  tall->Branch("ZHeightErr",& OneChZHeightErr);
  tall->Branch("ZBaseLineErr",& OneChZBaseLineErr);



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
    OneChPhiPosErr=SiPMPhiPosErr[iCh];
    OneChPhiPosDesign=SiPMPhiPosDesign[iCh];
    OneChPhiDataQual=SiPMPhiDataQual[iCh];
    OneChPhiTopWidthErr=SiPMPhiTopWidthErr[iCh];
    OneChPhiSigmaErr=SiPMPhiSigmaErr[iCh];
    OneChPhiHeightErr=SiPMPhiHeightErr[iCh];
    OneChPhiBaseLineErr=SiPMPhiBaseLineErr[iCh];
	//SiPMAllPhiPosGap[iCh]=OneChPhiPos-OneChPhiPosDesign;
	//Z Position
	OneChZPos=SiPMZPos[iCh];

	OneChZPosDesign=SiPMZPosDesign[iCh];
	//	OneChZPosGap=OneChZPos-OneChZPosDesign;
    OneChZDataQual=SiPMZDataQual[iCh];
	//SiPMAllZPosGap[iCh]=OneChZPos-OneChZPosDesign;
    OneChZPosErr=SiPMZPosErr[iCh];
    OneChZTopWidthErr=SiPMZTopWidthErr[iCh];
    OneChZSigmaErr=SiPMZSigmaErr[iCh];
    OneChZHeightErr=SiPMZHeightErr[iCh];
    OneChZBaseLineErr=SiPMZBaseLineErr[iCh];
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;

}

void UCIRunAnalysis(int run, Bool_t PhiScan) {
  // TCanvas* canvas1 =new TCanvas("canvas1","fitting",900,600);
  /*-----Define rec tree to be read----*/
  TFile* frec;
  if(PhiScan == true ){
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/PhiScan/xray%06d_graphs.root", run),"READ");
  }else{	
	frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/ZScan/xray%06d_graphs.root", run),"READ");
  }

  std::cout<<"Read rec "<<run<<std::endl;

  // canvas1->Divide(6,6);
  TF1* FitFunc[32];
  TGraphErrors* grScaler[32];
  for(int i=0;i<32;i++){
	grScaler[i]=(TGraphErrors*)frec->Get(Form("mppc%d",i));
	TString graphname=grScaler[i]->GetTitle();
	stringstream ss;
	//char graphname[64]=grScaler[i]->GetTitle();
	//	std::cout<<graphname<<std::endl;
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

	SiPMPhiPosDesign[MPPCch]=GraphPhiPos;
	SiPMZPosDesign[MPPCch]=GraphZPos;
	TString fitfuncname;
	fitfuncname.Form("fit%d",i);
	FitFunc[i] = new TF1(fitfuncname,TwoGaus,-300,300,5);
	Double_t FitErr[5];
	Double_t MeasPos;
	Double_t TopWidth;
	Double_t sigma;
	Double_t Height;
	Double_t Baseline;
	if(PhiScan==true){
	  FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
	  FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,60);
	  FitFunc[i]->SetParLimits(1,0,5);
	  FitFunc[i]->SetParLimits(2,0,5);
	  FitFunc[i]->SetParLimits(3,0,1000);
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
	  SiPMPhiPosErr[MPPCch]=FitErr[0];
	  SiPMPhiTopWidthErr[MPPCch]=FitErr[1];
	  SiPMPhiSigmaErr[MPPCch]=FitErr[2];
	  SiPMPhiHeightErr[MPPCch]=FitErr[3];
	  SiPMPhiBaseLineErr[MPPCch]=FitErr[4];

	}else{
	  SiPMZPos[MPPCch]=MeasPos;	 
	  SiPMZDataQual[MPPCch]=DataQual(FitErr,PhiScan,MeasPos,GraphZPos);
	  SiPMZPosErr[MPPCch]=FitErr[0];
	  SiPMZTopWidthErr[MPPCch]=FitErr[1];
	  SiPMZSigmaErr[MPPCch]=FitErr[2];
	  SiPMZHeightErr[MPPCch]=FitErr[3];
	  SiPMZBaseLineErr[MPPCch]=FitErr[4];
	}

	//std::cout<< "Position:  "<<MeasPos<<"+-"<<MeasPosErr<<std::endl;
	//std::cout<< "Top Width:  "<< TopWidth <<"  sigma:  "<<sigma<<"  Baseline:  "<<Baseline<<std::endl;
	//canvas1->cd(i+1);
	//grScaler[i]->Draw("ap");
	//FitFunc[i]->SetLineColor(kRed);
	//FitFunc[i]->Draw("same");
  }



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
                                  
