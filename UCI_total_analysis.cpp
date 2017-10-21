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
Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t TwoGaus(Double_t *x,Double_t *par);
Double_t SiPMZPos[nMPPC];
Double_t SiPMZPosErr[nMPPC];
Bool_t SiPMZDataQual[nMPPC];
Double_t SiPMPhiPos[nMPPC];
Double_t SiPMPhiPosErr[nMPPC];
Bool_t SiPMPhiDataQual[nMPPC];
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
  //Double_t OneChPhiPosDesign;
  Double_t OneChPhiPosErr;
  //Double_t OneChPhiPosGap;
  Bool_t OneChDataQualPhi;

  Double_t OneChZPos;
  //Double_t OneChZPosDesign;
  Double_t OneChZPosErr;
  //Double_t OneChZPosGap;
  Bool_t OneChDataQualZ;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("PhiPos",&OneChPhiPos);
  // tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiPosErr",& OneChPhiPosErr);
  tall->Branch("DataQualPhi",& OneChDataQualPhi);

  tall->Branch("ZPos",&OneChZPos);
  //tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZPosErr",& OneChZPosErr);
  tall->Branch("DataQualZ",& OneChDataQualZ);  

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
    //OneChPhiPosDesign=SiPMPhiPosDesign[iCh];
    OneChDataQualPhi=SiPMPhiDataQual[iCh];
	//SiPMAllPhiPosGap[iCh]=OneChPhiPos-OneChPhiPosDesign;
	//Z Position
	OneChZPos=SiPMZPos[iCh];
    OneChZPosErr=SiPMZPosErr[iCh];
	//    OneChZPosDesign=SiPMAllZPosDesign[iCh];
	//	OneChZPosGap=OneChZPos-OneChZPosDesign;
    OneChDataQualZ=SiPMZDataQual[iCh];
	//SiPMAllZPosGap[iCh]=OneChZPos-OneChZPosDesign;
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
	TString fitfuncname;
	fitfuncname.Form("fit%d",i);
	FitFunc[i] = new TF1(fitfuncname,TwoGaus,-300,300,5);
	Double_t MeasPos;
	Double_t MeasPosErr;
	Double_t TopWidth;
	Double_t sigma;
	Double_t Baseline;
	if(PhiScan==true){
	  FitFunc[i]->SetRange(GraphPhiPos-10,GraphPhiPos+10);
	  FitFunc[i]->SetParameters(GraphPhiPos,0.7,0.2,80,40);
	  FitFunc[i]->SetParLimits(1,0,2);
	  FitFunc[i]->SetParLimits(2,0,2);
	  FitFunc[i]->SetParLimits(3,0,200);
	  FitFunc[i]->SetParLimits(4,0,100);
	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
	  


	}else{
	  FitFunc[i]->SetRange(GraphZPos-30,GraphZPos+30);
	  FitFunc[i]->SetParameters(GraphZPos,10,10,80,40);
	  FitFunc[i]->SetParLimits(1,0,15);
	  FitFunc[i]->SetParLimits(2,0,25);
	  FitFunc[i]->SetParLimits(3,0,200);
	  FitFunc[i]->SetParLimits(4,0,100);
	  grScaler[i]->Fit(Form("fit%d",i),"MNQ");
	}
	
	MeasPos = FitFunc[i]->GetParameter(0);
	MeasPosErr = FitFunc[i]->GetParError(0);
	TopWidth = FitFunc[i]->GetParameter(1);
	sigma = FitFunc[i]->GetParameter(2);
	Baseline = FitFunc[i]->GetParameter(4);

	if(PhiScan==true){
	  SiPMPhiPos[MPPCch]=MeasPos;
	  SiPMPhiPosErr[MPPCch]=MeasPosErr;
	  SiPMPhiDataQual[MPPCch]=true;
	}else{
	  SiPMZPos[MPPCch]=MeasPos;
	  SiPMZPosErr[MPPCch]=MeasPosErr;
	  SiPMZDataQual[MPPCch]=true;
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
                                  
