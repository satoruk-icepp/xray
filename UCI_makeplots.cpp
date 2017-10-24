#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "InnerGeometry.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];

void UCI_makeplots(){
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  TCanvas* canvas4=new TCanvas("canvas4","Gap Correlation",600,600);
  TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_TG.root","READ");
  TTree *xrayac = (TTree*)frec->Get("uci");

  TGraph* grChZPos=new TGraph();
  TGraph* grGapCor=new TGraph();
  TH1D* WidthHist=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",400,0,20);
  //WidthHist->SetStats(0); //非表示


  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};


  Int_t ChNum;

  Double_t PhiPos;
  Double_t PhiPosErr;
  Bool_t PhiDataQual;
  Double_t PhiPosDesign;
  Double_t ZPos;
  Double_t ZPosErr;
  Double_t ZPosDesign;
  Bool_t ZDataQual;

  xrayac->SetBranchAddress("ChNum",&ChNum);

  xrayac->SetBranchAddress("PhiPos",&PhiPos);
  xrayac->SetBranchAddress("PhiPosErr",&PhiPosErr);
  xrayac->SetBranchAddress("PhiDataQual",&PhiDataQual);
  xrayac->SetBranchAddress("PhiPosDesign",&PhiPosDesign);

  xrayac->SetBranchAddress("ZPos",&ZPos);
  xrayac->SetBranchAddress("ZPosErr",&ZPosErr);
  xrayac->SetBranchAddress("ZDataQual",&ZDataQual);
  xrayac->SetBranchAddress("ZPosDesign",&ZPosDesign);

  Int_t N=xrayac->GetEntries();
  std::cout<<"All channels: "<<N<<std::endl;
  Int_t cntz=0;
  Int_t cntphi=0;
  Double_t tmpzpos;
  Bool_t former=false;
  Double_t theta= -TMath::Pi()*0.0/180;

  for(int iCh=0;iCh<nMPPC;iCh++){
	xrayac->GetEntry(iCh);
	Int_t chline=floor(iCh/NLine);
	//std::cout<<"channel:  "<<iCh<<"  line:  "<<chline;
	Double_t XPosRealDesign;
	Double_t YPosRealDesign;
	Double_t ZPosRealDesign;
	Double_t PhiPosRealDesign;
	for(int i=0;i<4;i++){
	  if(chline>=CFRPOrigin[i]&&chline<CFRPOrigin[i+1]){
		//std::cout<<"CFRP:"<<i<<std::endl;
		ZPosDesign=ZPosDesign+CFRPGap[i];
		/*Double_t tmpy=YPosDesign;
		Double_t tmpz=ZPosDesign;
		Double_t LocalR=sqrt(tmpy*tmpy+tmpz*tmpz);
		Double_t LocalPhi = TMath::ATan2(tmpy,tmpz);
		XPosRealDesign=XPosDesign;
		YPosRealDesign=LocalR*sin(LocalPhi+theta);
		ZPosRealDesign=LocalR*cos(LocalPhi+theta);
		PhiPosRealDesign=XYZ2Phi(XPosRealDesign,YPosRealDesign,ZPosRealDesign);
		break;*/
	  }
	}
	if(PhiDataQual==true){
	  PhiPosGapAllch[iCh]=PhiPos-PhiPosDesign;
	  //AllPhiPosErr[iCh]=PhiPosErr;
	}else{
	  PhiPosGapAllch[iCh]=-100;
	}
	if(ZDataQual==true&&std::abs(ZPos)<120){
	  if(former==true){
		WidthHist->Fill((ZPos-tmpzpos)/cos(theta));
		//	std::cout<<"factor: "<<cos(theta)<<std::endl;
		if((ZPos-tmpzpos)/cos(theta)>17.0){
		  std::cout<<"iCh:  "<<iCh<<"  row:  "<<iCh/NLine<<"  line:  "<<iCh%NLine<<std::endl;		
		}

	  }
	  if(iCh%44!=21){
		former=true;
	  }else{
		former=false;
	  }
	  tmpzpos=ZPos;
	  ZPosGapAllch[iCh]=ZPos-ZPosDesign;
	  //AllZPosErr[iCh]=ZPosErr;
	  grChZPos->SetPoint(iCh,ChNum,ZPos-ZPosDesign);
	}else{
	  ZPosGapAllch[iCh]=-100;
	  former=false;
	}
	//grGapCor->SetPoint(iCh,AllZPosGap[iCh],AllPhiPosGap[iCh]);

  }


  canvas1->cd();
  TPaveText *ptPhi = new TPaveText(.2,.925,.8,.975);
  ptPhi->AddText("Gap between Designed Value in the Phi Direction");
  ptPhi->Draw();
  InnerGeometry(PhiPosGapAllch,-0.5,0.5);
  
  canvas2->cd();
  TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
  ptZ->AddText("Gap between Designed Value in the Z Direction");
  ptZ->Draw();
  InnerGeometry(ZPosGapAllch,-10.0,0.0);
  
  canvas3->cd();
  grChZPos->SetTitle("Gap from the design value;channel;Gap[mm]");
  grChZPos->Draw("ap");
  /*canvas4->cd();
  grGapCor->GetXaxis()->SetLimits(-10,10);
  grGapCor->SetMinimum(-10);
  grGapCor->SetMaximum(10);
  grGapCor->Draw("ap");
  */
  canvas5->cd();
  canvas5->SetGrid(0,0);
  gStyle->SetFuncColor(kRed);
  WidthHist->Fit("gaus");
  WidthHist->Draw();
}

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
