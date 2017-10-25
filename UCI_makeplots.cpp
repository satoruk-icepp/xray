#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "InnerGeometry.h"
#include "DataQual.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];
Double_t ZChiSqAllch[nMPPC];
Bool_t ZMeasuredAllch[nMPPC];

Bool_t DataQual(Double_t *FitErr, Bool_t PhiScan, Double_t Position, Double_t Design);

void UCI_makeplots(){
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  TCanvas* canvas4=new TCanvas("canvas4","Gap Correlation",600,600);
  TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_tg.root","READ");
  TTree *txray = (TTree*)frec->Get("uci");

  TGraph* grChZPos=new TGraph();
  TGraph* grGapCor=new TGraph();
  TH1D* WidthHist=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",400,0,20);
  //WidthHist->SetStats(0); //非表示

  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  Int_t ChNum;

  Double_t PhiPos;
  Double_t PhiPosDesign;
  Double_t PhiErr[5];
  Double_t PhiChiSq;
  Bool_t PhiMeasured;

  Double_t ZPos;
  Double_t ZPosDesign;
  Double_t ZErr[5];
  Double_t ZChiSq;
  Bool_t ZMeasured;

  txray->SetBranchAddress("ChNum",&ChNum);

  txray->SetBranchAddress("PhiPos",&PhiPos);
  txray->SetBranchAddress("PhiPosDesign",&PhiPosDesign);
  txray->SetBranchAddress("PhiChiSq",&PhiChiSq);
  txray->SetBranchAddress("PhiMeasured",&PhiMeasured);

  txray->SetBranchAddress("ZPos",&ZPos);
  txray->SetBranchAddress("ZPosDesign",&ZPosDesign);
  txray->SetBranchAddress("ZChiSq",&ZChiSq);
  txray->SetBranchAddress("ZMeasured",&ZMeasured);

  for(int i=0;i<5;i++){
	txray->SetBranchAddress(Form("PhiFitErr%d",i),&PhiErr[i]);
	txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
  }

  Int_t N=txray->GetEntries();
  std::cout<<"All channels: "<<N<<std::endl;
  Double_t tmpzpos;
  Bool_t former=false;

  for(int iCh=0;iCh<nMPPC;iCh++){
	txray->GetEntry(iCh);

	Int_t chline=floor(iCh/NLine);
	for(int i=0;i<4;i++){
	  if(chline>=CFRPOrigin[i]&&chline<CFRPOrigin[i+1]){
		ZPosDesign=ZPosDesign+CFRPGap[i];
	  }
	}
	Bool_t ZDataQual=false;
	if(DataQual(ZErr,false,ZPos,ZPosDesign)==true){
	  ZDataQual=true;
	}

	if(ZMeasured==true&&ZDataQual==true){
	  if(former==true){
		WidthHist->Fill(ZPos-tmpzpos);
		//	std::cout<<"factor: "<<cos(theta)<<std::endl;
		if(ZPos-tmpzpos>17.0){
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
	  ZChiSqAllch[iCh]=ZChiSq;
	  ZMeasuredAllch[iCh]=ZMeasured;
	  grChZPos->SetPoint(iCh,ChNum,ZPos-ZPosDesign);
	}else{
	  ZPosGapAllch[iCh]=-100;
	  ZChiSqAllch[iCh]=-100;
	  former=false;
	}

	Bool_t PhiDataQual=false;
	//if(PhiChiSq<5000){
	if(DataQual(PhiErr,true,PhiPos,PhiPosDesign)==true){
	  PhiDataQual=true;
	}
	//	}

	if(PhiMeasured==true&&PhiDataQual==true){
	  PhiPosGapAllch[iCh]=PhiPos-PhiPosDesign;
	}else{
	  PhiPosGapAllch[iCh]=-100;
	  former=false;
	}

	//  if(ZChiSq<0||ZChiSq>1000){
	// 	 std::cout<<"channel: "<<iCh<<" Chi Square: "<<ZChiSq<<std::endl;	

	// }
	//grGapCor->SetPoint(iCh,AllZPosGap[iCh],AllPhiPosGap[iCh]);

  }


  canvas1->cd();
  TPaveText *ptPhi = new TPaveText(.2,.925,.8,.975);
  ptPhi->AddText("#phi_{calc}-#phi_{design}");
  ptPhi->Draw();
  InnerGeometry(PhiPosGapAllch,-0.5,0.5);
  
  canvas2->cd();
  TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
  ptZ->AddText("Z_{calc}-Z_{design}");
  ptZ->Draw();
  //InnerGeometryArrow(ZPosGapAllch,ZMeasuredAllch,-10.0,0.0);
  InnerGeometry(ZPosGapAllch,-10.0,0.0);
  
  canvas3->cd();
  grChZPos->SetTitle("Gap from the design value;channel;Gap[mm]");
  grChZPos->Draw("ap");
  canvas4->cd();
  InnerGeometry(ZChiSqAllch,0,2000);
  /*
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
