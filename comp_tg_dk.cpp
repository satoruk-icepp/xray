#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "function.h"
//#include "InnerGeometry.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

Double_t ZPosXray[nMPPC];
Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];

void comp_tg_dk(){
  TCanvas* canvas1=new TCanvas("canvas1","Z diff",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Phi diff",600,600);
  // TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  // TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  // TCanvas* canvas4=new TCanvas("canvas4","Gap Correlation",600,600);
  // TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  TString datapath = "$(MEG2SYS)/analyzer/x-ray/";
  TString tgfilename = "xray_UCI_allch_tg.root";
  TString tgrootfile = datapath + tgfilename;
  TFile *ftg = new TFile(tgrootfile.Data(),"READ");
  TTree *ttg = (TTree*)ftg->Get("uci");

  TString itfilename = "xray_UCI_allch_dk.root";
  TString itrootfile = datapath+itfilename;
  TFile* fit = new TFile(itrootfile.Data(),"read");
  TTree* tit = (TTree*)fit->Get("uci");
  TH1D* Zdiff =new TH1D("Z difference","Z_{tg}-Z_{it};Z_{tg}-Z_{it}[mm];Channels",400,-2,2);
  TH1D* Phidiff =new TH1D("Phi difference","#phi_{tg}-#phi_{it};#phi_{tg}-#phi_{it}[deg];Channels",400,-0.2,0.2);
  /*
  TGraph* grChZPos=new TGraph();
  TGraph* grXrayitCor=new TGraph();
  TH1D* WidthHist=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",400,0,20);*/
  //WidthHist->SetStats(0); //非表示
  

  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  /*xray tree */
  Int_t tgChNum;

  Double_t tgPhiPos;
  //Double_t tgPhiPosErr;
  //Double_t PhiPosDesign;
  Double_t tgZPos;
  Bool_t ZMeasured;
  //Double_t ZPosErr;
  //Double_t ZPosDesign;

  ttg->SetBranchAddress("ChNum",&tgChNum);

  ttg->SetBranchAddress("PhiPos",&tgPhiPos);
  // ttg->SetBranchAddress("PhiPosErr",&PhiPosErr);
  // ttg->SetBranchAddress("PhiDataQual",&PhiDataQual);
  // ttg->SetBranchAddress("PhiPosDesign",&PhiPosDesign);

  ttg->SetBranchAddress("ZPos",&tgZPos);
  ttg->SetBranchAddress("ZMeasured",&ZMeasured);
  // ttg->SetBranchAddress("ZPosErr",&ZPosErr);
  // ttg->SetBranchAddress("ZDataQual",&ZDataQual);
  // ttg->SetBranchAddress("ZPosDesign",&ZPosDesign);

  /*faro tree*/
  Int_t itChNum;
  // Double_t itXPos;
  // Double_t itYPos;
  Double_t itPhiPos;
  Double_t itZPos;

  // Bool_t itDataQual;
  tit->SetBranchAddress("ChNum",&itChNum);
  tit->SetBranchAddress("PhiPos",&itPhiPos);
  // tit->SetBranchAddress("XPos",&itXPos);
  // tit->SetBranchAddress("YPos",&itYPos);
  tit->SetBranchAddress("ZPos",&itZPos);
  // tit->SetBranchAddress("DataQual",&itDataQual);

  for(int iCh=0;iCh<nMPPC;iCh++){
	ttg->GetEntry(iCh);
	tit->GetEntry(iCh);
	Int_t chline=floor(iCh/NLine);

	if(ZMeasured==true&&std::abs(tgZPos)<150){
	  Zdiff->Fill(tgZPos-itZPos);
	  Phidiff->Fill(tgPhiPos-itPhiPos);
	  // ZPosGapAllch[iCh]=ZPos-itZPos;
	  // AllZPosErr[iCh]=ZPosErr;
	  // grChZPos->SetPoint(iCh,ZPos,ZPos-itZPos);
	  // grXrayitCor->SetPoint(iCh,ZPos,itZPos);
	}// else{
	//   ZPosGapAllch[iCh]=-100;
	//   former=false;
	// }
	//grGapCor->SetPoint(iCh,AllZPosGap[iCh],AllPhiPosGap[iCh]);

  }


  /* canvas1->cd();
	 TPaveText *ptPhi = new TPaveText(.2,.925,.8,.975);
	 ptPhi->AddText("Gap between Designed Value in the Phi Direction");
	 ptPhi->Draw();
	 InnerGeometry(PhiPosGapAllch,-0.5,0.5);*/
  
  //  canvas2->cd();
  // TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
  // ptZ->AddText("Deviation between Designed Value in the Z Direction");
  // ptZ->Draw();
  // InnerGeometry(ZPosGapAllch,-10.0,10.0);
  
  // canvas3->cd();
  // grChZPos->SetMarkerStyle(22);
  // grChZPos->SetMarkerColor(2);
  // grChZPos->SetTitle("Deviation from the design value;Z_{Xray}[mm];Deviation[mm]");
  // grChZPos->Draw("ap");
  // canvas4->cd();
  // //grXrayit->GetXaxis()->SetLimits(-10,10);
  // //grXrayit->SetMinimum(-10);
  // //grXrayit->SetMaximum(10);
  // grXrayitCor->Draw("ap");

  
	canvas1->cd();
	canvas1->SetGrid(0,0);
	gStyle->SetFuncColor(kRed);
	Zdiff->Fit("gaus");
	Zdiff->Draw();

	canvas2->cd();
	gStyle->SetFuncColor(kRed);
	Phidiff->Fit("gaus");
	Phidiff->Draw();
}
