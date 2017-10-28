#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "function.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

Double_t ZPosXray[nMPPC];
Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];

void make3Dplot(){
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  gStyle->SetPalette(1);
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  //TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  //TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  //TCanvas* canvas4=new TCanvas("canvas4","Gap Correlation",600,600);
  //TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  TGraph2D* geometry =new TGraph2D();

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  TString farofilename = "XECFAROMPPC.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  TTree* tfaro = (TTree*)ffaro->Get("faro");

  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  /*faro tree*/
  Int_t FaroChNum;
  Double_t FaroXPos;
  Double_t FaroYPos;
  Double_t FaroZPos;
  Bool_t FaroDataQual;
  tfaro->SetBranchAddress("channel",&FaroChNum);
  tfaro->SetBranchAddress("XPos",&FaroXPos);
  tfaro->SetBranchAddress("YPos",&FaroYPos);
  tfaro->SetBranchAddress("ZPos",&FaroZPos);
  tfaro->SetBranchAddress("DataQual",&FaroDataQual);
  Int_t Nwf=tfaro->GetEntries();
  for(int iCh=0;iCh<Nwf;iCh++){
	tfaro->GetEntry(iCh);
	if(FaroDataQual==true){
	  geometry->SetPoint(geometry->GetN(),FaroZPos,FaroXPos,FaroYPos);
	}


  }

  canvas1->cd();
  geometry->SetTitle("Faro + Leica Result;Z Position[mm];X Position[mm];Y Position[mm]");
  geometry->SetMaximum(800);
  geometry->SetMinimum(-800);
  geometry->GetXaxis()->SetLimits(-400,400);
  geometry->GetYaxis()->SetLimits(-800,800); 
  geometry->SetMarkerStyle(20);
  geometry->SetMarkerSize(0.2);
  geometry->Draw("pcol");
  /*
  
	canvas2->cd();
	TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
	ptZ->AddText("Deviation between Designed Value in the Z Direction");
	ptZ->Draw();
	InnerGeometry(ZPosGapAllch,-10.0,10.0);
  
	canvas3->cd();
	grChZPos->SetMarkerStyle(22);
	grChZPos->SetMarkerColor(2);
	grChZPos->SetTitle("Deviation from the design value;Z_{Xray}[mm];Deviation[mm]");
	grChZPos->Draw("ap");
	canvas4->cd();
	grGapCor->GetXaxis()->SetLimits(-10,10);
	grGapCor->SetMinimum(-10);
	grGapCor->SetMaximum(10);
	grGapCor->Draw("ap");
  */
  /*
	canvas5->cd();
	canvas5->SetGrid(0,0);
	gStyle->SetFuncColor(kRed);
	WidthHist->Fit("gaus");
	WidthHist->Draw();*/
}
