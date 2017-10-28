#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#include "DataQual.h"
#include "InnerGeometry.h"
#include "transcoordinate.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44
#define Nfitparam 5

void comp_xray_faro(){
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","deformation",600,600);

  TString xraydatapath = "$(MEG2SYS)/analyzer/macros/xec/xray/";
  TString xrayfilename = "xray_UCI_tg.root";
  TString xrayrootfile = xraydatapath + xrayfilename;
  TFile *fxray = new TFile(xrayrootfile.Data(),"READ");
  TTree *txray = (TTree*)fxray->Get("txray");

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  TString farofilename = "XECFAROMPPC.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  TTree* tfaro = (TTree*)ffaro->Get("faro");

  TGraph* grChZPos=new TGraph();
  TGraph* grXrayFaroCor=new TGraph();
  TGraph2D* grXFZ2D = new TGraph2D();
  TGraph2D* grXFPhi2D = new TGraph2D();
  TGraph* grXFPhi1DPhi= new TGraph();
  TH1D* WidthHist=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",400,0,20);
  //WidthHist->SetStats(0); //非表示


  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  /*xray tree */
  Int_t XRayChNum;
  Double_t PhiResult[Nfitparam];
  Double_t PhiErr[Nfitparam];
  Double_t PhiPosDesign;
  Bool_t PhiMeasured;

  Double_t ZResult[Nfitparam];
  Double_t ZErr[Nfitparam];
  Double_t ZPosDesign;
  Bool_t ZMeasured;

  txray->SetBranchAddress("ChNum",&XRayChNum);

  txray->SetBranchAddress("PhiPosDesign",&PhiPosDesign);
  txray->SetBranchAddress("PhiMeasured",&PhiMeasured);
  txray->SetBranchAddress("ZPosDesign",&ZPosDesign);
  txray->SetBranchAddress("ZMeasured",&ZMeasured);

  for(int i=0;i<5;i++){
	txray->SetBranchAddress(Form("PhiFitResult%d",i),&PhiResult[i]);
	txray->SetBranchAddress(Form("PhiFitErr%d",i),&PhiErr[i]);
	txray->SetBranchAddress(Form("ZFitResult%d",i),&ZResult[i]);
	txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
  }

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

  Int_t N=txray->GetEntries();
  std::cout<<"All channels: "<<N<<std::endl;

  for(int iCh=0;iCh<nMPPC;iCh++){
	txray->GetEntry(iCh);
	tfaro->GetEntry(iCh);
	Int_t chline=floor(iCh/NLine);
	Double_t XrayZPos = ZResult[0];
	Double_t XrayPhiPos = PhiResult[0];
	Double_t FaroPhiPos = XYZ2Phi(FaroXPos,FaroYPos,FaroZPos);
	Bool_t ZDataQual=false;
	if(DataQual(ZErr,false,XrayZPos,ZPosDesign)==true){
	  ZDataQual=true;
	}


	Bool_t PhiDataQual=false;
	//if(PhiChiSq<5000){
	if(DataQual(PhiErr,true,XrayPhiPos,PhiPosDesign)==true){
	  PhiDataQual=true;
	}

	if(PhiMeasured==true&&ZMeasured==true&&PhiDataQual&&ZDataQual==true&&std::abs(XrayZPos)<150&&FaroDataQual==true){
	  //grChZPos->SetPoint(iCh,ZPos,ZPos-FaroZPos);
	  //grXrayFaroCor->SetPoint(iCh,ZPos,FaroZPos);
	  grXFZ2D->SetPoint(grXFZ2D->GetN(),XrayZPos,XrayPhiPos,XrayZPos-FaroZPos);
	}

	if(PhiMeasured==true&&ZMeasured==true&&FaroDataQual==true&&PhiDataQual==true){
	  grXFPhi1DPhi->SetPoint(grXFPhi1DPhi->GetN(),XrayPhiPos,XrayPhiPos-FaroPhiPos);
	  grXFPhi2D->SetPoint(grXFPhi2D->GetN(),XrayZPos,XrayPhiPos,XrayPhiPos-FaroPhiPos);
	}
  }

  canvas1->cd();
  grXFPhi2D->SetMarkerStyle(20);
  grXFPhi2D->SetMinimum(0);
  grXFPhi2D->SetTitle("#Phi Deviation;Z_{Xray}[mm];#phi_{Xray}[deg];#phi_{Xray}-#phi_{Faro}[deg]");
  grXFPhi2D->Draw("pcol");
  
  canvas2->cd();
  grXFZ2D->SetMarkerStyle(20);
  grXFZ2D->SetTitle("Z Deviation;Z_{Xray}[mm];#phi_{Xray}[deg];Z_{Xray}-Z_{Faro}[mm]");
  grXFZ2D->Draw("pcol");

  canvas3->cd();
  grXFPhi1DPhi->SetTitle("#phi deviation;#phi_{X-ray}[deg];#phi_{X-ray}-#phi_{Faro}[deg]");
  grXFPhi1DPhi->SetMinimum(0);
  grXFPhi1DPhi->SetMarkerStyle(20);
  grXFPhi1DPhi->SetMarkerColor(kRed);
  grXFPhi1DPhi->Draw("ap");
}
