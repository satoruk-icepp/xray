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

Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];
Bool_t TrueAllch[nMPPC];
Bool_t ZValidAllch[nMPPC];
Bool_t PhiValidAllch[nMPPC];

void comp_xray_faro(){
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","deformation",600,600);
  TCanvas* canvas4=new TCanvas("canvas4","direct comparison",600,600);

  TString xraydatapath = "$(MEG2SYS)/analyzer/macros/xec/xray/";
  TString xrayfilename = "xray_UCI_corr_tg.root";
  TString xrayrootfile = xraydatapath + xrayfilename;
  TFile *fxray = new TFile(xrayrootfile.Data(),"READ");
  TTree *txray = (TTree*)fxray->Get("txray");

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  //TString farofilename = "XECFAROMPPC.root";
  //TString farofilename = "FAROMPPC_ip_mesh.root";
  //TString farofilename = "FAROMPPC_ip_mesh_shrink.root";
  TString farofilename = "FAROMPPC_ip_mesh_polar.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  //TTree* tfaro = (TTree*)ffaro->Get("faro");
  TTree* tfaro = (TTree*)ffaro->Get("tip");

  TFile* fout= new TFile("CompXF_shrink.root","recreate");


  TGraph* grChZPos = new TGraph();
  TGraph* grXrayFaroCor = new TGraph();
  TGraph2D* grXFZ2D = new TGraph2D();
  TGraph2D* grXFPhi2D = new TGraph2D();
  TGraph* grXFPhi1DPhi= new TGraph();
  TGraph* grChXrayZPos= new TGraph();
  TGraph* grChFaroZPos= new TGraph();
  TGraph2D* grFaro3D = new TGraph2D();
  TGraph2D* grXray3D = new TGraph2D();
  TGraph* grXrayMap = new TGraph();
  TGraph* grFaroMap = new TGraph();

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

  for(int i=0;i<4;i++){
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
  //tfaro->SetBranchAddress("DataQual",&FaroDataQual);

  Int_t N=txray->GetEntries();
  std::cout<<"All channels: "<<N<<std::endl;
  Int_t ZQualArray[4]={};
  Int_t PhiQualArray[4]={};
  for(int iCh=0;iCh<nMPPC;iCh++){
    txray->GetEntry(iCh);
    tfaro->GetEntry(iCh);
    Int_t chline=floor(iCh/NLine);
    Double_t XrayZPos = ZResult[0];
    Double_t XrayPhiPos = PhiResult[0];
    Double_t FaroPhiPos = XYZ2Phi(FaroXPos,FaroYPos,FaroZPos);
    TrueAllch[iCh]=true;

    Bool_t ZDataQual=JudgeQual(ZQualArray,false,ZMeasured,PhiMeasured,XrayZPos,XrayPhiPos,ZPosDesign,ZErr);
    Bool_t PhiDataQual=JudgeQual(PhiQualArray,true,ZMeasured,PhiMeasured,XrayZPos,XrayPhiPos,PhiPosDesign,PhiErr);

    if(PhiDataQual==true&&ZDataQual==true){
      //grXrayFaroCor->SetPoint(iCh,ZPos,FaroZPos);
      grXFZ2D->SetPoint(grXFZ2D->GetN(),FaroZPos,FaroPhiPos,XrayZPos-FaroZPos);
      ZPosGapAllch[iCh]=XrayZPos-FaroZPos;
      ZValidAllch[iCh]=true;
      grChXrayZPos->SetPoint(grChXrayZPos->GetN(),iCh,XrayZPos);

      grXFPhi1DPhi->SetPoint(grXFPhi1DPhi->GetN(),FaroPhiPos,XrayPhiPos-FaroPhiPos);
      grXFPhi2D->SetPoint(grXFPhi2D->GetN(),FaroZPos,FaroPhiPos,XrayPhiPos-FaroPhiPos);
      PhiPosGapAllch[iCh]=XrayPhiPos-FaroPhiPos;
      PhiValidAllch[iCh]=true;

      grFaroMap->SetPoint(grFaroMap->GetN(),FaroZPos,FaroPhiPos);
      grXrayMap->SetPoint(grXrayMap->GetN(),XrayZPos,XrayPhiPos);
    }
    grChFaroZPos->SetPoint(grChFaroZPos->GetN(),iCh,FaroZPos);
    grFaro3D->SetPoint(grFaro3D->GetN(),FaroXPos,FaroYPos,FaroZPos);
  }

  canvas1->cd();
  grXFPhi2D->SetMarkerStyle(20);
  grXFPhi2D->SetMinimum(-0.5);
  grXFPhi2D->SetMaximum(0.5);
  grXFPhi2D->SetTitle("#Phi Deviation;Z_{Xray}[mm];#phi_{Xray}[deg];#phi_{Xray}-#phi_{Faro}[deg]");
  //grXFPhi2D->Draw("pcol");

  InnerGeometry("#phi_{X-ray}-#phi_{Faro}[deg]",PhiPosGapAllch,PhiValidAllch,TrueAllch,-0.5,0.5);

  canvas2->cd();
  grXFZ2D->SetMarkerStyle(20);
  grXFZ2D->SetMaximum(2);
  grXFZ2D->SetMinimum(-2);
  grXFZ2D->SetTitle("Z Deviation;Z_{Xray}[mm];#phi_{Xray}[deg];Z_{Xray}-Z_{Faro}[mm]");
  //grXFZ2D->Draw("pcol");

  InnerGeometry("Z_{Xray}-Z_{Faro}[mm]",ZPosGapAllch,ZValidAllch,TrueAllch,-2,2);
  canvas3->cd();
  grXFPhi1DPhi->SetName("phidep");
  grXFPhi1DPhi->SetTitle("#phi deviation;#phi_{Faro}[deg];#phi_{X-ray}-#phi_{Faro}[deg]");
  //grXFPhi1DPhi->SetMinimum(0);
  grXFPhi1DPhi->SetMarkerStyle(20);
  grXFPhi1DPhi->SetMarkerColor(kRed);
  grXFPhi1DPhi->Draw("ap");

  canvas4->cd();
  grXrayMap->SetMarkerStyle(20);
  grXrayMap->SetMarkerColor(kRed);
  grFaroMap->SetMarkerStyle(20);
  grFaroMap->SetMarkerColor(kBlue);

  grXrayMap->Draw("ap");
  grFaroMap->Draw("p same");

  fout->cd();
  grXFPhi1DPhi->Write();
  fout->Close();

  //grFaro3D->Draw("p0");
}
