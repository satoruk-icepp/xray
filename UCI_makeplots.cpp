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
#define Nfitparam 5

Double_t ZPosDevAllch[nMPPC];
Double_t ZChiSqAllch[nMPPC];
Bool_t ZMeasuredAllch[nMPPC];
Bool_t ZValidAllch[nMPPC];
Double_t ZPosErrAllch[nMPPC];


Bool_t PhiMeasuredAllch[nMPPC];
Bool_t PhiValidAllch[nMPPC];
Double_t PhiPosDevAllch[nMPPC];
Double_t PhiChiSqAllch[nMPPC];


Bool_t AllTrue[nMPPC];


void UCI_makeplots(){
  // gStyle->SetTitleOffset( 2,"XYZ");
  // gStyle->SetTitleSize( 0.03,"XYZ");
  // gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,800);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,800);
  TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  TCanvas* canvas4=new TCanvas("canvas4","Chi Square",600,600);
  TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  //TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_tg.root","READ");
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/xray_UCI_corr_tg.root","READ");
  //TTree *txray = (TTree*)frec->Get("uci");
  TTree *txray = (TTree*)frec->Get("txray");
  //TTree *txray = (TTree*)frec->Get("xrayac");

  TGraph* grChZPos=new TGraph();
  TGraph* grGapCor=new TGraph();
  TH1D* WidthHist=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",200,0,20);
  TGraph2D* grZDev=new TGraph2D();
  TGraph2D* grPhiDev=new TGraph2D();
  //WidthHist->SetStats(0); //非表示

  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  Int_t ChNum;

  Double_t PhiResult[Nfitparam];
  Double_t PhiErr[Nfitparam];
  Double_t PhiPosDesign;
  Double_t PhiChiSq;
  Bool_t PhiMeasured;

  Double_t ZResult[Nfitparam];
  Double_t ZErr[Nfitparam];
  Double_t ZPosDesign;
  Double_t ZChiSq;
  Bool_t ZMeasured;

  txray->SetBranchAddress("ChNum",&ChNum);

  txray->SetBranchAddress("PhiPosDesign",&PhiPosDesign);
  txray->SetBranchAddress("PhiChiSq",&PhiChiSq);
  txray->SetBranchAddress("PhiMeasured",&PhiMeasured);

  txray->SetBranchAddress("ZPosDesign",&ZPosDesign);
  txray->SetBranchAddress("ZChiSq",&ZChiSq);
  txray->SetBranchAddress("ZMeasured",&ZMeasured);

  for(int i=0;i<5;i++){
    txray->SetBranchAddress(Form("PhiFitResult%d",i),&PhiResult[i]);
    txray->SetBranchAddress(Form("PhiFitErr%d",i),&PhiErr[i]);
    txray->SetBranchAddress(Form("ZFitResult%d",i),&ZResult[i]);
    txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
  }

  Int_t N=txray->GetEntries();
  std::cout<<"All channels: "<<N<<std::endl;
  Double_t tmpzpos;
  Bool_t former=false;

  for(int iCh=0;iCh<nMPPC;iCh++){
    txray->GetEntry(iCh);
    Double_t ZPos=ZResult[0];
    Double_t PhiPos=PhiResult[0];
    ZChiSqAllch[iCh]=ZChiSq;
    ZMeasuredAllch[iCh]=ZMeasured;
    ZPosErrAllch[iCh]=ZResult[2];
    PhiChiSqAllch[iCh]=PhiChiSq;
    PhiMeasuredAllch[iCh]=PhiMeasured;

    Int_t chline=floor(iCh/NLine);
    for(int i=0;i<4;i++){
      if(chline>=CFRPOrigin[i]&&chline<CFRPOrigin[i+1]){
        ZPosDesign=ZPosDesign+CFRPGap[i];
      }
    }
    Bool_t ZDataQual=false;
    if(FitQual(ZErr,false)==true){
      if(DevQual(ZPos,ZPosDesign,false)==true){
        ZDataQual=true;
      }
    }
    Bool_t PhiDataQual=false;
    if(FitQual(PhiErr,true)==true){
      if(DevQual(PhiPos,PhiPosDesign,true)==true){
        PhiDataQual=true;
      }else{
        std::cout<<"Calculated position is deviated: "<<iCh<<"  row:  "<<iCh/NLine<<"  line:  "<<iCh%NLine<<std::endl;
      }
    }else{
      std::cout<< "Fit Quality Cut: "<<iCh<<"  row:  "<<iCh/NLine<<"  line:  "<<iCh%NLine<<std::endl;
    }

    if(ZMeasured==true&&ZDataQual==true&&std::abs(ZResult[0])<120){
      if(former==true){
        WidthHist->Fill(ZResult[0]-tmpzpos);
        //	std::cout<<"factor: "<<cos(theta)<<std::endl;
        if(ZResult[0]-tmpzpos>17.0){
          std::cout<<"iCh:  "<<iCh<<"  row:  "<<iCh/NLine<<"  line:  "<<iCh%NLine<<std::endl;		
        }
      }
      if(iCh%44!=21){
        former=true;
      }else{
        former=false;
      }
      tmpzpos=ZResult[0];
      ZPosDevAllch[iCh]=ZPos-ZPosDesign;
      ZValidAllch[iCh]=true;
      grChZPos->SetPoint(grChZPos->GetN(),ChNum,ZResult[0]-ZPosDesign);
      grZDev->SetPoint(grZDev->GetN(),ZPosDesign,PhiPosDesign,ZResult[0]-ZPosDesign);
    }else{
      former=false;
    }

    if(PhiMeasured==true&&PhiDataQual==true&&std::abs(ZResult[0])<120&&ZDataQual==true){
      PhiPosDevAllch[iCh]=PhiPos-PhiPosDesign;
      PhiValidAllch[iCh]=true;
      grPhiDev->SetPoint(grPhiDev->GetN(),ZPosDesign,PhiPosDesign,PhiResult[0]-PhiPosDesign);
    }else{
      former=false;
    }
    AllTrue[iCh]=true;
    if(ZChiSq>1000){
      //	  std::cout <<"Chi Square is larger than 1000 : "<<iCh<<std::endl;
    }

  }


  canvas1->cd();
  //TPaveText *ptPhi = new TPaveText(.2,.925,.8,.975);
  // ptPhi->AddText("#phi_{calc}-#phi_{design}");
  // ptPhi->Draw();
  // InnerGeometry(PhiPosGapAllch,-0.5,0.5);
  grPhiDev->SetTitle("#phi Deviation;Z_{nom}[mm];#phi_{nom}[deg];#phi_{calc}-#phi_{nom}[deg]");
  grPhiDev->GetXaxis()->SetLimits(-300,300);
  grPhiDev->GetYaxis()->SetLimits(100,250);
  grPhiDev->SetMaximum(0.4);
  grPhiDev->SetMinimum(-0.1);
  grPhiDev->SetMarkerStyle(21);
  grPhiDev->SetMarkerSize(1);
  //grPhiDev->Draw("pcol");
  InnerGeometry(PhiPosDevAllch,PhiMeasuredAllch,PhiValidAllch,-0.5,0.5);
  // canvas2->cd();
  // TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
  // ptZ->AddText("Z_{calc}-Z_{design}");
  // ptZ->Draw();
  // //InnerGeometryArrow(ZPosDevAllch,ZvalidAllch,-10.0,0.0);
  // InnerGeometry(ZPosDevAllch,-10.0,0.0);

  // canvas3->cd();
  // grChZPos->SetTitle("Gap from the design value;channel;Gap[mm]");
  // grChZPos->Draw("ap");
  canvas2->cd();
  //InnerGeometry(ZChiSqAllch,0,2000);
  //grZDev->SetTitle("#Phi Deviation;Z_{nom};#phi_{nom};Z_{calc}-Z_{nom}");
  grZDev->SetTitle("Z Deviation;Z_{nom}[mm];#phi_{nom}[deg];Z_{calc}-Z_{nom}[mm]");
  grZDev->GetXaxis()->SetLimits(-300,300);
  grZDev->GetYaxis()->SetLimits(100,250);
  grZDev->SetMarkerStyle(20);
  grZDev->SetMaximum(0);
  grZDev->SetMinimum(-8);
  //  grZDev->Draw("colz");
  InnerGeometry(ZPosDevAllch,ZMeasuredAllch,ZValidAllch,-10.0,0.0);

  canvas3->cd();
  InnerGeometry(PhiChiSqAllch,PhiMeasuredAllch,AllTrue ,0.0,2000.0);

  canvas4->cd();
  InnerGeometry(ZPosErrAllch,ZMeasuredAllch,AllTrue ,0.0,20.0);
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
