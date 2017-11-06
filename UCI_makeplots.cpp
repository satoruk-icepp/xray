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

  const int NPCB=NRow*2;

  TGraph* grChZPos=new TGraph();
  TGraph* grPhiZDev=new TGraph();
  TGraph* grGapCor=new TGraph();
  TGraphErrors* grRowZPos[NPCB];
  TGraph* grDMPPC= new TGraph();;
  TGraph* grBeforeQC= new TGraph();
  TGraph* grAfterQC= new TGraph();
  TF1* FitLine[NPCB];
  TH1D* Dhist= new TH1D("distance","Distance between MPPCs fit result;Distance[mm];PCBs",160,14.5,15.5);
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

  for(int i=0;i<NPCB;i++){
    grRowZPos[i]= new TGraphErrors();
    FitLine[i]=new TF1(Form("line%d",i),"[0]*x+[1]");
  }

  Int_t ZQuallowZ=0;
  Int_t ZQualwellfit=0;
  Int_t ZQualsmalldev=0;
  Int_t ZMeasuredMPPCs=0;
  Int_t PhiQuallowZ=0;
  Int_t PhiQualwellfit=0;
  Int_t PhiQualsmalldev=0;
  Int_t PhiMeasuredMPPCs=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    txray->GetEntry(iCh);
    Double_t ZPos=ZResult[0];
    Double_t PhiPos=PhiResult[0];
    ZChiSqAllch[iCh]=ZChiSq;
    ZMeasuredAllch[iCh]=ZMeasured;
    ZPosErrAllch[iCh]=ZResult[2];
    PhiChiSqAllch[iCh]=PhiChiSq;
    PhiMeasuredAllch[iCh]=PhiMeasured;
    if (ZMeasured==true&&PhiMeasured==true) {
      grBeforeQC->SetPoint(grBeforeQC->GetN(),ZPos,PhiPos);
    }

    Int_t chline=floor(iCh/NLine);
    for(int i=0;i<4;i++){
      if(chline>=CFRPOrigin[i]&&chline<CFRPOrigin[i+1]){
        ZPosDesign=ZPosDesign+CFRPGap[i];
      }
    }
    Bool_t ZDataQual=false;
    if(ZMeasured==true){
      ++ZMeasuredMPPCs;
      if(LargeZQual(ZPos,ZMeasured)==true){
        ++ZQuallowZ;
        if(FitQual(ZErr,false)==true){
          ++ZQualwellfit;
          if(DevQual(ZPos,ZPosDesign,false)==true){
            ++ZQualsmalldev;
            ZDataQual=true;
          }
        }
      }
    }
    Bool_t PhiDataQual=false;
    if(PhiMeasured==true){
      ++PhiMeasuredMPPCs;
      if(LargeZQual(ZPos,ZMeasured)==true){
        ++PhiQuallowZ;
        if(FitQual(PhiErr,true)==true){
          ++PhiQualwellfit;
          if(DevQual(PhiPos,PhiPosDesign,true)==true){
            ++PhiQualsmalldev;
            PhiDataQual=true;
          }
        }
      }
    }
    if(ZDataQual==true&&PhiDataQual==true){
      grAfterQC->SetPoint(grAfterQC->GetN(),ZPos,PhiPos);
    }
    if(ZDataQual==true&&PhiDataQual==true){
      Int_t indexPCB;
      Int_t indexgraph;
      if(iCh%44<22){
        indexPCB=2*chline;
        indexgraph=grRowZPos[indexPCB]->GetN();
        grRowZPos[indexPCB]->SetPoint(indexgraph,iCh,ZResult[0]);
        grRowZPos[indexPCB]->SetPointError(indexgraph,0,ZErr[0]);
      }else{
        indexPCB=2*chline+1;
        indexgraph=grRowZPos[indexPCB]->GetN();
        grRowZPos[indexPCB]->SetPoint(indexgraph,iCh,ZResult[0]);
        grRowZPos[indexPCB]->SetPointError(indexgraph,0,ZErr[0]);
      }
      Double_t ZGap;
      if(former==true){
        ZGap=ZResult[0]-tmpzpos;
        if(iCh%2==0){
          WidthHist->Fill(ZGap);
        }
        //	std::cout<<"factor: "<<cos(theta)<<std::endl;
        if(ZGap>17.0){
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
      grChZPos->SetPoint(grChZPos->GetN(),ChNum,ZPos);
      grZDev->SetPoint(grZDev->GetN(),ZPosDesign,PhiPosDesign,ZResult[0]-ZPosDesign);
      if(PhiPos>120&&PhiPos<240){
        grPhiZDev->SetPoint(grPhiZDev->GetN(),PhiPos,ZPos-ZPosDesign);
      }
    }else{
      former=false;
    }

    if(PhiDataQual==true){
      PhiPosDevAllch[iCh]=PhiPos-PhiPosDesign;
      PhiValidAllch[iCh]=true;
      grPhiDev->SetPoint(grPhiDev->GetN(),ZPosDesign,PhiPosDesign,PhiResult[0]-PhiPosDesign);
    }else{
      former=false;
    }
    AllTrue[iCh]=true;


  }

  for(int i=0;i<NPCB;i++){
    grRowZPos[i]->Fit(Form("line%d",i),"Q");
    Double_t distance=FitLine[i]->GetParameter(0);
    grDMPPC->SetPoint(grDMPPC->GetN(),i,distance);
    Dhist->Fill(distance);
    std::cout<<"PCB: "<<i<<" distance: "<<distance<<" data points: "<<grRowZPos[i]->GetN()<<std::endl;
    TString grtitle;
    TString grtoptitle;
    if(i%2==0){
      grtoptitle.Form("PCB Row:%d,%s;",i/2,"US");
    }else{
      grtoptitle.Form("PCB Row:%d,%s;",i/2,"DS");
    }
    grtitle= grtoptitle+"Channel;Z Position[mm]";
    grRowZPos[i]->SetTitle(grtitle);
  }
  std::cout<<"Z Quality Cut"<<std::endl;
  std::cout<<"Measured: "<<ZMeasuredMPPCs<<std::endl;
  std::cout<<"Low Z: "<<ZQuallowZ<<std::endl;
  std::cout<<"Well-fit: "<<ZQualwellfit<<std::endl;
  std::cout<<"Small deviation: "<<ZQualsmalldev<<std::endl;

  std::cout<<"Phi Quality Cut"<<std::endl;
  std::cout<<"Measured: "<<PhiMeasuredMPPCs<<std::endl;
  std::cout<<"Low Z: "<<PhiQuallowZ<<std::endl;
  std::cout<<"Well-fit: "<<PhiQualwellfit<<std::endl;
  std::cout<<"Small deviation: "<<PhiQualsmalldev<<std::endl;

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
  grDMPPC->SetTitle("Distance between MPPCs;PCB;Distance[mm]");
  grDMPPC->SetMarkerStyle(20);
  //grDMPPC->Draw("ap");
  grBeforeQC->GetXaxis()->SetLimits(-160,160);
  grBeforeQC->SetMaximum(250);
  grBeforeQC->SetMinimum(110);
  grBeforeQC->SetMarkerStyle(21);
  grBeforeQC->SetMarkerSize(0.5);
  grBeforeQC->SetTitle("MPPC position;Z Position[mm];Phi Position[deg]");
  grBeforeQC->Draw("ap");
  grAfterQC->SetMarkerStyle(21);
  grAfterQC->SetMarkerSize(0.5);
  grAfterQC->SetMarkerColor(kRed);
  grAfterQC->Draw("same p");
  //InnerGeometry(PhiPosDevAllch,PhiMeasuredAllch,PhiValidAllch,-0.5,0.5);
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
  //InnerGeometry(ZPosDevAllch,ZMeasuredAllch,ZValidAllch,-10.0,0.0);
  Dhist->Draw();
  canvas3->cd();
  Int_t vispcb=47;
  grRowZPos[vispcb]->SetMarkerStyle(20);
  grRowZPos[vispcb]->SetMarkerColor(kRed);
  grRowZPos[vispcb]->Draw("ap");
  //InnerGeometry(PhiChiSqAllch,PhiMeasuredAllch,AllTrue ,0.0,2000.0);

  canvas4->cd();
  grPhiZDev->Draw("ap");
  //  InnerGeometry(ZPosErrAllch,ZMeasuredAllch,AllTrue ,0.0,20.0);
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
