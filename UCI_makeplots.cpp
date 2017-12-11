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
#define NCFRP 4
#define Nfitparam 5

Double_t ZPosDevAllch[nMPPC];
Double_t ZChiSqAllch[nMPPC];
Bool_t ZMeasuredAllch[nMPPC];
Bool_t ZValidAllch[nMPPC];
Double_t ZPosErrAllch[nMPPC];
Double_t ZPosAllch[NRow][NLine];

Bool_t PhiMeasuredAllch[nMPPC];
Bool_t PhiValidAllch[nMPPC];
Double_t PhiPosDevAllch[nMPPC];
Double_t PhiChiSqAllch[nMPPC];
Double_t PhiPosAllch[NRow][NLine];

Bool_t QualAllch[NRow][NLine];
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
   TCanvas* canvas6=new TCanvas("canvas6","Z deviation map",600,600);
   TCanvas* canvas7=new TCanvas("canavs7","Phi deviaiton map",600,600);
   //TFile *frec = new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/xray_raw_tg.root","READ");
   TFile *frec = new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/xray_UCI_corr_tg.root","READ");
   //TTree *txray = (TTree*)frec->Get("uci");
   TTree *txray = (TTree*)frec->Get("txray");
   //TTree *txray = (TTree*)frec->Get("xrayac");

   const int NPCB=NRow*2;

   TGraph* grChZPos=new TGraph();
   TGraph* grPhiZDev=new TGraph();
   TGraph* grGapCor=new TGraph();
   TGraph* grChDisOdd=new TGraph();
   TGraph* grChDisEven=new TGraph();
   TGraphErrors* grRowZPos[NPCB];
   TGraphErrors* grColumnZPos[NLine][4];
   TGraph* grDMPPC= new TGraph();
   TGraph* grBeforeQC= new TGraph();
   TGraph* grAfterQC= new TGraph();
   TGraph* grDAQRatePhi=new TGraph();
   TGraph* grDAQRateCorr=new TGraph();
   TF1* FitLine[NPCB];
   TH1D* Dhist= new TH1D("distance","Distance between MPPCs fit result;Distance[mm];PCBs",160,14.5,15.5);
   TH1D* WidthHistOdd=new TH1D("Z Spacing","#Delta Z;#Delta Z[mm];Channels",200,0,20);
   TH1D* WidthHistEven=new TH1D("Distance","Distance from the next MPPC;Distance[mm];Channels",200,0,20);
   TH1D* PhiDevHistEven= new TH1D("Phi Spacing Even","#Delta Phi;#Delta #phi[deg];channels",200,0,5);
   TH1D* PhiDevHistOdd= new TH1D("Phi Spacing Odd","#Delta Phi;#Delta #phi[deg];channels",200,0,5);
   TGraph2D* grZDev=new TGraph2D();
   TGraph2D* grPhiDev=new TGraph2D();
   //WidthHistOdd->SetStats(0); //非表示

   /*Define CFRP*/
   Int_t CFRPOrigin[NCFRP+1]={0,24,47,70,93};
   Double_t CFRPGap[NCFRP]={0,0.6,2.5,4.4};
   //Double_t CFRPGap[NCFRP]={0,0,0,0};

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

   for (int i = 0; i < NLine; i++) {
      for (int j = 0; j < NCFRP; j++) {
         grColumnZPos[i][j]=new TGraphErrors();
      }

   }

   Int_t ZQualArray[4]={};
   Int_t PhiQualArray[4]={};

   for(int iCh=0;iCh<nMPPC;iCh++){
      txray->GetEntry(iCh);
      Double_t ZPos=ZResult[0];
      Double_t PhiPos=PhiResult[0];
      Int_t tmprow=iCh/NLine;
      Int_t tmpcolumn=iCh/NLine;
      ZChiSqAllch[iCh]=ZChiSq;
      ZMeasuredAllch[iCh]=ZMeasured;
      ZPosErrAllch[iCh]=ZResult[2];
      PhiChiSqAllch[iCh]=PhiChiSq;
      PhiMeasuredAllch[iCh]=PhiMeasured;
      if (ZMeasured==true&&PhiMeasured==true) {
         grBeforeQC->SetPoint(grBeforeQC->GetN(),ZPos,PhiPos);
      }

      Int_t chrow=floor(iCh/NLine);
      Int_t chcolumn=iCh-NLine*chrow;
      for(int i=0;i<NCFRP;i++){
         if(chrow>=CFRPOrigin[i]&&chrow<CFRPOrigin[i+1]){
            ZPosDesign=ZPosDesign+CFRPGap[i];
         }
      }

      Bool_t ZDataQual=JudgeQual(ZQualArray,false,ZMeasured,PhiMeasured,ZPos,PhiPos,ZPosDesign,ZErr);
      Bool_t PhiDataQual=JudgeQual(PhiQualArray,true,ZMeasured,PhiMeasured,ZPos,PhiPos,PhiPosDesign,PhiErr);

      if(ZDataQual==true&&PhiDataQual==true){
         grDAQRateCorr->SetPoint(grDAQRateCorr->GetN(),ZResult[3],PhiResult[3]);
         grDAQRatePhi->SetPoint(grDAQRatePhi->GetN(),PhiPos,PhiResult[3]);
         grAfterQC->SetPoint(grAfterQC->GetN(),ZPos,PhiPos);
         for (int j = 0; j < NCFRP; j++) {
            if(chrow>=CFRPOrigin[j]&&chrow<CFRPOrigin[j+1]){
               grColumnZPos[chcolumn][j]->SetPoint(grColumnZPos[chcolumn][j]->GetN(),ZPos,PhiPos);
            }
         }
         DataQual[tmprow][tmpcolumn]=true;
         ZPosAllch[tmprow][tmpcolumn]=ZPos;
         PhiPosAllch[tmprow][tmpcolumn]=PhiPos;
      }
      if(ZDataQual==true){
         Int_t indexPCB;
         Int_t indexgraph;

         if(iCh%44<22){
            indexPCB=2*chrow;
            indexgraph=grRowZPos[indexPCB]->GetN();
            grRowZPos[indexPCB]->SetPoint(indexgraph,iCh,ZResult[0]);
            grRowZPos[indexPCB]->SetPointError(indexgraph,0,ZErr[0]);
         }else{
            indexPCB=2*chrow+1;
            indexgraph=grRowZPos[indexPCB]->GetN();
            grRowZPos[indexPCB]->SetPoint(indexgraph,iCh,ZResult[0]);
            grRowZPos[indexPCB]->SetPointError(indexgraph,0,ZErr[0]);
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

   for(int i=0;i<NRow;i++){
      for(int j=0;j<NLine;j++){
         Double_t ZGap;
         Double_t PhiGap;
         if(DataQual[i][j]==true){
            if(DataQual[i-1][j]==true){
               ZGap=ZPosAllch[i][j]-ZPosAllch[i-1][j];
               if(j%2==0&&j!=22){
                  WidthHistEven->Fill(ZGap);
               }else if(j%2==1){
                  WidthHistOdd->Fill(ZGap);
               }
            }
            if(j>CFRPOrigin[0]&&j<CFRPOrigin[1]){
               if(DataQual[i][j-1]==true){
                  PhiGap=PhiPosAllch[i][j]-PhiPosAllch[i][j-1];
                  if(j%2==0){
                     PhiDevHistEven->Fill(PhiGap);
                  }else{
                     PhiDevHistOdd->Fill(PhiGap);
                  }
               }
            }
         }
      }
   }
   for(int i=0;i<NPCB;i++){
      grRowZPos[i]->Fit(Form("line%d",i),"Q");
      Double_t distance=FitLine[i]->GetParameter(0);
      grDMPPC->SetPoint(grDMPPC->GetN(),i,distance);
      Dhist->Fill(distance);
      //std::cout<<"PCB: "<<i<<" distance: "<<distance<<" data points: "<<grRowZPos[i]->GetN()<<std::endl;
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

   canvas1->cd();

   canvas1->Divide(2,2);
   Int_t VisColumn=21;
   TF1* lin[NCFRP];
   for (int i = 0; i < NCFRP; i++) {
      canvas1->cd(i+1);
      lin[i]= new TF1(Form("pol1%d",i),"1./[0]*(x-[1])",-150,150);
      grColumnZPos[VisColumn][i]->SetTitle(Form("CFRP %d;Z Position[mm];#phi Position[deg]",i+1));
      grColumnZPos[VisColumn][i]->SetMarkerStyle(20);
      grColumnZPos[VisColumn][i]->SetMarkerColor(kBlue);
      grColumnZPos[VisColumn][i]->Fit(Form("pol1%d",i),"N");
      grColumnZPos[VisColumn][i]->Draw("ap");
      Double_t slope = lin[i]->GetParameter(0);
      Double_t slopeerr = lin[i]->GetParError(0);
      Double_t cut = lin[i]->GetParameter(1);
      Double_t cuterr = lin[i]->GetParError(1);
      Double_t Zatphi150= slope*150.+cut;
      Double_t Zatphi150err= slopeerr*150.+cuterr;
      Double_t Zatphi180= slope*180.+cut;
      Double_t Zatphi180err= slopeerr*180.+cuterr;
      Double_t Zatphi210= slope*210.+cut;
      Double_t Zatphi210err= slopeerr*210.+cuterr;
      std::cout<<" Phi150: "<<Zatphi150<<"+-"<<Zatphi150err<<std::endl;
      std::cout<<" Phi180: "<<Zatphi180<<"+-"<<Zatphi180err<<std::endl;
      std::cout<<" Phi210: "<<Zatphi210<<"+-"<<Zatphi210err<<std::endl;
      lin[i]->Draw("same");
   }

   std::cout<<"debug"<<lin[0]->Eval(0)<<std::endl;

   canvas2->cd();
   //TPaveText *ptPhi = new TPaveText(.2,.925,.8,.975);
   // ptPhi->AddText("#phi_{calc}-#phi_{design}");
   // ptPhi->Draw();
   // InnerGeometry(PhiPosGapAllch,-0.5,0.5);
   // grPhiDev->SetTitle("#phi Deviation;Z_{nom}[mm];#phi_{nom}[deg];#phi_{calc}-#phi_{nom}[deg]");
   // grPhiDev->GetXaxis()->SetLimits(-300,300);
   // grPhiDev->GetYaxis()->SetLimits(100,250);
   // grPhiDev->SetMaximum(0.4);
   // grPhiDev->SetMinimum(-0.1);
   // grPhiDev->SetMarkerStyle(21);
   // grPhiDev->SetMarkerSize(1);
   // //grPhiDev->Draw("pcol");
   // grDMPPC->SetTitle("Distance between MPPCs;PCB;Distance[mm]");
   // grDMPPC->SetMarkerStyle(20);
   // //grDMPPC->Draw("ap");
   // grBeforeQC->GetXaxis()->SetLimits(-160,160);
   // grBeforeQC->SetMaximum(250);
   // grBeforeQC->SetMinimum(110);
   // grBeforeQC->SetMarkerStyle(21);
   // grBeforeQC->SetMarkerSize(0.5);
   // grBeforeQC->SetTitle("MPPC position;Z Position[mm];Phi Position[deg]");
   // grBeforeQC->Draw("ap");
   // grAfterQC->SetMarkerStyle(21);
   // grAfterQC->SetMarkerSize(0.5);
   // grAfterQC->SetMarkerColor(kRed);
   // grAfterQC->Draw("same p");
   // for(int i=0;i<NCFRP;i++){
   //   lin[i]->Draw("same");
   // }
   grDAQRatePhi->SetMarkerStyle(22);
   grDAQRatePhi->SetMarkerColor(kRed);
   grDAQRatePhi->Draw("ap");
   // grDAQRateCorr->SetMarkerStyle(22);
   // grDAQRateCorr->SetMarkerColor(kRed);
   // grDAQRateCorr->Draw("ap");
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


   //InnerGeometry(ZChiSqAllch,0,2000);
   //grZDev->SetTitle("#Phi Deviation;Z_{nom};#phi_{nom};Z_{calc}-Z_{nom}");
   grZDev->SetTitle("Z Deviation;Z_{nom}[mm];#phi_{nom}[deg];Z_{calc}-Z_{nom}[mm]");
   grZDev->GetXaxis()->SetLimits(-300,300);
   grZDev->GetYaxis()->SetLimits(100,250);
   grZDev->SetMarkerStyle(20);
   grZDev->SetMaximum(0);
   grZDev->SetMinimum(-8);
   //  grZDev->Draw("colz");



   //


   //Dhist->Draw();
   canvas3->cd();
   Int_t vispcb=47;
   grRowZPos[vispcb]->SetMarkerStyle(20);
   grRowZPos[vispcb]->SetMarkerColor(kRed);
   //grRowZPos[vispcb]->Draw("ap");
   //InnerGeometry(ZChiSqAllch,ZMeasuredAllch,AllTrue ,0.0,1000.0);
   grChDisOdd->SetTitle("Distance between MPPCs;#ChNum;Distance[mm]");
   grChDisOdd->SetMarkerStyle(20);
   grChDisOdd->SetMarkerColor(kRed);
   grChDisOdd->Draw("ap");
   grChDisEven->SetMarkerStyle(20);
   grChDisEven->SetMarkerColor(kBlue);
   grChDisEven->Draw("p same");
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
   // canvas5->SetGrid(0,0);
   // gStyle->SetFuncColor(kRed);
   WidthHistOdd->SetLineColor(kBlue);
   // WidthHistOdd->Fit("gaus");
   WidthHistOdd->Draw();
   WidthHistEven->SetLineColor(kRed);
   WidthHistEven->Draw("same");

   canvas6->cd();
   TString ZIGTitle="Z_{calc}-Z_{nom}[mm]";
   InnerGeometry(ZIGTitle,ZPosDevAllch,ZMeasuredAllch,ZValidAllch,-10.0,0.0);
   canvas7->cd();
   TString PhiIGTitle="#phi_{calc}-#phi_{nom}[deg]";
   InnerGeometry(PhiIGTitle,PhiPosDevAllch,PhiMeasuredAllch,PhiValidAllch,-0.5,0.5);
}
