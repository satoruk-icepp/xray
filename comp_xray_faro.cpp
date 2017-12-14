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

void Rotate3D(Double_t *oripos, Double_t *rotpos,Int_t axis);
void PlusOffset(Double_t *position,Double_t *offset);
void EulerRotation(Double_t alpha,Double_t beta, Double_t gamma,Double_t *position);
TMatrixD Rotation(Int_t axis, Double_t angle);
void comp_xray_faro(){
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","deformation",600,600);
  TCanvas* canvas4=new TCanvas("canvas4","direct comparison",600,600);
  TCanvas* canvas5=new TCanvas("canvas5","Z deviation",600,600);

  TString xraydatapath = "$(MEG2SYS)/analyzer/macros/xec/xray/";
  TString xrayfilename = "xray_UCI_corr_tg.root";
  TString xrayrootfile = xraydatapath + xrayfilename;
  TFile *fxray = new TFile(xrayrootfile.Data(),"READ");
  TTree *txray = (TTree*)fxray->Get("txray");

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  //TString farofilename = "XECFAROMPPC.root";
  //TString farofilename = "FAROMPPC_ip_mesh.root";
  //TString farofilename = "FAROMPPC_ip_mesh_shrink.root";
  TString farofilename = "FAROMPPC_ip_mesh_polar_shrink.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  //TTree* tfaro = (TTree*)ffaro->Get("faro");
  TTree* tfaro = (TTree*)ffaro->Get("tip");

  TFile* fout= new TFile("CompXF_shrink.root","recreate");


  TGraph* grChZPos = new TGraph();
  TGraph* grXrayFaroCor = new TGraph();
  TGraph* grXFPhi1DPhi= new TGraph();
  TGraph* grXFZ1DPhi = new TGraph();
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
    Double_t FaroPos[3]={FaroXPos,FaroYPos,FaroZPos};
    Double_t Offset[3]={0,0,0.5};
    // Double_t FaroPosAR[3];
    // Rotate3D(FaroPos,FaroPosAR,0);
    PlusOffset(FaroPos,Offset);
    EulerRotation(-0.08,0.08,0,FaroPos);
    Int_t chline=floor(iCh/NLine);
    Double_t XrayZPos = ZResult[0];
    Double_t XrayPhiPos = PhiResult[0];
    Double_t FaroXPosAR=FaroPos[0];
    Double_t FaroYPosAR=FaroPos[1];
    Double_t FaroZPosAR=FaroPos[2];
    Double_t FaroPhiPos = XYZ2Phi(FaroXPosAR,FaroYPosAR,FaroZPosAR);
    TrueAllch[iCh]=true;

    Bool_t ZDataQual=JudgeQual(ZQualArray,false,ZMeasured,PhiMeasured,XrayZPos,XrayPhiPos,ZPosDesign,ZErr);
    Bool_t PhiDataQual=JudgeQual(PhiQualArray,true,ZMeasured,PhiMeasured,XrayZPos,XrayPhiPos,PhiPosDesign,PhiErr);

    if(PhiDataQual==true&&ZDataQual==true){
      //grXrayFaroCor->SetPoint(iCh,ZPos,FaroZPos);

      grChXrayZPos->SetPoint(grChXrayZPos->GetN(),iCh,XrayZPos);

      grXFPhi1DPhi->SetPoint(grXFPhi1DPhi->GetN(),FaroPhiPos,XrayPhiPos-FaroPhiPos);
      grXFZ1DPhi->SetPoint(grXFZ1DPhi->GetN(),FaroPhiPos,XrayZPos-FaroZPosAR);


      grFaroMap->SetPoint(grFaroMap->GetN(),FaroZPosAR,FaroPhiPos);
      grXrayMap->SetPoint(grXrayMap->GetN(),XrayZPos,XrayPhiPos);
    }
    grChFaroZPos->SetPoint(grChFaroZPos->GetN(),iCh,FaroZPosAR);
  }

  canvas3->cd();
  grXFPhi1DPhi->SetName("phidep");
  grXFPhi1DPhi->SetTitle("#Delta #phi=#phi_{Xray}-#phi_{Faro}[deg];#phi_{Faro}[deg];#Delta #phi[deg]");
  //grXFPhi1DPhi->SetMinimum(0);
  grXFPhi1DPhi->SetMarkerStyle(7);
  grXFPhi1DPhi->SetMarkerColor(kRed);
  grXFPhi1DPhi->Draw("ap");

  canvas4->cd();
  grXrayMap->SetTitle("MPPC Position;Z[mm];#phi[deg]");
  grXrayMap->SetMarkerStyle(7);
  grXrayMap->SetMarkerColor(kRed);
  grFaroMap->SetMarkerStyle(7);
  grFaroMap->SetMarkerColor(kBlue);

  grXrayMap->Draw("ap");
  grFaroMap->Draw("p same");

  fout->cd();
  grXFPhi1DPhi->Write();
  fout->Close();

  canvas5->cd();
  grXFZ1DPhi->SetTitle("#Delta Z=Z_{Xray}-Z_{Faro}[mm];#phi_{Faro}[deg];#Delta Z[mm]");
  grXFZ1DPhi->SetMaximum(3);
  grXFZ1DPhi->SetMinimum(-3);
  grXFZ1DPhi->SetMarkerStyle(7);
  grXFZ1DPhi->SetMarkerColor(kRed);
  grXFZ1DPhi->Draw("ap");

  //grFaro3D->Draw("p0");
  Double_t position[3]={1,0,0};

  // std::cout<<"matrix 0 0: " << Rotation(0,0)[0][0]<<std::endl;
  EulerRotation(10,0,0,position);
  std::cout<<"vector: "<<position[0]<<position[1]<<position[2]<<std::endl;
}

void Rotate3D(Double_t *oripos, Double_t *rotpos,Int_t axis){
   Double_t theta=0.08*TMath::DegToRad();


   Double_t Xrot[3][3]={
      {1,0,0},
      {0,TMath::Cos(theta),-TMath::Sin(theta)},
      {0,TMath::Sin(theta),TMath::Cos(theta)}
   };
   Double_t Yrot[3][3]={
      {TMath::Cos(theta),0,TMath::Sin(theta)},
      {0,1,0},
      {-TMath::Sin(theta),0,TMath::Cos(theta)}
   };
   Double_t Zrot[3][3]={
      {TMath::Cos(theta),-TMath::Sin(theta),0},
      {TMath::Sin(theta),TMath::Cos(theta),0},
      {0,0,1}
   };


   for(int i=0;i<3;i++){
      Double_t tmprot=0;
      for(int j=0;j<3;j++){
         tmprot+=Xrot[i][j]*oripos[j];
      }
      rotpos[i]=tmprot;
   }
}

void PlusOffset(Double_t *position,Double_t *offset){
   for(int i=0;i<3;i++){
      position[i]+=offset[i];
   }
}

void EulerRotation(Double_t alpha,Double_t beta, Double_t gamma,Double_t  *position){
   TVectorD *vector = new TVectorD(3);
   for(int i=0;i<3;i++){
      (*vector)[i]=position[i];
   }
   *vector*=Rotation(2,alpha*TMath::DegToRad());
   *vector*=Rotation(0,beta*TMath::DegToRad());
   *vector*=Rotation(2,gamma*TMath::DegToRad());
   for(int i=0;i<3;i++){
      position[i]=(*vector)[i];
   }
}

TMatrixD Rotation(Int_t axis, Double_t angle){
   TMatrixD matrix(3,3);
   switch(axis){
   case 0: //X axis
      matrix[0][0]=1;
      matrix[0][1]=0;
      matrix[0][2]=0;
      matrix[1][0]=0;
      matrix[1][1]=TMath::Cos(angle);
      matrix[1][2]=-TMath::Sin(angle);
      matrix[2][0]=0;
      matrix[2][1]=TMath::Sin(angle);
      matrix[2][2]=TMath::Cos(angle);
      break;
   case 1:
      matrix[0][0]=TMath::Cos(angle);
      matrix[0][1]=0;
      matrix[0][2]=-TMath::Sin(angle);
      matrix[1][0]=0;
      matrix[1][1]=1;
      matrix[1][2]=0;
      matrix[2][0]=TMath::Sin(angle);
      matrix[2][1]=0;
      matrix[2][2]=TMath::Cos(angle);
      break;
   case 2:
      matrix[0][0]=TMath::Cos(angle);
      matrix[0][1]=-TMath::Sin(angle);
      matrix[0][2]=0;
      matrix[1][0]=TMath::Sin(angle);
      matrix[1][1]=TMath::Cos(angle);
      matrix[1][2]=0;
      matrix[2][0]=0;
      matrix[2][1]=0;
      matrix[2][2]=1;
      break;
   }

   return matrix;
}
