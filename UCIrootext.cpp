#include "transcoordinate.h"
#include "DataQual.h"
#define Nfitparam 5
#define Nmppc 4092

void UCIrootext(){
  TCanvas* canvas1=new TCanvas("canvas1","position",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","deviation", 600,600);
  TCanvas* canvas3=new TCanvas("canvas3","phi", 600,600);
  TCanvas* canvas4=new TCanvas("canvas4","deviation", 1200,600);
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/xray_UCI_corr_tg.root","READ");
  TTree *txray = (TTree*)frec->Get("txray");
  TFile* fzin= new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/zscan_allruns_combined.root","read");
  TFile* fphiin= new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/phiscan_allruns_ntuple.root","read");
  TNtuple* zscan=  (TNtuple*)fzin->Get("ntup_comb");
  TNtuple* phiscan= (TNtuple*)fphiin->Get("ntup");
  TGraphErrors* grXrayZPhi= new TGraphErrors();
  TGraphErrors* grFaroZPhi= new TGraphErrors();
  TGraphErrors* grXrayZPhiTokyo= new TGraphErrors();
  // TH1F* phidev= new TH1F("phi deviation","phi deviation",100, -2,2);
  TGraphErrors* grPhiDev= new TGraphErrors();
  TH1D* hZDev_UCITokyo=new TH1D("z deviation","#Delta Z[mm]=Z_{UCI}-Z_{Tokyo};#Delta Z[mm];channels",200,-5,5);
  TH1D* hPhiDev_UCITokyo=new TH1D("phi deviation","#Delta #phi=#phi_{UCI}-#phi_{Tokyo}[deg];#Delta #phi[deg];channels",200,-1,1);
  Double_t Zcalcallch[Nmppc];
  Double_t Phicalcallch[Nmppc];
  Int_t Validallch[Nmppc]={};
  Float_t Zcalc;
  Float_t Zch;
  Float_t Phicalc;
  Float_t Phich;
  zscan->SetBranchAddress("zbkg",&Zcalc);
  zscan->SetBranchAddress("mppc",&Zch);
  phiscan->SetBranchAddress("phibkg",&Phicalc);
  phiscan->SetBranchAddress("mppc",&Phich);




  Int_t Nzin= zscan->GetEntries();
  Int_t Nphiin= phiscan->GetEntries();
  TH1F* distevmod= new TH1F("#Delta Z","space from the next MPPC;#Delta Z[mm];channels",100,13,17);
  TH1F* distodmev= new TH1F("a","space from the next MPPC;Z_{even}-Z_{odd}[mm];channels",100,13,17);
  for (int i = 0; i < Nzin; i++) {
    zscan->GetEntry(i);
    Validallch[(Int_t)Zch]+=1;
    if(TMath::Abs(Zcalc)<120.){
      Validallch[(Int_t)Zch]+=1;
      // std::cout<<"Z: "<<std::abs(Zcalc)<<std::endl;
    }
    Zcalcallch[(Int_t)Zch]=Zcalc;
  }
  for (int i = 0; i < Nphiin; i++) {
    phiscan->GetEntry(i);
    Validallch[(Int_t)Phich]+=1;
    Phicalcallch[(Int_t)Phich]=Phicalc;
  }
  for (int i = 0; i < Nmppc; i++) {
    if (Validallch[i]==3) {
      grXrayZPhi->SetPoint(grXrayZPhi->GetN(),Zcalcallch[i],Phicalcallch[i]);
      if(Validallch[i-1]==3){
        if (i%2==0) {
          distevmod->Fill(Zcalcallch[i]-Zcalcallch[i-1]);
        }else{
          distodmev->Fill(Zcalcallch[i]-Zcalcallch[i-1]);
        }
      }
    }
  }

  canvas2->cd();
  distevmod->SetLineColor(kRed);
  distevmod->Draw();
  distodmev->SetLineColor(kBlue);
  distodmev->Draw("same");

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  //TString farofilename = "XECFAROMPPC.root";
  //TString farofilename = "FAROMPPC_ip_mesh.root";
  //TString farofilename = "FAROMPPC_ip_mesh_shrink.root";
  TString farofilename = "FAROMPPC_ip_mesh_polar.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  //TTree* tfaro = (TTree*)ffaro->Get("faro");
  TTree* tfaro = (TTree*)ffaro->Get("tip");

  Int_t FaroChNum;
  Double_t FaroXPos;
  Double_t FaroYPos;
  Double_t FaroZPos;
  tfaro->SetBranchAddress("channel",&FaroChNum);
  tfaro->SetBranchAddress("XPos",&FaroXPos);
  tfaro->SetBranchAddress("YPos",&FaroYPos);
  tfaro->SetBranchAddress("ZPos",&FaroZPos);

  for (int i = 0; i < Nmppc; i++) {
    tfaro->GetEntry(i);
    Double_t FaroPhiPos = XYZ2Phi(FaroXPos,FaroYPos,FaroZPos);
    grFaroZPhi->SetPoint(grFaroZPhi->GetN(),FaroZPos,FaroPhiPos);
    if (Validallch[i]==3) {
      // phidev->Fill(Phicalcallch[i]-FaroPhiPos);
      grPhiDev->SetPoint(grPhiDev->GetN(),FaroPhiPos,Phicalcallch[i]-FaroPhiPos);
    }
  }


  canvas3->cd();
  grPhiDev->Draw("ap");

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
  Int_t ZQualArray[4]={};
  Int_t PhiQualArray[4]={};
  for (int i = 0; i < Nmppc; i++) {
    txray->GetEntry(i);
    Double_t ZPos=ZResult[0];
    Double_t PhiPos=PhiResult[0];
    Bool_t ZDataQual=JudgeQual(ZQualArray,false,ZMeasured,PhiMeasured,ZPos,PhiPos,ZPosDesign,ZErr);
    Bool_t PhiDataQual=JudgeQual(PhiQualArray,true,ZMeasured,PhiMeasured,ZPos,PhiPos,PhiPosDesign,PhiErr);
    if(ZDataQual==true&&PhiDataQual==true){
      grXrayZPhiTokyo->SetPoint(grXrayZPhiTokyo->GetN(),ZPos,PhiPos);
      if (Validallch[i]==3) {
        hZDev_UCITokyo->Fill(Zcalcallch[i]-ZPos);
        hPhiDev_UCITokyo->Fill(Phicalcallch[i]-PhiPos);
      }

    }
  }

  // canvas1->cd();
  // grXrayZPhi->SetTitle("MPPC position;Z[mm];#phi[deg]");
  // grXrayZPhi->SetMarkerStyle(7);
  // grXrayZPhi->SetMarkerSize(1);
  // grXrayZPhi->SetMarkerColor(kRed);
  // grXrayZPhi->Draw("ap");

  canvas1->cd();
  grXrayZPhiTokyo->SetMarkerStyle(7);
  grXrayZPhiTokyo->SetMarkerColor(kRed);
  grXrayZPhiTokyo->Draw("ap");

  canvas1->cd();
  grFaroZPhi->SetMarkerStyle(7);
  grFaroZPhi->SetMarkerColor(kBlue);
  grFaroZPhi->Draw("same p");



  canvas4->Divide(2,1);
  canvas4->cd(1);
  hZDev_UCITokyo->Draw();
  canvas4->cd(2);
  hPhiDev_UCITokyo->Draw();

}
