#define Nmppc 4092
#include "transcoordinate.h"

void UCIrootext(){
  TCanvas* canvas1=new TCanvas("canvas1","position",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","deviation", 600,600);
    TCanvas* canvas3=new TCanvas("canvas3","phi", 600,600);
  TFile* fzin= new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/zscan_allruns_combined.root","read");
  TFile* fphiin= new TFile("$(MEG2SYS)/analyzer/macros/xec/xray/phiscan_allruns_ntuple.root","read");
  TNtuple* zscan=  (TNtuple*)fzin->Get("ntup_comb");
  TNtuple* phiscan= (TNtuple*)fphiin->Get("ntup");
  TGraphErrors* grXrayZPhi= new TGraphErrors();
  TGraphErrors* grFaroZPhi= new TGraphErrors();
  // TH1F* phidev= new TH1F("phi deviation","phi deviation",100, -2,2);
  TGraphErrors* grPhiDev= new TGraphErrors();
  Double_t Zcalcallch[Nmppc];
  Double_t Phicalcallch[Nmppc];
  Int_t Validallch[Nmppc]={};
  Float_t Zcalc;
  Float_t Zch;
  Float_t Phicalc;
  Float_t Phich;
  zscan->SetBranchAddress("zcalc",&Zcalc);
  zscan->SetBranchAddress("mppc",&Zch);
  phiscan->SetBranchAddress("phicalc",&Phicalc);
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
  canvas1->cd();
  grXrayZPhi->SetTitle("MPPC position;Z[mm];#phi[deg]");
  grXrayZPhi->SetMarkerStyle(7);
  grXrayZPhi->SetMarkerSize(1);
  grXrayZPhi->SetMarkerColor(kRed);
  grXrayZPhi->Draw("ap");
  canvas2->cd();
  distevmod->SetLineColor(kRed);
  distevmod->Draw();
  distodmev->SetLineColor(kBlue);
  distodmev->Draw("same");

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  //TString farofilename = "XECFAROMPPC.root";
  //TString farofilename = "FAROMPPC_ip_mesh.root";
  //TString farofilename = "FAROMPPC_ip_mesh_shrink.root";
  TString farofilename = "FAROMPPC_ip_mesh_polar_shrink.root";
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
  canvas1->cd();
  grFaroZPhi->SetMarkerStyle(7);
  grFaroZPhi->SetMarkerColor(kBlue);
  grFaroZPhi->Draw("same p");

  canvas3->cd();
  grPhiDev->Draw("ap");
}
