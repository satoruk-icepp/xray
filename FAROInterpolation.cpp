#include "function.h"

#define NLine 44
#define NRow 93

void FAROInterpolation(){
  TString datapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  TString filename = "XECFAROMPPC_well_fitted.root";
  TString filepath= datapath + filename;
  TFile* fin = new TFile(filepath.Data(),"READ");
  TTree* tin = (TTree*)fin->Get("faro");
  //
  Int_t MPPCchannel;
  Double_t MPPCXPos;
  Double_t MPPCYPos;
  Double_t MPPCZPos;
  Bool_t DataQual;
  tin->SetBranchAddress("channel",&MPPCchannel);
  tin->SetBranchAddress("XPos",&MPPCXPos);
  tin->SetBranchAddress("YPos",&MPPCYPos);
  tin->SetBranchAddress("ZPos",&MPPCZPos);
  tin->SetBranchAddress("DataQual",&DataQual);

  TGraph2D* total3Dplot=new TGraph2D();;
  TGraph* Row3DplotZX[NRow][2];
  TGraph* Row3DplotZY[NRow][2];
  TGraph2D* Column3Dplot[NLine];
  TF1* FitLine[NRow][2];

  Int_t Nwf=tin->GetEntries();
  std::cout<<"total channel: "<<Nwf<<std::endl;
  int NMppcRow[NRow][2]={{}};
  for(int Row=0; Row<NRow;Row++){
    for(int dir=0;dir<2;dir++){
      Row3DplotZX[Row][dir]= new TGraph();
      Row3DplotZY[Row][dir]= new TGraph();
    }
    //
  }
  for(int Column=0;Column<NLine;Column++){
  Column3Dplot[Column]=new TGraph2D();
  }
  for(int i=0;i<Nwf;i++ ){
    tin->GetEntry(i);
    Int_t Row = MPPCchannel/NLine;
    Int_t Column = MPPCchannel%NLine;
    Int_t side;
    if(Column<22){
      side=0;
    }else{
      side =1;
    }
    NMppcRow[Row][side]+=1;
    Row3DplotZX[Row][side]->SetPoint(Row3DplotZX[Row][side]->GetN(),MPPCZPos,MPPCXPos);
    Row3DplotZY[Row][side]->SetPoint(Row3DplotZY[Row][side]->GetN(),MPPCZPos,MPPCYPos);
    Column3Dplot[Column]->SetPoint(Column3Dplot[Column]->GetN(),MPPCXPos,MPPCYPos,MPPCZPos);
    total3Dplot->SetPoint(total3Dplot->GetN(),MPPCXPos,MPPCYPos,MPPCZPos);
  }
  for(int Row=0; Row<NRow;Row++){
    TF1* f[2];
    for(int i=0;i<2;i++){
      f[i]=new TF1(Form("fit%d",i),"[0]*x+[1]");
    }
    for(int side=0;side<2;side++){
      std::cout<<"Row: "<<Row<<" side: "<<side<<" MPPC: "<<NMppcRow[Row][side]<<std::endl;
      if(NMppcRow[Row][side]>1){
        Row3DplotZX[Row][side]->Fit(Form("fit%d",0),"NQ");
        Row3DplotZY[Row][side]->Fit(Form("fit%d",1),"NQ");
        Double_t xs=f[0]->GetParameter(0);        
        Double_t xc=f[0]->GetParameter(1);
        Double_t ys=f[1]->GetParameter(0);
        Double_t yc=f[1]->GetParameter(1);
        Double_t phi= TMath::RadToDeg()*TMath::ATan(ys/xs);
        Double_t theta = TMath::RadToDeg()*TMath::ACos(1/sqrt(1+xs*xs+ys*ys));
        std::cout <<"Row: "<<Row<<" side: "<<side<<" theta: "<<theta<<" phi: "<<phi<<std::endl;
      }
    }
  }
  TCanvas* canvas1=new TCanvas("canvas1","visual",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","total",600,600);
  TCanvas* canvas3=new TCanvas("canvas3","Column",600,600);
  canvas1->cd();
  Int_t VisRow=77;
  Int_t Visside=0;

  Row3DplotZX[VisRow][Visside]->SetMarkerStyle(20);
  Row3DplotZX[VisRow][Visside]->SetMarkerSize(0.5);
  Row3DplotZX[VisRow][Visside]->Draw("ap");
  canvas2->cd();
  total3Dplot->SetMarkerStyle(20);
  total3Dplot->Draw("pcol");

canvas3->cd();
Int_t VisColumn=7;
Column3Dplot[VisColumn]->SetMarkerStyle(20);
Column3Dplot[VisColumn]->Draw("pcol");

  //FitLine[VisRow][0]->SetLineColor(kRed); 
  //FitLine[VisRow][0]->Draw("same");
  //std::cout << "par0: "<<par0<<" par1: "<<par1<<std::endl;
}
