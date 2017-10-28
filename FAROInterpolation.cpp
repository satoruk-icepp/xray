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

  TGraph2D* Row3Dplot[NRow];
  TF2* FitLine[NRow];

  Int_t Nwf=tin->GetEntries();
  std::cout<<"total channel: "<<Nwf<<std::endl;
  int NMppcRow[NRow]={};
  for(int Row=0; Row<NRow;Row++){
	Row3Dplot[Row]= new TGraph2D();
	FitLine[Row] = new TF2(Form("fit%d",Row),TDLine,-800,800,-800,800,3);
  }
  for(int i=0;i<Nwf;i++ ){
	tin->GetEntry(i);
	Int_t Row = MPPCchannel/NLine;
	Int_t Column = MPPCchannel%NLine;
	NMppcRow[Row]+=1;
	Row3Dplot[Row]->SetPoint(Row3Dplot[Row]->GetN(),MPPCXPos,MPPCYPos,MPPCZPos);
  }
  for(int Row=0; Row<NRow;Row++){
	//	std::cout<<"Row: "<<Row<<" MPPC: "<<NMppcRow[Row]<<std::endl;
	if(Row3Dplot[Row]!=0){
	  Row3Dplot[Row]->Fit(Form("fit%d",Row),"NQ");
	}
  }
  TCanvas* canvas1=new TCanvas("canvas1","visual",600,600);
  canvas1->cd();
  Int_t VisRow=7;
  Row3Dplot[VisRow]->SetMaximum(300);
  Row3Dplot[VisRow]->SetMinimum(-300);
  Row3Dplot[VisRow]->GetXaxis()->SetLimits(-800,0);
  Row3Dplot[VisRow]->GetYaxis()->SetLimits(-800,800);
  Row3Dplot[VisRow]->SetMarkerStyle(20);
  Row3Dplot[VisRow]->SetMarkerSize(0.5);
  Row3Dplot[VisRow]->Draw("p0");
  // FitLine[VisRow]->Draw("same"); 
}
