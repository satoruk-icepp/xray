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

  TGraph* Row3DplotXZ[NRow];
  TGraph* Row3DplotYZ[NRow];
  TF1* FitLine[NRow][2];

  Int_t Nwf=tin->GetEntries();
  std::cout<<"total channel: "<<Nwf<<std::endl;
  int NMppcRow[NRow]={};
  for(int Row=0; Row<NRow;Row++){
	Row3DplotXZ[Row]= new TGraph();
	Row3DplotYZ[Row]= new TGraph();
	for(int dir=0;dir<2;dir++){
	  FitLine[Row][dir] = new TF1(Form("fit%d_dir%d",Row,dir),"[0]*x+[1]");
	}
	//
  }
  for(int i=0;i<Nwf;i++ ){
	tin->GetEntry(i);
	Int_t Row = MPPCchannel/NLine;
	Int_t Column = MPPCchannel%NLine;
	NMppcRow[Row]+=1;
	Row3DplotXZ[Row]->SetPoint(Row3DplotXZ[Row]->GetN(),MPPCXPos,MPPCZPos);
	Row3DplotYZ[Row]->SetPoint(Row3DplotYZ[Row]->GetN(),MPPCYPos,MPPCZPos);
  }
  for(int Row=0; Row<NRow;Row++){
	//	std::cout<<"Row: "<<Row<<" MPPC: "<<NMppcRow[Row]<<std::endl;
	if(NMppcRow[Row]!=0){
	  Row3DplotXZ[Row]->Fit(Form("fit%d_dir%d",Row,0),"NQ");
	  Row3DplotYZ[Row]->Fit(Form("fit%d_dir%d",Row,1),"NQ");
	}
  }
  TCanvas* canvas1=new TCanvas("canvas1","visual",600,600);
  canvas1->cd();
  Int_t VisRow=7;
  Row3DplotXZ[VisRow]->SetMaximum(300);
  Row3DplotXZ[VisRow]->SetMinimum(-300);
  Row3DplotXZ[VisRow]->GetXaxis()->SetLimits(-800,0);
  Row3DplotXZ[VisRow]->SetMarkerStyle(20);
  Row3DplotXZ[VisRow]->SetMarkerSize(0.5);
  Row3DplotXZ[VisRow]->Draw("ap");
  FitLine[VisRow][0]->SetLineColor(kRed); 
   FitLine[VisRow][0]->Draw("same"); 
}
