#include "RawRunAnalysis.cpp"

void RawDataCollect(){
  //TCanvas* canvas1=new TCanvas("canvas1","map",600,600);

  TString foutpath = "$(MEG2SYS)/analyzer/macros/xec/xray/";
  TString foutname = "xray_raw_tg.root";
  TString foutroot = foutpath + foutname;
  TFile *fall =new TFile(foutroot.Data(),"RECREATE");
  TTree *tall =new TTree("txray","txray");
  //TGraphErrors* grPhiPosAllch =new TGraphErrors();
  std::cout<<"Writing to : "<<foutname<<std::endl;
  Int_t ChNum;

  Double_t XPos;
  Double_t YPos;

  Double_t PhiPosDesign;
  Double_t PhiFitResult[5];
  Double_t PhiFitErr[5];
  Double_t PhiChiSq;
  Bool_t PhiMeasured;

  Double_t ZPosDesign;
  Double_t ZFitResult[5];
  Double_t ZFitErr[5];
  Double_t ZChiSq;
  Bool_t ZMeasured;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("XPos",&XPos);
  tall->Branch("YPos",&YPos);

  tall->Branch("PhiPosDesign",&PhiPosDesign);
  tall->Branch("PhiChiSq",&PhiChiSq);
  tall->Branch("PhiMeasured",&PhiMeasured);

  tall->Branch("ZPosDesign",&ZPosDesign);
  tall->Branch("ZChiSq",&ZChiSq);
  tall->Branch("ZMeasured",&ZMeasured);

  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("PhiFitResult%d",i),&PhiFitResult[i]);
	tall->Branch(Form("ZFitResult%d",i),&ZFitResult[i]);
	tall->Branch(Form("PhiFitErr%d",i),&PhiFitErr[i]);
	tall->Branch(Form("ZFitErr%d",i),&ZFitErr[i]);
  }

  for(int i=0;i<PhiRunNum;i++){
    RawRunAnalysis(PhiRunList[i],true,false);
  }

  for(int i=0;i<ZRunNum;i++){
    RawRunAnalysis(ZRunList[i],false,false);
  }


  Int_t cnt=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;

	XPos=XPosAllch[iCh];
	YPos=YPosAllch[iCh];
	//Phi Position
    PhiPosDesign=PhiPosDesignAllch[iCh];
    PhiChiSq=PhiChiSqAllch[iCh];
	PhiMeasured=PhiMeasuredAllch[iCh];
	//Z Position
    ZPosDesign=ZPosDesignAllch[iCh];
    ZChiSq=ZChiSqAllch[iCh];
	ZMeasured=ZMeasuredAllch[iCh];
	for(int i=0;i<5;i++){
	  PhiFitResult[i]=PhiFitResultAllch[iCh][i];
	  ZFitResult[i]=ZFitResultAllch[iCh][i];
	  PhiFitErr[i]=PhiFitErrAllch[iCh][i];
	  ZFitErr[i]=ZFitErrAllch[iCh][i];
	}
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;
}
