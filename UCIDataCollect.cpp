#include "RunList.h"
#include "UCIRunAnalysis.cpp"

Int_t PhiRunNum=sizeof(PhiRunList)/sizeof(PhiRunList[0]);
Int_t ZRunNum=sizeof(ZRunList)/sizeof(ZRunList[0]);

void UCIDataCollect(){
  TString foutpath = "$(MEG2SYS)/analyzer/x-ray/";
  TString foutname = "xray_UCI_tg.root";
  TString foutroot = foutpath + foutname;
  TFile *fall =new TFile(foutroot.Data(),"RECREATE");
  TTree *tall =new TTree("txray","txray");
  std::cout<<"Writing to : "<<foutname<<std::endl;
  Int_t ChNum;

  Double_t OneChPhiPosDesign;
  Double_t OneChPhiResult[Nfitparam];
  Double_t OneChPhiErr[Nfitparam];
  Double_t OneChPhiChiSq;
  Bool_t OneChPhiMeasured;

  Double_t OneChZPosDesign;
  Double_t OneChZResult[Nfitparam];
  Double_t OneChZErr[Nfitparam];
  Double_t OneChZChiSq;
  Bool_t OneChZMeasured;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiChiSq",&OneChPhiChiSq);
  tall->Branch("PhiMeasured",&OneChPhiMeasured);
  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("PhiFitResult%d",i),&OneChPhiResult[i]);
	tall->Branch(Form("PhiFitErr%d",i),&OneChPhiErr[i]);
  }

  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZChiSq",&OneChZChiSq);
  tall->Branch("ZMeasured",&OneChZMeasured);
  for(int i=0;i<Nfitparam;i++){
	tall->Branch(Form("ZFitResult%d",i),&OneChZResult[i]);
	tall->Branch(Form("ZFitErr%d",i),&OneChZErr[i]);
  }

  for(int i=0;i<PhiRunNum;i++){
    UCIRunAnalysis(PhiRunList[i],true,false);
  }

  for(int i=0;i<ZRunNum;i++){
    UCIRunAnalysis(ZRunList[i],false,false);
  }

  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;
	//Phi Position
    OneChPhiPosDesign=SiPMPhiPosDesign[iCh];
	OneChPhiChiSq=SiPMPhiChiSq[iCh];
	OneChPhiMeasured=SiPMPhiMeasured[iCh];
	for(int i=0;i<Nfitparam;i++){
	  OneChPhiResult[i]=SiPMPhiResult[iCh][i];
	  OneChPhiErr[i]=SiPMPhiErr[iCh][i];
	}

	//Z Position
	OneChZPosDesign=SiPMZPosDesign[iCh];
	OneChZChiSq=SiPMZChiSq[iCh];
	OneChZMeasured=SiPMZMeasured[iCh];
	for(int i=0;i<Nfitparam;i++){
	  OneChZResult[i]=SiPMZResult[iCh][i];
	  OneChZErr[i]=SiPMZErr[iCh][i];
	}
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  return;
}
