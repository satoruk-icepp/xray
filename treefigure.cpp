
void treefigure(){
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_DK.root","READ");
  TTree *txray = (TTree*)frec->Get("uci");
  Double_t ZErr[5];
  Double_t ZPosErr;
  Double_t ZTopWidthErr;
  Double_t ZSigmaErr;
  Double_t ZHeightErr;
  Double_t ZBaseLineErr;
  txray->SetBranchAddress("ZPosErr",&ZPosErr);
  txray->SetBranchAddress("ZTopWidthErr",&ZTopWidthErr);
  txray->SetBranchAddress("ZSigmaErr",&ZSigmaErr);
  txray->SetBranchAddress("ZHeightErr",&ZHeightErr);
  txray->SetBranchAddress("ZBaseLineErr",&ZBaseLineErr);

  TH1D* ZErrhist[5];
  for(int i=0;i<5;i++){
	TString HistName;
	TString HistTitle;
	TString Unit;
	switch(i){
	case 0:
	  Unit="mm";
	case 1:
	}
	TString XLabel;
	TString YLabel="channels";
	ZErrhist[i]=new TH1D(Form("par[%d]",i),Form("par[%d];par[%d]",i,i));
  }

  Int_t Nch=txray->GetEntries();
  for(int i=0;i<Nch;i++){
	txray->GetEntry(i);
	  }

}
