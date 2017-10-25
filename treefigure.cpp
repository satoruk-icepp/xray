void treefigure(){
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_tg.root","READ");
  TTree *txray = (TTree*)frec->Get("uci");
  TCanvas* canvasZ= new TCanvas("canvasZ","canvasZ",600,600);
  TCanvas* canvasPhi= new TCanvas("canvasPhi","canvasPhi",600,600);
  Bool_t ZMeasured;
  Bool_t PhiMeasured;
  Double_t ZErr[5];
  Double_t PhiErr[5];
  txray->SetBranchAddress("ZMeasured",&ZMeasured);
  txray->SetBranchAddress("PhiMeasured",&PhiMeasured);
  for(int i=0;i<5;i++){
	txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
	txray->SetBranchAddress(Form("PhiFitErr%d",i),&PhiErr[i]);
  }

  TH1D* ZErrHist[5];
  TH1D* PhiErrHist[5];

  Double_t ZQCut[5]={0.5,2,1,100,0.5};
  Double_t ZHistMax[5]={1,4,2,300,1};
  Double_t PhiQCut[5]={0.05,0.2,0.1,8,0.4};
  Double_t PhiHistMax[5]={0.1,0.4,0.2,16,0.8};
  for(int i=0;i<5;i++){
	TString Parameter;
	Parameter.Form("Fit Error %d",i);
	TString ZHistName="Z"+Parameter;
	TString PhiHistName="Phi"+Parameter;
	TString ZHistTitle;
	TString PhiHistTitle;
	TString Unit;
	Int_t NBin=100;
	Double_t HistMin=0;
	switch(i){
	case 0:
	  Unit="mm/deg";
	  break;
	case 1:
	  Unit="mm/deg";
	  break;
	case 2:
	  Unit="mm/deg";
	  break;
	case 3:
	  Unit="Hz";
	  break;
	case 4:
	  Unit="Hz";
	  break;
	default:
	  Unit="";
	  break;
	}
	TString XLabel=Parameter+"["+Unit+"]";
	TString YLabel="channels";
	ZHistTitle=ZHistName+";"+XLabel+";"+YLabel;
	PhiHistTitle=PhiHistName+";"+XLabel+";"+YLabel;
	ZErrHist[i]=new TH1D(ZHistName,ZHistTitle,NBin,HistMin,ZHistMax[i]);
	PhiErrHist[i]=new TH1D(PhiHistName,PhiHistTitle,NBin,HistMin,PhiHistMax[i]);
  }

  Int_t Npass=0;
  Int_t Nch=txray->GetEntries();
  for(int i=0;i<Nch;i++){
	txray->GetEntry(i);
	for(int i=0;i<5;i++){
	  if(ZMeasured==true){
		ZErrHist[i]->Fill(ZErr[i]);
	  }
	  if(PhiMeasured==true){
		PhiErrHist[i]->Fill(PhiErr[i]);
	  }
	}
  }

  canvasZ->Divide(2,3);
  canvasPhi->Divide(2,3);
  for(int i=0;i<5;i++){
	canvasZ->cd(i+1);
	ZErrHist[i]->Draw();
	Double_t Zymax=ZErrHist[i]->GetMaximum();
	TLine* Zline=new TLine(ZQCut[i],0,ZQCut[i],Zymax);
	Zline->SetLineColor(kRed);
	Zline->Draw();
	canvasPhi->cd(i+1);
	PhiErrHist[i]->Draw();
	Double_t Phiymax=PhiErrHist[i]->GetMaximum();
	TLine* Philine=new TLine(PhiQCut[i],0,PhiQCut[i],Phiymax);
	Philine->SetLineColor(kRed);
	Philine->Draw();
  }

}
