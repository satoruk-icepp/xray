void treefigure(){
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_TG.root","READ");
  TTree *txray = (TTree*)frec->Get("uci");
  TCanvas* canvasZ= new TCanvas("canvasZ","canvasZ",600,600);
  TCanvas* canvasPhi= new TCanvas("canvasPhi","canvasPhi",600,600);

  Double_t ZErr[5];
  Double_t PhiErr[5];
  for(int i=0;i<5;i++){
	txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
	txray->SetBranchAddress(Form("PhiFitErr%d",i),&PhiErr[i]);
  }

  TH1D* ZErrHist[5];
  TH1D* PhiErrHist[5];

  Double_t ZQCut[5];
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
	Double_t HistMax;
	switch(i){
	case 0:
	  Unit="mm/deg";
	  HistMax=1;
	  ZQCut[i]=0.5;
	  break;
	case 1:
	  Unit="mm/deg";
	  HistMax=2;
	  ZQCut[i]=1;
	  break;
	case 2:
	  Unit="mm/deg";
	  HistMax=5;
	  ZQCut[i]=2;
	  break;
	case 3:
	  Unit="Hz";
	  HistMax=10;
	  ZQCut[i]=2;
	  break;
	case 4:
	  Unit="Hz";
	  HistMax=0.8;
	  ZQCut[i]=0.5;
	  break;
	default:
	  Unit="";
	  HistMax=1;
	  break;
	}
	TString XLabel=Parameter+"["+Unit+"]";
	TString YLabel="channels";
	ZHistTitle=ZHistName+";"+XLabel+";"+YLabel;
	PhiHistTitle=PhiHistName+";"+XLabel+";"+YLabel;
	ZErrHist[i]=new TH1D(ZHistName,ZHistTitle,NBin,HistMin,HistMax);
	PhiErrHist[i]=new TH1D(PhiHistName,PhiHistTitle,NBin,HistMin,HistMax);
  }

  Int_t Nch=txray->GetEntries();
  for(int i=0;i<Nch;i++){
	txray->GetEntry(i);
	for(int i=0;i<5;i++){
	  if(ZErr[i]!=0){
		ZErrHist[i]->Fill(ZErr[i]);
	  }
	  if(PhiErr[i]!=0){
		PhiErrHist[i]->Fill(PhiErr[i]);
	  }
	}
  }
  canvasZ->Divide(2,3);
  canvasPhi->Divide(2,3);
  for(int i=0;i<5;i++){
	canvasZ->cd(i+1);
	ZErrHist[i]->Draw();
	Double_t ymax=ZErrHist[i]->GetMaximum();
	TLine* line=new TLine(ZQCut[i],0,ZQCut[i],ymax);
	line->SetLineColor(kRed);
	line->Draw();
	canvasPhi->cd(i+1);
	PhiErrHist[i]->Draw();
  }

}
