void treefigure(){
  TFile *frec = new TFile("$(MEG2SYS)/analyzer/x-ray/xray_UCI_allch_DK.root","READ");
  TTree *txray = (TTree*)frec->Get("uci");
  TCanvas* canvasZ= new TCanvas("canvasZ","canvasZ",600,600);
  //TCanvas* canvasPhi= new TCanvas("canvasPhi","canvasPhi",600,600);

  Double_t ZErr[5];
  for(int i=0;i<5;i++){
	txray->SetBranchAddress(Form("ZFitErr%d",i),&ZErr[i]);
  }

  TH1D* ZErrHist[5];
  Double_t QCut[5];
  for(int i=0;i<5;i++){
	TString Parameter;
	Parameter.Form("Fit Error %d",i);
	TString HistName=Parameter;
	TString HistTitle;
	TString Unit;
	Int_t NBin=100;
	Double_t HistMin=0;
	Double_t HistMax;
	switch(i){
	case 0:
	  Unit="mm";
	  HistMax=1;
	  QCut[i]=0.5;
	  break;
	case 1:
	  Unit="mm";
	  HistMax=2;
	  QCut[i]=1;
	  break;
	case 2:
	  Unit="mm";
	  HistMax=5;
	  QCut[i]=2;
	  break;
	case 3:
	  Unit="Hz";
	  HistMax=10;
	  QCut[i]=2;
	  break;
	case 4:
	  Unit="Hz";
	  HistMax=0.8;
	  QCut[i]=0.5;
	  break;
	default:
	  Unit="";
	  HistMax=1;
	  break;
	}
	TString XLabel=Parameter+"["+Unit+"]";
	TString YLabel="channels";
	HistTitle=Parameter+";"+XLabel+";"+YLabel;
	ZErrHist[i]=new TH1D(HistName,HistTitle,NBin,HistMin,HistMax);
  }

  Int_t Nch=txray->GetEntries();
  for(int i=0;i<Nch;i++){
	txray->GetEntry(i);
	for(int i=0;i<5;i++){
	  if(ZErr[i]!=0){
		ZErrHist[i]->Fill(ZErr[i]);
	  }
	}
  }
  canvasZ->Divide(2,3);
  for(int i=0;i<5;i++){
	canvasZ->cd(i+1);
	ZErrHist[i]->Draw();
	Double_t ymax=ZErrHist[i]->GetMaximum();
	TLine* line=new TLine(QCut[i],0,QCut[i],ymax);
	line->SetLineColor(kRed);
	line->Draw();
  }

}
