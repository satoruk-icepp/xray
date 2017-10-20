#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#endif
#define nMPPC 4092
#define kNchforXray 272
#define NRow 93
#define NLine 44

void OneRunAnalysis(int run, int scan);
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */);
Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t TwoGaus(Double_t *x,Double_t *par);
/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

Double_t SiPMAllXPos[nMPPC];
Double_t SiPMAllYPos[nMPPC];
Double_t SiPMAllPhiPos[nMPPC];
Double_t SiPMAllPhiPosErr[nMPPC];
Double_t SiPMAllPhiPosDesign[nMPPC];
Double_t SiPMAllPhiPosGap[nMPPC];
Double_t SiPMAllZPos[nMPPC];
Double_t SiPMAllZPosErr[nMPPC];
Double_t SiPMAllZPosDesign[nMPPC];
Double_t SiPMAllZPosGap[nMPPC];
Bool_t SiPMAllDataQualPhi[nMPPC];
Bool_t SiPMAllDataQualZ[nMPPC];

Int_t PhiRunList[]={305000,305001,305002,305003,305004,305005,305006,//phi_configA
		    305031,305032,305033,305034,305035,305036,305037,//phi_configB
		    305051,305062,305063,305064,305069,305070,305071,//phi_configC
					305085,/*305133,*/305086,305087,305098,305108,//phi_configD
		    305080,305079,305248,305077,305074,//phi_configE
		    305120,305121,305125,305127,305128,//phi_configF
		    305145,305144,305143,305142,305141,//phi_configG
		    305132,305134,305135,305136,305137,//phi_configH
					305235,305234,305233,305232,305231,305243,305242,/*305244,*///phi_configI
		    305208,305209,305217,305223,305224,//phi_configJ
		    305159,305158,305150,305149,305148,//phi_configK
		    305179,305187,305191,305192,305193,//phi_configL
};

Int_t ZRunList[]={305311,305312,305313,305314,305315,305316,305317,305318,//z_configA
				  305322,305323,305324,305330,305331,305332,305333,305334,305341,//z_configB
				  305350,305351,305352,305353,305354,305355,305357,305358,//z_configC
				  305299,305300,305301,305302,305303,305304,305305,305306,//z_configD
				  305363,305365,305366,305367,305368,305370,305371,305372,//z_configE
				  305381,305382,305383,305384,305385,305389,305390,305391,//z_configF
				  305395,305396,305397,305398,305399,305400,305401,305402,//z_configG
				  305413,305414,305415,305416,305417,305418,305419,305420,//z_configH
				  305439,305440,305441,305442,305448,305451,305452,305453,//z_configI
				  305456,305457,305458,305459,305460,305461,305462,305466,//z_configJ
				  305472,305473,305474,305475,305476,305477,305478,//z_configK
				  305512,305513,305514,305515,305516,305517,305518//z_configL
};

Int_t PhiRunNum=sizeof(PhiRunList)/sizeof(PhiRunList[0]);
Int_t ZRunNum=sizeof(ZRunList)/sizeof(ZRunList[0]);

void xray_general(){
  //TCanvas* canvas1=new TCanvas("canvas1","map",600,600);

  TFile *fall =new TFile("$(MEG2SYS)/analyzer/x-ray/xray_allch_dg.root","RECREATE");
  TTree *tall =new TTree("xrayac","xrayac");
  //TGraphErrors* grPhiPosAllch =new TGraphErrors();
  Int_t ChNum;

  Double_t OneChXPos;
  Double_t OneChYPos;

  Double_t OneChPhiPos;
  Double_t OneChPhiPosDesign;
  Double_t OneChPhiPosErr;
  Double_t OneChPhiPosGap;
  Bool_t OneChDataQualPhi;

  Double_t OneChZPos;
  Double_t OneChZPosDesign;
  Double_t OneChZPosErr;
  Double_t OneChZPosGap;
  Bool_t OneChDataQualZ;

  tall->Branch("ChNum",&ChNum);

  tall->Branch("XPos",&OneChXPos);
  tall->Branch("YPos",&OneChYPos);

  tall->Branch("PhiPos",&OneChPhiPos);
  tall->Branch("PhiPosDesign",&OneChPhiPosDesign);
  tall->Branch("PhiPosErr",& OneChPhiPosErr);
  tall->Branch("DataQualPhi",& OneChDataQualPhi);

  tall->Branch("ZPos",&OneChZPos);
  tall->Branch("ZPosDesign",&OneChZPosDesign);
  tall->Branch("ZPosErr",& OneChZPosErr);
  tall->Branch("DataQualZ",& OneChDataQualZ);

  for(int i=0;i<PhiRunNum;i++){
    OneRunAnalysis(PhiRunList[i],0);
  }

  for(int i=0;i<ZRunNum;i++){
    OneRunAnalysis(ZRunList[i],1);
  }


  Int_t cnt=0;
  for(int iCh=0;iCh<nMPPC;iCh++){
    ChNum=iCh;

	OneChXPos=SiPMAllXPos[iCh];
	OneChYPos=SiPMAllYPos[iCh];
	//Phi Position
    OneChPhiPos=SiPMAllPhiPos[iCh];
    OneChPhiPosErr=SiPMAllPhiPosErr[iCh];
    OneChPhiPosDesign=SiPMAllPhiPosDesign[iCh];
    OneChDataQualPhi=SiPMAllDataQualPhi[iCh];
	SiPMAllPhiPosGap[iCh]=OneChPhiPos-OneChPhiPosDesign;
	//Z Position
	OneChZPos=SiPMAllZPos[iCh];
    OneChZPosErr=SiPMAllZPosErr[iCh];
    OneChZPosDesign=SiPMAllZPosDesign[iCh];
	OneChZPosGap=OneChZPos-OneChZPosDesign;
    OneChDataQualZ=SiPMAllDataQualZ[iCh];
	SiPMAllZPosGap[iCh]=OneChZPos-OneChZPosDesign;
	/*if(OneChDataQualPhi==true){
	//grPhiPosAllch->SetPoint(cnt,iCh,OneChPhiPos-OneChPhiPosDesign);
	//grPhiPosAllch->SetPointError(cnt,0,OneChPhiPosErr);
	cnt++;
	}*/
    tall->Fill();
  }
  fall->cd();
  tall->Write();
  fall->Close();
  //grPhiPosAllch->Draw("ap");
  //canvas1->cd();
  //InnerGeometry(SiPMAllPhiPosGap,5,-5);
  return;
}



//________________________________________________________________________________
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */) {
  return x == 0.0 && y == 0.0 ? 0.0 : 180 + TMath::ATan2(-y, -x) *180/ TMath::Pi();
}

Int_t arraysearch(std::vector<Double_t> array, Double_t value){
  Bool_t found=false;
  Int_t sizearr=array.size();
  int i=-1;
  if(sizearr!=0){
    for (i = 0; i < sizearr; i++) {
      if (array[i]==value) {
	found=true;
	break;
      }
    }
  }
  if(found==false){
    i=-1;  
  }
  return i;
}


Double_t TwoGaus(Double_t *x,Double_t *par){
  Double_t xx=x[0];
  Double_t mean=par[0];
  Double_t topwidth=par[1];
  Double_t meanfg=par[0]-par[1]/2;
  Double_t meansg=par[0]+par[1]/2;
  Double_t sigma=par[2];
  Double_t height=par[3];
  Double_t f;
  if(xx<meanfg){
	f=height*TMath::Gaus(xx,meanfg,sigma,true);
  }else if(xx>=meanfg&&xx<meansg){
	f= height*TMath::Gaus(meanfg,meanfg,sigma,true);
  }else{
	f= height*TMath::Gaus(xx,meansg,sigma,true);
  }
  return f;
}

Double_t ArbFunc(Double_t *x,Double_t *par){
  Float_t xx=x[0];
  //Double_t f = par[0]/sqrt(2*TMath::Pi()*par[1])*TMath::Exp(-pow(xx-par[2],2)/(2*par[1]));
  //par[0]:mean
  //par[1]:topwidth
  //par[2]:zerowidth-topwidth
  //par[3]:height
  Double_t topstart=par[0]-par[1]/2;
  Double_t topend=par[0]+par[1]/2;
  Double_t redge=topstart-par[2]/2;
  Double_t fedge=topend+par[2]/2;
  Double_t height=par[3];
  Double_t f;
  if(xx<redge||xx>=fedge){
    f=0;
  }else if(xx>=redge&&xx<topstart){
    Double_t ratio=std::abs((xx-redge)/(topstart-redge));
    f= height*ratio;
  }else if(xx>=topstart&&xx<topend){
    f=height;  
  }else{
    Double_t ratio=std::abs((xx-topend)/(fedge-topend));
    f=height*(1-ratio);
  }
  return f;
}

void OneRunAnalysis(int run, int scan) {
  Int_t SiPMchNum;
  Double_t SiPMPos;
  Double_t SiPMPosErr;
  Double_t SiPMPosDesign;
  Double_t ScalerHeight;
  Double_t ScalerTopWidth;
  Double_t ScalerGroundWidth;

  /*-----Define rec tree to be read----*/
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/rec%06d.root", run),"READ");
  TTree *rec = (TTree*)frec->Get("rec");
  Int_t nEvent = rec->GetEntries();
  std::cout<<"Read rec "<<run<<". "<<nEvent<<" scaler events."<<std::endl;

  /*-----Read PM RunHeader----*/
  TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");
  MEGXECPMRunHeader *pmrh = 0;

  Int_t valid[kNchforXray] ={};
  Int_t MPPCindex[kNchforXray] ={};
  Float_t MPPCXYZ[kNchforXray][3] ={};
  Float_t MPPCPhi[kNchforXray]={};
  Int_t ScalerCh =0;
  Double_t theta= -TMath::Pi()*0.0/180;
  for (Int_t iPM = 0; iPM < nMPPC; iPM++) {
    pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
    if (pmrh->GetDRSAddress() >= 0) {
      ScalerCh = pmrh->GetDRSChipID() * 8 + pmrh->GetDRSChannelOnChip();
      //define whether measurement is valid or not
      if (ScalerCh < 0 || ScalerCh > kNchforXray) continue;
      valid[ScalerCh] = 1;
      MPPCindex[ScalerCh] = iPM;
      MPPCXYZ[ScalerCh][0] = pmrh->GetXYZAt(0) * 10;
	  Double_t tmpy=pmrh->GetXYZAt(1) * 10;
	  Double_t tmpz=pmrh->GetXYZAt(2) * 10;
	  Double_t LocalR=sqrt(tmpy*tmpy+tmpz*tmpz);
	  Double_t LocalPhi = TMath::ATan2(tmpy,tmpz);
      MPPCXYZ[ScalerCh][1] = LocalR*sin(LocalPhi+theta);
	  MPPCXYZ[ScalerCh][2] = LocalR*cos(LocalPhi+theta);
      MPPCPhi[ScalerCh] = XYZ2Phi(MPPCXYZ[ScalerCh][0], MPPCXYZ[ScalerCh][1], MPPCXYZ[ScalerCh][2]);
    }
  }
  /*-----Define rec tree to be read----*/
  TClonesArray* recTRGScaler = new TClonesArray("MEGTRGScaler");
  TBranch *branchrectrgscaler = 0;
  branchrectrgscaler = rec->GetBranch("trgscaler.");
  if (branchrectrgscaler) branchrectrgscaler->SetAddress(&recTRGScaler);
  MEGXRAYData* recXRAYData = new MEGXRAYData();
  TBranch *branchrecxraydata = 0;
  branchrecxraydata = rec->GetBranch("xraydata.");
  if (branchrecxraydata) branchrecxraydata->SetAddress(&recXRAYData);

  gStyle->SetTitleW(2);
  TGraph* grBeamPosition= new TGraph();
  TGraph* grScaler[kNchforXray];
  TGraph* grScalerABLS[kNchforXray];
  TF1* FitFunc[kNchforXray];
  Double_t BaseLine[kNchforXray];
  Double_t TopLine[kNchforXray];
  TH1I* ScalerHist[kNchforXray];
  //TF1* BaseGaus[kNchforXray];
  Double_t ScalerAverage[kNchforXray];

  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    grScaler[iScalerCh] = new TGraph();
    grScalerABLS[iScalerCh] = new TGraph();
    TString histname,histtitle/*,basefuncname*/,fitfuncname;
    histname.Form("histo%d",iScalerCh);
    histtitle.Form("histo title%d;Trigger Rate[Hz];Data Points",iScalerCh);
    //basefuncname.Form("base%d",iScalerCh);
    fitfuncname.Form("fit%d",iScalerCh);
    ScalerHist[iScalerCh] = new TH1I(histname,histtitle,30,0,180);
    //BaseGaus[iScalerCh] = new TF1(basefuncname,"[0]*TMath::Gaus(x,[1],[2])");
    FitFunc[iScalerCh] = new TF1(fitfuncname,TwoGaus,-300,300,4);


  }

  /*-----Baseline determination----*/
  Bool_t HighScaler=false;
  Int_t Scaler_buf;
  Float_t BeamZ_buf, BeamPhi_buf;
  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {
    branchrectrgscaler->GetEntry(iEvent);
    branchrecxraydata->GetEntry(iEvent);
    // read the data quality
    if (recXRAYData->GetIsGood()) continue;
    Int_t nTRGScalerEntries = recTRGScaler->GetEntriesFast();
    if (nTRGScalerEntries < kNchforXray) return;
    for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      Scaler_buf = ((MEGTRGScaler*)(recTRGScaler->At(iScalerCh)))->GetScaler();
      ScalerHist[iScalerCh]->Fill(Scaler_buf);
	  ScalerAverage[iScalerCh]+=Scaler_buf;
    }
  }

  /*Fitting Section*/
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
	ScalerAverage[iScalerCh]=ScalerAverage[iScalerCh]/nEvent;
    if (valid[iScalerCh]){
      Double_t histmean=ScalerHist[iScalerCh]->GetMean();
      Double_t histRMS=ScalerHist[iScalerCh]->GetRMS();
      //BaseGaus[iScalerCh]->SetParameters(100,histmean,histRMS); 
      //BaseGaus[iScalerCh]->SetParLimits(1,histmean-5,histmean+5);
      //ScalerHist[iScalerCh]->Fit(Form("base%d",iScalerCh),"NQ","",histmean-10,histmean+10);
      //Double_t baseline=BaseGaus[iScalerCh]->GetParameter(1);
	  Double_t baseline=histmean;
      //Double_t baseRMS=BaseGaus[iScalerCh]->GetParameter(2);
      BaseLine[iScalerCh]=baseline;
    }
  }

  /*Scaler Graph After Baseline Subtraction*/
  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {
    branchrectrgscaler->GetEntry(iEvent);
    branchrecxraydata->GetEntry(iEvent);
    // read the data quality
    if (recXRAYData->GetIsGood()) continue;
    for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      Scaler_buf = ((MEGTRGScaler*)(recTRGScaler->At(iScalerCh)))->GetScaler();
      BeamZ_buf = recXRAYData->GetBeamZ();
      BeamPhi_buf = recXRAYData->GetBeamPhi();
      if (scan == 0) {
		//condition: there is no large gap beam and MPPC(phi direction)
		if (TMath::Abs(BeamZ_buf - MPPCXYZ[iScalerCh][2]) < 18. ){
		  grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamZ_buf, Scaler_buf);
		  grScalerABLS[iScalerCh]->SetPoint(grScalerABLS[iScalerCh]->GetN(), BeamPhi_buf, Scaler_buf-BaseLine[iScalerCh]);
		}
      } 
      else if (scan == 1) {//phi scan
		//condition: there is no large gap beam and MPPC(Z direction)
		if (TMath::Abs(BeamPhi_buf - MPPCPhi[iScalerCh]) < 1.5 && TMath::Abs(MPPCXYZ[iScalerCh][2])< 157 ){
		  grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamZ_buf, Scaler_buf);
		  grScalerABLS[iScalerCh]->SetPoint(grScalerABLS[iScalerCh]->GetN(), BeamZ_buf, Scaler_buf-BaseLine[iScalerCh]);
		}
      }
    }
  }
 
  /*Writing to File...*/
  //  fout->cd();
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    if (valid[iScalerCh]&&grScaler[iScalerCh]->GetN()>10){
      SiPMchNum=MPPCindex[iScalerCh];
	  SiPMAllXPos[SiPMchNum]=MPPCXYZ[iScalerCh][0];
	  SiPMAllYPos[SiPMchNum]=MPPCXYZ[iScalerCh][1];
      if(scan==0){
		FitFunc[iScalerCh]->SetRange(MPPCPhi[iScalerCh]-15,MPPCPhi[iScalerCh]+15);
		FitFunc[iScalerCh]->SetParLimits(0,MPPCPhi[iScalerCh]-5,MPPCPhi[iScalerCh]+5);
		FitFunc[iScalerCh]->SetParLimits(1,0,2);
		FitFunc[iScalerCh]->SetParLimits(2,0,2);
		FitFunc[iScalerCh]->SetParLimits(3,0,80);
		FitFunc[iScalerCh]->SetParameters(MPPCPhi[iScalerCh],0.6,0.1,30);
		grScalerABLS[iScalerCh]->Fit(Form("fit%d",iScalerCh),"MNQ");
		//Fit Quality
		Double_t TmpFitErr;
		Double_t ErrorLimit[4]={1,1,1,10};
		Bool_t FineFit=true;
		Double_t ErrSquare=0;
		for(int i=0;i<4;i++){
		  TmpFitErr=FitFunc[iScalerCh]->GetParError(i);
		  if(TmpFitErr>ErrorLimit[i]){
			FineFit=false;
		  }
		}

		SiPMPosDesign=MPPCPhi[iScalerCh];
		SiPMAllPhiPosDesign[SiPMchNum]=MPPCPhi[iScalerCh];
		if(FineFit==true&&SiPMAllDataQualPhi[SiPMchNum]==false){
		  SiPMAllPhiPos[SiPMchNum]=FitFunc[iScalerCh]->GetParameter(0);
		  SiPMAllPhiPosErr[SiPMchNum]=FitFunc[iScalerCh]->GetParError(0);
		}
		if(abs(SiPMAllPhiPos[SiPMchNum]-SiPMAllPhiPosDesign[SiPMchNum])<10&&FineFit==true){
		  SiPMAllDataQualPhi[SiPMchNum]=true;	
		}else{
		  SiPMAllDataQualPhi[SiPMchNum]=false;	
		}
      }else if(scan==1){
		FitFunc[iScalerCh]->SetRange(MPPCXYZ[iScalerCh][2]-50,MPPCXYZ[iScalerCh][2]+50);
		FitFunc[iScalerCh]->SetParLimits(0,MPPCXYZ[iScalerCh][2]-15,MPPCXYZ[iScalerCh][2]+15);	
		FitFunc[iScalerCh]->SetParLimits(1,0,15);
		FitFunc[iScalerCh]->SetParLimits(2,0,25);
		FitFunc[iScalerCh]->SetParLimits(3,0,80);
		FitFunc[iScalerCh]->SetParameters(MPPCXYZ[iScalerCh][2],10,18,40);
		grScalerABLS[iScalerCh]->Fit(Form("fit%d",iScalerCh),"MNQ");
		SiPMPosDesign=MPPCXYZ[iScalerCh][2];
		SiPMAllZPosDesign[SiPMchNum]=MPPCXYZ[iScalerCh][2];
		Double_t TmpFitErr;
		Double_t ErrorLimit[4]={1,1,1,10};
		Bool_t FineFit=true;
		Double_t ErrSquare=0;
		for(int i=0;i<4;i++){
		  TmpFitErr=FitFunc[iScalerCh]->GetParError(i);
		  if(TmpFitErr>ErrorLimit[i]){
			FineFit=false;
		  }
		}

		if(FineFit==true&&SiPMAllDataQualZ[SiPMchNum]==false){
		  SiPMAllZPos[SiPMchNum]=FitFunc[iScalerCh]->GetParameter(0);
		  SiPMAllZPosErr[SiPMchNum]=FitFunc[iScalerCh]->GetParError(0);
		}
		if(abs(SiPMAllZPos[SiPMchNum]-SiPMAllZPosDesign[SiPMchNum])<10&&FineFit==true){
		  SiPMAllDataQualZ[SiPMchNum]=true;	
		}else{
		  SiPMAllDataQualZ[SiPMchNum]=false;	
		}
      }
      SiPMPos=FitFunc[iScalerCh]->GetParameter(0);
      SiPMPosErr=FitFunc[iScalerCh]->GetParError(0);
      ScalerHeight=FitFunc[iScalerCh]->GetParameter(3);
      ScalerTopWidth=FitFunc[iScalerCh]->GetParameter(1);
      ScalerGroundWidth=ScalerTopWidth+FitFunc[iScalerCh]->GetParameter(2);
    }
  }

  return;
}
