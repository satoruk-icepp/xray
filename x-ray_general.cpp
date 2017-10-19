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

void xray_ana(Int_t run, Int_t scan);
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */);
Int_t arraysearch(std::vector<Double_t> array, Double_t value);
Double_t ArbFunc(Double_t *x,Double_t *par);
Double_t SiPMAllZPos[nMPPC];
Double_t SiPMAllPhiPos[nMPPC];
Double_t SiPMAllZPosErr[nMPPC];
Double_t SiPMAllPhiPosErr[nMPPC];
Double_t SiPMAllPhiPosDesign[nMPPC];

/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/
const Int_t PhiRunNum=1;
Int_t PhiRunList[PhiRunNum]={305181};


void xray_general(){
  for(int i=0;i<PhiRunNum;i++){
    BLMeas(PhiRunList[i],0);
  }
}

void BLMeas(int run, int scan) {
  TFile *fout =new TFile(Form("xray_%06d_%01d.root",run,scan),"RECREATE");
  TTree *tout =new TTree("xray","xray");
  Int_t SiPMchNum;  
  Double_t SiPMPos;
  Double_t SiPMPosErr;
  Double_t SiPMPosDesign;
  Double_t ScalerHeight;
  Double_t ScalerTopWidth;
  Double_t ScalerGroundWidth;
  tout->Branch("SiPMPos",&SiPMPos);
  tout->Branch("SiPMchNum",&SiPMchNum);
  tout->Branch("SiPMPosErr",&SiPMPosErr);
  tout->Branch("SiPMPosDesign",&SiPMPosDesign);
  tout->Branch("ScalerHeight",&ScalerHeight);
  tout->Branch("ScalerTopWidth",&ScalerTopWidth);
  tout->Branch("ScalerGroundWidth",&ScalerGroundWidth);

  /*-----Define rec tree to be read----*/
  TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/x-ray/rec%06d.root", run),"READ");
  TTree *rec = (TTree*)frec->Get("rec");
  Int_t nEvent = rec->GetEntries();
  std::cout<<"Read rec "<<run<<". "<<nEvent<<" scaler events."<<std::endl;

  /*-----Read PM RunHeader----*/
  TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");
  MEGXECPMRunHeader *pmrh = 0;

  Int_t valid[kNchforXray] ={};
  /*position and index of each MPPC*/
  Int_t MPPCindex[kNchforXray] ={};
  Float_t MPPCXYZ[kNchforXray][3] ={};
  Float_t MPPCPhi[kNchforXray]={};
  Int_t ScalerCh =0;
  for (Int_t iPM = 0; iPM < nMPPC; iPM++) {
    pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
    if (pmrh->GetDRSAddress() >= 0) {
      //std::cout<<pmrh->GetDRSChipID()<<" "<<pmrh->GetDRSChannelOnChip()<<std::endl;
      ScalerCh = pmrh->GetDRSChipID() * 8 + pmrh->GetDRSChannelOnChip();

      //define whether measurement is valid or not
      if (ScalerCh < 0 || ScalerCh > kNchforXray) continue;
      valid[ScalerCh] = 1;
      MPPCindex[ScalerCh] = iPM;
      MPPCXYZ[ScalerCh][0] = pmrh->GetXYZAt(0) * 10;
      MPPCXYZ[ScalerCh][1] = pmrh->GetXYZAt(1) * 10;
      MPPCXYZ[ScalerCh][2] = pmrh->GetXYZAt(2) * 10;
      //converting the Cartesian coordinates to cylinderical coordinates
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
  //beam position
  TGraph* grBeamPosition= new TGraph();
  //Scaler graph
  TGraph* grScaler[kNchforXray];
  //Scaler graph after subtraction
  TGraph* grScalerABLS[kNchforXray];
  TF1* FitFunc[kNchforXray];
  Double_t BaseLine[kNchforXray];
  Double_t TopLine[kNchforXray];
  TH1I* ScalerHist[kNchforXray];
  TF1* BaseGaus[kNchforXray];
  TF1* TopGaus[kNchforXray];

  std::vector<std::vector<Int_t > > Multiplicity;
  std::vector<std::vector<Int_t > > ScalerSum;
  std::vector<std::vector<Double_t > > ScalerAve;
  std::vector<std::vector<Double_t > > StoreBeamZ;
  std::vector<std::vector<Double_t > > StoreBeamPhi;
  Multiplicity.resize(kNchforXray);
  ScalerSum.resize(kNchforXray);
  ScalerAve.resize(kNchforXray);
  StoreBeamZ.resize(kNchforXray);
  StoreBeamPhi.resize(kNchforXray);

  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    grScaler[iScalerCh] = new TGraph();
    grScalerABLS[iScalerCh] = new TGraph();
    TString histname,histtitle,basefuncname,topfuncname,fitfuncname;
    histname.Form("histo%d",iScalerCh);
    histtitle.Form("histo title%d;Trigger Rate[Hz];Data Points",iScalerCh);
    basefuncname.Form("base%d",iScalerCh);
    topfuncname.Form("top%d",iScalerCh);
    fitfuncname.Form("fit%d",iScalerCh);
    ScalerHist[iScalerCh] = new TH1I(histname,histtitle,30,0,180);
    BaseGaus[iScalerCh] = new TF1(basefuncname,"[0]*TMath::Gaus(x,[1],[2])");
    TopGaus[iScalerCh] = new TF1(topfuncname,"[0]*TMath::Gaus(x,[1],[2])");
    FitFunc[iScalerCh] = new TF1(fitfuncname,ArbFunc,-10,300,4);
    FitFunc[iScalerCh]->SetParLimits(1,0,15);
    FitFunc[iScalerCh]->SetParLimits(2,0,25);
    FitFunc[iScalerCh]->SetParLimits(3,0,80);

  }

  /*-----Baseline determination----*/
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
      BeamZ_buf = recXRAYData->GetBeamZ();
      BeamPhi_buf = recXRAYData->GetBeamPhi();
      Scaler_buf = ((MEGTRGScaler*)(recTRGScaler->At(iScalerCh)))->GetScaler();
      ScalerHist[iScalerCh]->Fill(Scaler_buf);
      Int_t as=arraysearch(StoreBeamZ[iScalerCh],BeamZ_buf);
      if(as==-1){
	ScalerSum[iScalerCh].push_back(Scaler_buf);
	Multiplicity[iScalerCh].push_back(0);
	if(scan==0){
	  StoreBeamPhi[iScalerCh].push_back(BeamPhi_buf);
	}else if(scan==1){
	  StoreBeamZ[iScalerCh].push_back(BeamZ_buf);
	}
      }else{
	Int_t mp=Multiplicity[iScalerCh][as];
	Multiplicity[iScalerCh][as]=mp+1;
	ScalerSum[iScalerCh][as]=ScalerSum[iScalerCh][as]+Scaler_buf;
      }
    }
  }

  /*Fitting Section*/
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    if (valid[iScalerCh]){
      Double_t histmean=ScalerHist[iScalerCh]->GetMean();
      Double_t histRMS=ScalerHist[iScalerCh]->GetRMS();
      BaseGaus[iScalerCh]->SetParameters(100,histmean,histRMS); 
      ScalerHist[iScalerCh]->Fit(Form("base%d",iScalerCh),"NQ","",20,120);
      Double_t baseline=BaseGaus[iScalerCh]->GetParameter(1);
      Double_t baseRMS=BaseGaus[iScalerCh]->GetParameter(2);
      BaseLine[iScalerCh]=baseline;
      TopGaus[iScalerCh]->SetParameters(100,baseline+40,10); 
      ScalerHist[iScalerCh]->Fit(Form("top%d",iScalerCh),"NQ","",baseline+2*baseRMS,180);
      TopLine[iScalerCh]=(TopGaus[iScalerCh]->GetParameter(1))-BaseLine[iScalerCh];
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
 
  /*Title description*/
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    if (!valid[iScalerCh]) {
      grScaler[iScalerCh]->SetTitle("No MPPC connected.");
    } else {
      grScaler[iScalerCh]->SetTitle(Form("WD: %d, MPPC: %04d (@phi: %2.1lf deg, z: %2.1lf mm)", iScalerCh, MPPCindex[iScalerCh], MPPCPhi[iScalerCh], MPPCXYZ[iScalerCh][2]));
    }
  }

  /*Writing to File...*/
  fout->cd();
  for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
    if (valid[iScalerCh]&&grScaler[iScalerCh]->GetN()>10){
      SiPMchNum=MPPCindex[iScalerCh];
      if(scan==0){
	FitFunc[iScalerCh]->SetParameters(MPPCPhi[iScalerCh],1.2,1.5,40);
	grScalerABLS[iScalerCh]->Fit(Form("fit%d",iScalerCh),"NQ");
	SiPMPosDesign=MPPCPhi[iScalerCh];
	SiPMAllPhiPos[SiPMchNum]=FitFunc[iScalerCh]->GetParameter(0);
	SiPMAllPhiPosErr[SiPMchNum]=FitFunc[iScalerCh]->GetParError(0);
	SiPMAllPhiPosDesign[SiPMchNum]=MPPCPhi[iScalerCh];
      }else if(scan==1){
	FitFunc[iScalerCh]->SetParameters(MPPCXYZ[iScalerCh][2],10,18,40);
	grScalerABLS[iScalerCh]->Fit(Form("fit%d",iScalerCh),"NQ");
	SiPMPosDesign=MPPCXYZ[iScalerCh][2];
      }
      SiPMPos=FitFunc[iScalerCh]->GetParameter(0);
      SiPMPosErr=FitFunc[iScalerCh]->GetParError(0);
      ScalerHeight=FitFunc[iScalerCh]->GetParameter(3);
      ScalerTopWidth=FitFunc[iScalerCh]->GetParameter(1);
      ScalerGroundWidth=ScalerTopWidth+FitFunc[iScalerCh]->GetParameter(2);
      std::cout<<"iScalerCh:  "<<iScalerCh<<"  SiPM channel:  "<<SiPMchNum<<"  position:  "<<SiPMPos<<"+-"<<SiPMPosErr<<std::endl;
      grScaler[iScalerCh]->Write(Form("WDch%03d",iScalerCh));
      grScalerABLS[iScalerCh]->Write(Form("WDch%03d_ABLS",iScalerCh));
      ScalerHist[iScalerCh]->Write();
      tout->Fill();
    }
  }
  fout->cd();
  tout->Write();
  fout->Close();
  return;
}


//________________________________________________________________________________
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */) {
   return x == 0.0 && y == 0.0 ? 0.0 : 180 * degree + TMath::ATan2(-y, -x) * radian;
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
