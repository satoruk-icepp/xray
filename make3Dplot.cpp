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

Double_t ZPosXray[nMPPC];
Double_t ZPosGapAllch[nMPPC];
Double_t PhiPosGapAllch[nMPPC];


void InnerGeometry(Double_t PropertyAllSiPM[nMPPC],Double_t Min, Double_t Max);
int colordtm(Double_t value,Double_t vmin, Double_t vmax);

void make3Dplot(){
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  gStyle->SetPalette(1);
  TCanvas* canvas1=new TCanvas("canvas1","Phi gap",600,600);
  //TCanvas* canvas2=new TCanvas("canvas2","Z gap",600,600);
  //TCanvas* canvas3=new TCanvas("canvas3","Z position",600,600);
  //TCanvas* canvas4=new TCanvas("canvas4","Gap Correlation",600,600);
  //TCanvas* canvas5 = new TCanvas("canvas5","neighbor",600,600);

  TGraph2D* geometry =new TGraph2D();

  TString farodatapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  TString farofilename = "XECFAROMPPC.root";
  TString farorootfile = farodatapath+farofilename;
  TFile* ffaro = new TFile(farorootfile.Data(),"read");
  TTree* tfaro = (TTree*)ffaro->Get("faro");

  /*Define CFRP*/
  Int_t CFRPOrigin[5]={0,24,47,70,93};
  Double_t CFRPGap[4]={0,0.6,2.1,3.7};

  /*faro tree*/
  Int_t FaroChNum;
  Double_t FaroXPos;
  Double_t FaroYPos;
  Double_t FaroZPos;
  Bool_t FaroDataQual;
  tfaro->SetBranchAddress("channel",&FaroChNum);
  tfaro->SetBranchAddress("XPos",&FaroXPos);
  tfaro->SetBranchAddress("YPos",&FaroYPos);
  tfaro->SetBranchAddress("ZPos",&FaroZPos);
  tfaro->SetBranchAddress("DataQual",&FaroDataQual);

  for(int iCh=0;iCh<nMPPC;iCh++){
	tfaro->GetEntry(iCh);
	if(FaroDataQual==true){
	}
	geometry->SetPoint(iCh,FaroZPos,FaroXPos,FaroYPos);

  }

  canvas1->cd();
  geometry->SetTitle("Faro Geometry;Z Position[mm];X Position[mm];Y Position[mm]");
  geometry->SetMaximum(800);
  geometry->SetMinimum(-800);
  geometry->GetXaxis()->SetLimits(-400,400);
  geometry->GetYaxis()->SetLimits(-800,800); 

  //geometry->SetMarkerSize(5);
  geometry->Draw("pcol");
  /*
  
	canvas2->cd();
	TPaveText *ptZ = new TPaveText(.2,.925,.8,.975);
	ptZ->AddText("Deviation between Designed Value in the Z Direction");
	ptZ->Draw();
	InnerGeometry(ZPosGapAllch,-10.0,10.0);
  
	canvas3->cd();
	grChZPos->SetMarkerStyle(22);
	grChZPos->SetMarkerColor(2);
	grChZPos->SetTitle("Deviation from the design value;Z_{Xray}[mm];Deviation[mm]");
	grChZPos->Draw("ap");
	canvas4->cd();
	grGapCor->GetXaxis()->SetLimits(-10,10);
	grGapCor->SetMinimum(-10);
	grGapCor->SetMaximum(10);
	grGapCor->Draw("ap");
  */
  /*
	canvas5->cd();
	canvas5->SetGrid(0,0);
	gStyle->SetFuncColor(kRed);
	WidthHist->Fit("gaus");
	WidthHist->Draw();*/
}

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

void InnerGeometry(Double_t PropertyAllSiPM[nMPPC],Double_t Min,Double_t Max){

  Double_t CenterBoxWidth=0.8;
  Double_t CenterBoxHeight=0.8;
  Double_t CenterBoxOriginX=0.85;
  Double_t CenterBoxOriginY=(1+CenterBoxHeight)/2;
  Double_t SiPMBoxWidth=CenterBoxWidth/NLine;
  Double_t SiPMBoxHeight=CenterBoxHeight/NRow;
  for(int i = 0; i<NRow; i++){
	for(int j = 0; j<NLine; j++){
	  Double_t OneBoxXmax=CenterBoxOriginX-j*SiPMBoxWidth;
	  Double_t OneBoxYmax=CenterBoxOriginY-i*SiPMBoxHeight;
	  Double_t OneBoxXmin=OneBoxXmax-SiPMBoxWidth;
	  Double_t OneBoxYmin=OneBoxYmax-SiPMBoxHeight;
	  TBox *b= new TBox(OneBoxXmin,OneBoxYmin,OneBoxXmax,OneBoxYmax);
	  //TString chnum;
	  //chnum.Form("%d",i*NRow+j);
	  //TText *t=new TText(OneBoxXmin,OneBoxYmin,chnum);
	  Int_t Color=colordtm(PropertyAllSiPM[i*NLine+j],Min,Max);
	  b->SetFillColor(Color);
	  b->Draw();
	  //t->Draw();
	}
  }

  Double_t CBHeight=0.8;
  Double_t CBWidth=0.05;
  Double_t CBOriginX=0.95;
  Double_t CBOriginY=0.1;
  Int_t Ncolor=100;
  Double_t OneCBHeight=CBHeight/Ncolor;
  for(int i=0;i<Ncolor;i++){
	TBox *b= new TBox(CBOriginX,CBOriginY+i*OneCBHeight,CBOriginX+CBWidth,CBOriginY+(i+1)*OneCBHeight);
	Int_t CBColor=colordtm(i,0,100);
	b->SetFillColor(CBColor);
	b->Draw();
  }
  TString MinNum,MaxNum;
  MinNum.Form("%f",Min);
  MaxNum.Form("%f",Max);
  TText *MinText=new TText(0.9,0.05,Form("%.2lf [mm/deg]",Min));
  TText *MaxText=new TText(0.9,0.9,Form("%.2lf [mm/deg]",Max));
  MinText->SetTextSize(0.02);
  MaxText->SetTextSize(0.02);
  MinText->Draw();
  MaxText->Draw();
}


  int colordtm(Double_t value,Double_t vmin, Double_t vmax){
	if (value<vmin) {
	  value=vmin;
	}
	if (value>vmax) {
	  value=vmax;
	}
	double ratio=std::abs(value-vmin)/std::abs(vmax-vmin);
	double cmax=100;
	double cmin=51;
	double crange=cmax-cmin;
	double dbcolor=cmin+crange*ratio;
	int color=floor(dbcolor);
	return color;
  }
