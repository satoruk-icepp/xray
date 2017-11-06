#include <iostream>
#include <string>
#include <fstream>

using namespace std;

Double_t Rad2Deg(Double_t Rad);
Double_t PlaneEq(Double_t *x,Double_t *par);

void ReferencePlots(){
  string datapath="/meg/home/kobayashi_s/meg2/analyzer/macros/xec/survey/XECTransform/";
  string TrackerData = "XECLEICAreference.dat";
  string FAROData = "XECFAROreferenceOnCOBRA.dat";
  string TrackerPath = datapath+TrackerData;
  string FAROPath = datapath+FAROData;
  //  std::cout<<TrackerPath<<std::endl;
  std::ifstream ifLEICA(TrackerPath);
  string strLEICA;
  if(ifLEICA.fail()) {
    cerr << "File do not exist.\n";
    exit(0);
  }
  TGraph2D *grLEICA = new TGraph2D();
  //TGraph* ProjLEICA=new TGraph();

  std::vector<std::vector<Double_t>> LEICAPos;
  std::vector<std::vector<Double_t>> FAROPos;
  LEICAPos.resize(12);
  FAROPos.resize(12);
  Int_t Linc=0;
  while(getline(ifLEICA, strLEICA)) {
    Double_t Pos[3];
    sscanf(strLEICA.data(), "%lf\t%lf\t%lf", &Pos[0], &Pos[1], &Pos[2]);
    grLEICA->SetPoint(grLEICA->GetN(),Pos[0],Pos[1],Pos[2]);
    //ProjUS->SetPoint(ProjUS->GetN(),XPos,YPos);
//std::cout<<"debug"<<std::endl;
    LEICAPos[Linc].resize(3);
    for(int i=0;i<3;i++){
      LEICAPos[Linc][i]=Pos[i];
    }
    ++Linc;
  }

  std::ifstream ifFARO(FAROPath);
  string strFARO;
  if(ifFARO.fail()) {
    cerr << "File do not exist.\n";
    exit(0);
  }
  TGraph2D *grFARO = new TGraph2D();
  Int_t Finc=0;
  while(getline(ifFARO, strFARO)) {
    Double_t Pos[3];
    sscanf(strFARO.data(), "%lf\t%lf\t%lf", &Pos[0], &Pos[1], &Pos[2]);
    grFARO->SetPoint(grFARO->GetN(),Pos[0],Pos[1],Pos[2]);
    FAROPos[Finc].resize(3);
    for(int i=0;i<3;i++){
      FAROPos[Finc][i]=Pos[i];
    }
    ++Finc;

  }
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1 = new TCanvas("canvas1", "Upstream",600,600);
  TCanvas* canvas2 = new TCanvas("canvas2", "Downstream",600,600);
  TCanvas* canvas3 = new TCanvas("canvas3", "Overall",600,600);
  //TCanvas* canvas4 = new TCanvas("canvas4", "Projection",600,600);

for(int i=0;i<10;i++){
  std::cout<<"difference x:"<<LEICAPos[i][0]-FAROPos[i][0]<<" y: "<<LEICAPos[i][1]-FAROPos[i][1]<<" z: "<<LEICAPos[i][2]-FAROPos[i][2]<<std::endl;
}
  //TF2 *plainUS = new TF2("planeUS",PlaneEq,-1500,1500,-1500,1500,3);
  // TF2 *plainFARO = new TF2("planeFARO", PlaneEq,-1500,1500,-1500,1500,3);
  // plainUS->SetParameters(0,0,-600);
  // plainFARO->SetParameters(0,0,600);

  //Up stream
  canvas1->cd();
  //  grLEICA->Fit("planeUS","N");
  grLEICA->SetMaximum(-595);
  grLEICA->SetMinimum(-605);
  grLEICA->GetXaxis()->SetLimits(-1500,0);
  grLEICA->GetYaxis()->SetLimits(-1500,1500);
  grLEICA->SetMarkerStyle(20);
  grLEICA->SetMarkerSize(1.0);
  grLEICA->SetMarkerColor(kBlue);
  TString UpTitle="UpStream side face;";
  TString DownTitle="DownStream side face;";
  TString AllTitle="Faro and Leica;";
  TString AxisTitle="X Position[mm];Y Position[mm];Z Position[m]";

  grLEICA->Draw("ap");

  canvas2->cd();
  // grFARO->Fit("planeFARO","N");
  grFARO->SetMaximum(0);
  grFARO->SetMinimum(-1000);
  grFARO->GetXaxis()->SetLimits(-1500,0);
  grFARO->GetYaxis()->SetLimits(-1500,1500);
  grFARO->SetMarkerStyle(20);
  grFARO->SetMarkerSize(1.0);
  grFARO->SetMarkerColor(kRed);

  // plainFARO->SetTitle(DownTitle+AxisTitle);

  // plainFARO->Draw("surf");
  grFARO->Draw("ap");
  // Double_t ThetaFARO = plainFARO->GetParameter(0);
  // Double_t PhiFARO = plainFARO->GetParameter(1);
  // Double_t DegThetaFARO = Rad2Deg(ThetaFARO);
  // Double_t DegPhiFARO = Rad2Deg(PhiFARO);
  // std::cout<< "Tilt Angle Theta(about Y axis):  "<<DegThetaFARO<<std::endl;
  // std::cout<< "Rotation Angle Phi(about Z axis):  "<<DegPhiFARO<<std::endl;

  canvas3->cd();
  grLEICA->SetTitle(AllTitle+AxisTitle);
  // plainFARO->SetRange (-1500,-1500,-1500,1500,1500,1500);
  // plainFARO->Draw("surf");

  grLEICA->Draw("p0");
  grFARO->Draw("same p0");
  // Double_t planepar[3];
  // for(int i=0;i<3;i++){
  // 	planepar[i]=	plainUS->GetParameter(i);
  // }

  // Double_t xy[2]={0,0};

  // std::cout<<PlaneEq(xy,planepar)<<std::endl;


  // canvas4->cd();
  // ProjUS->Draw("ap");

  /*
     auto chi2Function = [&](const Double_t *par, Double_t planepar[3]) {
  //minimisation function computing the sum of squares of residuals
  // looping at the graph points
  Int_t np = grLEICA->GetN();
  Double_t f = 0;
  Double_t *x = grLEICA->GetX();
  Double_t *y = grLEICA->GetY();
  Double_t *z = grLEICA->GetZ();

  for (Int_t i=0;i<np;i++) {
  Double_t u = x[i] - par[0];
  Double_t v = y[i] - par[1];
  Double_t coor[2]={x[i],y[i]};
  Double_t z0 = PlaneEq(coor,*planepar);
  Double_t w = z[i] - z0;
  Double_t dr = par[2] - std::sqrt(u*u+v*v+w*w);
  f += dr*dr;
  }
  return f;
  };*/
}


Double_t PlaneEq(Double_t *x,Double_t *par){
  Double_t xx=x[0];
  Double_t yy=x[1];
  Double_t f = -tan(par[0])*(cos(par[1])*xx-sin(par[1])*yy)+par[2];
  return f;
}

Double_t Rad2Deg(Double_t Rad){
  Double_t degree =Rad*180/TMath::Pi();
  Int_t range = floor(degree/360);
  return degree -range*360;
}
/*
   Double_t Plane(Double_t *x, Double_t *par){
   Double_t xx=x[0];
   Double_t yy=x[1];
   Double_t z = par[0]*xx+par[1]*yy+par[2];
   return z;

   }*/

/*Double_t XECCylinder(Double_t *x,Double_t *par){
  Float_t xx=x[0];
  Float_t yy=x[1];
//Double_t f = par[0]/sqrt(2*TMath::Pi()*par[1])*TMath::Exp(-pow(xx-par[2],2)/(2*par[1]));
//par[0]:x coordinate of the center point
//par[1]:y coordinate of the center point
//par[2]:z coordinate of the center point
//par[3]:radius of the circle
Double_t OriginX = par[0];
Double_t OriginY = par[1];
Double_t OriginZ = par[2];
Double_t Theta = par[3];
Double_t Phi = par[4];
Double_t Radius = par[5];
Double_t z;
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



}*/
