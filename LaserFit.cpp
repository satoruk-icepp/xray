#include <iostream>
#include <string>
#include <fstream>

using namespace std;

Double_t Rad2Deg(Double_t Rad);
Double_t PlaneEq(Double_t *x,Double_t *par);

void LaserFit(){
  string datapath="/meg/home/kobayashi_s/meg2/analyzer/x-ray/data/";
  string USFileName = "LaserDataUS.txt";
  string DSFileName = "LaserDataDS.txt";
  string USFilePath = datapath+USFileName;
  string DSFilePath = datapath+DSFileName;
  std::ifstream ifUS(USFilePath);
  string strUS;
  if(ifUS.fail()) {
	cerr << "File do not exist.\n";
	exit(0);
  }
  TGraph2D *grUS = new TGraph2D();
  TGraph* ProjUS=new TGraph();
  while(getline(ifUS, strUS)) {
	Float_t XPos;
	Float_t YPos;
	Float_t ZPos;
	sscanf(strUS.data(), "%f\t%f\t%f", &XPos, &YPos, &ZPos);
	grUS->SetPoint(grUS->GetN(),XPos,YPos,ZPos);
	ProjUS->SetPoint(ProjUS->GetN(),XPos,YPos);

  }

  std::ifstream ifDS(DSFilePath);
  string strDS;
  if(ifDS.fail()) {
	cerr << "File do not exist.\n";
	exit(0);
  }
  TGraph2D *grDS = new TGraph2D();
  while(getline(ifDS, strDS)) {
	Float_t XPos;
	Float_t YPos;
	Float_t ZPos;
	sscanf(strDS.data(), "%f\t%f\t%f", &XPos, &YPos, &ZPos);
	grDS->SetPoint(grDS->GetN(),XPos,YPos,ZPos);

  }
  gStyle->SetTitleOffset( 2,"XYZ");
  gStyle->SetTitleSize( 0.03,"XYZ");
  gStyle->SetLabelSize( 0.03,"XYZ");
  TCanvas* canvas1 = new TCanvas("canvas1", "Upstream",600,600);
  TCanvas* canvas2 = new TCanvas("canvas2", "Downstream",600,600);
  TCanvas* canvas3 = new TCanvas("canvas3", "Overall",600,600);
  TCanvas* canvas4 = new TCanvas("canvas4", "Projection",600,600);


  TF2 *plainUS = new TF2("planeUS",PlaneEq,-1500,1500,-1500,1500,3);
  TF2 *plainDS = new TF2("planeDS", PlaneEq,-1500,1500,-1500,1500,3);
  plainUS->SetParameters(0,0,-600);
  plainDS->SetParameters(0,0,600);

  //Up stream
  canvas1->cd(); 
  grUS->Fit("planeUS","N");
  grUS->SetMaximum(1500);
  grUS->SetMinimum(-1500);
  grUS->GetXaxis()->SetLimits(-1500,1500);
  grUS->GetYaxis()->SetLimits(-1500,1500);
  grUS->SetMarkerStyle(22);
  grUS->SetMarkerColor(kBlue);
  TString UpTitle="UpStream side face;";
  TString DownTitle="DownStream side face;";
  TString AllTitle="Both faces;";
  TString AxisTitle="X Position[mm];Y Position[mm];Z Position[m]";
  plainUS->SetTitle(UpTitle+AxisTitle);
  plainUS->Draw("surf");
  grUS->Draw("same p");
  Double_t ThetaUS = plainUS->GetParameter(0);
  Double_t PhiUS = plainUS->GetParameter(1);
  Double_t DegThetaUS = Rad2Deg(ThetaUS);
  Double_t DegPhiUS = Rad2Deg(PhiUS);
  std::cout<< "Tilt Angle Theta(about Y axis):  "<<DegThetaUS<<std::endl;
  std::cout<< "Rotation Angle Phi(about Z axis):  "<<DegPhiUS<<std::endl;
  //Down Stream
  canvas2->cd();
  grDS->Fit("planeDS","N");
  grDS->SetMaximum(1500);
  grDS->SetMinimum(-1500);
  grDS->GetXaxis()->SetLimits(-1500,1500);
  grDS->GetYaxis()->SetLimits(-1500,1500);
  grDS->SetMarkerStyle(22);
  grDS->SetMarkerColor(kRed);
  plainDS->SetTitle(DownTitle+AxisTitle);

  plainDS->Draw("surf");
  grDS->Draw("same p");
  Double_t ThetaDS = plainDS->GetParameter(0);
  Double_t PhiDS = plainDS->GetParameter(1);
  Double_t DegThetaDS = Rad2Deg(ThetaDS);
  Double_t DegPhiDS = Rad2Deg(PhiDS);
  std::cout<< "Tilt Angle Theta(about Y axis):  "<<DegThetaDS<<std::endl;
  std::cout<< "Rotation Angle Phi(about Z axis):  "<<DegPhiDS<<std::endl;

  canvas3->cd();
  grUS->SetTitle(AllTitle+AxisTitle);
  plainDS->SetRange (-1500,-1500,-1500,1500,1500,1500);
  plainDS->Draw("surf");

  grUS->Draw("p same");
  grDS->Draw("p same");
  Double_t planepar[3];
  for(int i=0;i<3;i++){
	planepar[i]=	plainUS->GetParameter(i);
  }

  Double_t xy[2]={0,0};

  std::cout<<PlaneEq(xy,planepar)<<std::endl;


  canvas4->cd();
  ProjUS->Draw("ap");

  /*
	auto chi2Function = [&](const Double_t *par, Double_t planepar[3]) {
	//minimisation function computing the sum of squares of residuals
	// looping at the graph points
	Int_t np = grUS->GetN();
	Double_t f = 0;
	Double_t *x = grUS->GetX();
	Double_t *y = grUS->GetY();
	Double_t *z = grUS->GetZ();

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
