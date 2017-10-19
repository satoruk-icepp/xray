#include <iostream>
#include <string>
#include <fstream>

using namespace std;

Double_t Rad2Deg(Double_t Rad);

void LaserFit(){
  string datapath="/meg/home/kobayashi_s/meg2/analyzer/x-ray/data/";
  string USFileName=datapath+ "LaserDataUS.txt";
  std::ifstream ifUS(USFileName);
  string str;
  if(ifUS.fail()) {
	cerr << "File do not exist.\n";
	exit(0);
  }
  TGraph2D *grUS = new TGraph2D();
  while(getline(ifUS, str)) {
	Float_t XPos;
	Float_t YPos;
	Float_t ZPos;
	sscanf(str.data(), "%f\t%f\t%f", &XPos, &YPos, &ZPos);
	grUS->SetPoint(grUS->GetN(),XPos,YPos,ZPos);

  }
  TCanvas* canvas1 = new TCanvas("canvas1", "Upstream",600,600);
  TCanvas* canvas2 = new TCanvas("canvas2", "Downstream",600,600);
 
  TGraph2D *grDS = new TGraph2D("$MEG2SYS/analyzer/x-ray/data/LaserDataDS.txt");

  TString PlaneName= "-tan([0])*(sin([1])*x-cos([1])*y)+[2]";

  TF2 *plainUS = new TF2("planeUS",PlaneName,-1500,1500,-1500,1500);
  plainUS->SetParameters(0,0,-600);
  TF2 *plainDS = new TF2("planeDS", PlaneName,-1500,1500,-1500,1500);
  plainDS->SetParameters(0,0,600);
  //TF1 *dummy = new TF1("dummy","x+cos(1)");
  //TCanvas* dum=new TCanvas("dummy","dumy",500,500);
  //dum->cd();
  //dummy->Draw();
  //Up stream
  canvas1->cd();
  grUS->Fit("planeUS","N");
  grUS->SetMaximum(1500);
  grUS->SetMinimum(-1500);
  grUS->GetXaxis()->SetLimits(-1500,1500);
  grUS->GetYaxis()->SetLimits(-1500,1500);
  grUS->SetMarkerStyle(22);
  plainUS->SetTitle("UpStream side face;Z Position[mm];X Position[mm];Y Position[m]");

  plainUS->Draw("surf");
  grUS->Draw("same p0");
  Double_t ThetaUS = plainUS->GetParameter(0);
  Double_t PhiUS = plainUS->GetParameter(1);
  Double_t DegThetaUS = Rad2Deg(ThetaUS);
  Double_t DegPhiUS = Rad2Deg(PhiUS);
  std::cout<< "Tilt Angle Theta(about Y axis):  "<<DegThetaUS<<std::endl;
  std::cout<< "Rotation Angle Phi(about Z axis):  "<<DegPhiUS<<std::endl;
  //Down Stream
  canvas2->cd();
  grDS->GetXaxis()->SetLimits(-1500,1500);
  grDS->Fit("planeDS","N");
  grDS->SetMarkerStyle(22);
  grDS->SetTitle("DownStream side face;X Position[mm];Y Position[mm]");
  grDS->Draw("ap");
  //plainDS->Draw("same");
}

Double_t Rad2Deg(Double_t Rad){
Double_t degree =Rad*180/TMath::Pi();
 Int_t range = floor(degree/360);
  return 
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
