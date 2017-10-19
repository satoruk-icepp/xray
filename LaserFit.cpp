
void LaserFit(){
  TCanvas* canvas1 = new TCanvas("canvas1", "Upstream",600,600);
  TCanvas* canvas2 = new TCanvas("canvas2", "Downstream",600,600);
  TGraph2D *grUS = new TGraph2D("$MEG2SYS/analyzer/x-ray/data/LaserDataUS.txt");
  TGraph2D *grDS = new TGraph2D("$MEG2SYS/analyzer/x-ray/data/LaserDataDS.txt");

  TF2 *plainUS = new TF2("planeUS", "-tan([0])*(sin([1])*x-cos([1])*y)+[2]",-1500,1500,-1500,1500);
  TF2 *plainDS = new TF2("planeDS", "-tan([0])*(sin([1])*x-cos([1])*y)+[2]",-1500,1500,-1500,1500);
  //Up stream
  canvas1->cd();
  grUS->Fit("planeUS","N");
  grUS->SetMaximum(1500);
  grUS->SetMinimum(-1500);
  grUS->SetMarkerStyle(22);
  grDS->SetTitle("UpStream side face;X Position[mm];Y Position[mm]");

  grUS->Draw("ap");
  //plainUS->Draw("same");
  //Down Stream
  canvas2->cd();
  grDS->GetXaxis()->SetLimits(-1500,1500);
  grDS->Fit("planeDS","N");
  grDS->SetMarkerStyle(22);
  grDS->SetTitle("DownStream side face;X Position[mm];Y Position[mm]");
  grDS->Draw("ap");
  //plainDS->Draw("same");
}

Double_t Plane(Double_t *x, Double_t *par){
  Double_t xx=x[0];
  Double_t yy=x[1];
  Double_t z = par[0]*xx+par[1]*yy+par[2];
  return z;

}

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
