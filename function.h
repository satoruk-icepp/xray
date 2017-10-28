Double_t TwoGaus(Double_t *x,Double_t *par){
  Double_t xx=x[0];
  Double_t mean=par[0];
  Double_t topwidth=par[1];
  Double_t meanfg=par[0]-par[1]/2;
  Double_t meansg=par[0]+par[1]/2;
  Double_t sigma=par[2];
  Double_t height=par[3];
  Double_t offset=par[4];
  Double_t f;
  if(xx<meanfg){
	f=offset+height*TMath::Gaus(xx,meanfg,sigma,true);
  }else if(xx>=meanfg&&xx<meansg){
	f=offset+ height*TMath::Gaus(meanfg,meanfg,sigma,true);
  }else{
	f=offset+ height*TMath::Gaus(xx,meansg,sigma,true);
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
  Double_t offset = par[4];
  Double_t f;
  if(xx<redge||xx>=fedge){
	f = offset;
  }else if(xx>=redge&&xx<topstart){
	Double_t ratio=std::abs((xx-redge)/(topstart-redge));
	f = offset+height*ratio;
  }else if(xx>=topstart&&xx<topend){
	f = offset+height;  
  }else{
	Double_t ratio=std::abs((xx-topend)/(fedge-topend));
	f = offset+height*(1-ratio);
  }
  return f;
}
                                  
Double_t TDLine(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Double_t yy = x[1];
  Double_t zz = par[0]*xx+par[1]*yy+par[2];
  return zz;
}
