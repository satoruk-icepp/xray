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
	  // t->Draw();
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


