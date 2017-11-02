#define NRow 93
#define NLine 44

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

void InnerGeometry(Double_t *PropertyAllSiPM,Bool_t *MeasuredSiPM,Bool_t *validSiPM,Double_t Min,Double_t Max){
  TArrow* Zaxisar= new TArrow(0.7,0.05,0.6,0.05);
  TText* Zaxisdesc= new TText(0.7,0.025,"Z");
  Zaxisar->Draw();
  Zaxisdesc->Draw();

  TArrow* Phiaxisar= new TArrow(0.75,0.2,0.75,0.3);
  TText* Phiaxisdesc= new TText(0.7,0.15,"phi");
  Phiaxisar->Draw();
  Phiaxisdesc->Draw();


  TText* UStext= new TText(0.7,0.5,"US");
  TText* DStext= new TText(0,0.5,"DS");
  UStext->Draw();
  DStext->Draw();


  Double_t CenterBoxWidth=0.6;
  Double_t CenterBoxHeight=0.8;
  Double_t CenterBoxOriginX=0.7;
  Double_t CenterBoxOriginY=(1+CenterBoxHeight)/2;
  Double_t SiPMBoxWidth=CenterBoxWidth/NLine;
  Double_t SiPMBoxHeight=CenterBoxHeight/NRow;
  for(int i = 0; i<NRow; i++){
	for(int j = 0; j<NLine; j++){
	  Int_t ChNum=i*NLine+j;
	  Double_t OneBoxXmax=CenterBoxOriginX-j*SiPMBoxWidth;
	  Double_t OneBoxYmax=CenterBoxOriginY-i*SiPMBoxHeight;
	  Double_t OneBoxXmin=OneBoxXmax-SiPMBoxWidth;
	  Double_t OneBoxYmin=OneBoxYmax-SiPMBoxHeight;
	  TBox *b= new TBox(OneBoxXmin,OneBoxYmin,OneBoxXmax,OneBoxYmax);
	  //TString chnum;
	  //chnum.Form("%d",i*NRow+j);
	  //TText *t=new TText(OneBoxXmin,OneBoxYmin,chnum);
	  Int_t Color=1;
	  Int_t FillStyle=1001;
	  if(validSiPM[ChNum]==true){
		Color=colordtm(PropertyAllSiPM[ChNum],Min,Max);
	  }
	  if(MeasuredSiPM[ChNum]==false){
		Color=15;
		FillStyle=3244;
	  }
	  b->SetFillColor(Color);
	  b->Draw();
	  // t->Draw();
	}
  }

  Double_t CBHeight=0.8;
  Double_t CBWidth=0.05;
  Double_t CBOriginX=0.85;
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
  TText *MinText=new TText(0.80,0.05,Form("%.2lf",Min));
  TText *MaxText=new TText(0.80,0.9,Form("%.2lf",Max));
  MinText->SetTextSize(0.03);
  MaxText->SetTextSize(0.03);
  MinText->Draw();
  MaxText->Draw();
}


void InnerGeometryArrow(Double_t *PropertyAllSiPM,Bool_t *valid ,Double_t Min,Double_t Max){
  TArrow* axisar= new TArrow(0.05,0.05,0.1,0.05);
  TText* axisdesc= new TText(0,0,"Z Axis");
  axisar->Draw();
  axisdesc->Draw();
  
  Double_t CenterBoxWidth=0.8;
  Double_t CenterBoxHeight=0.8;
  Double_t CenterBoxOriginX=0.85;
  Double_t CenterBoxOriginY=(1+CenterBoxHeight)/2;
  Double_t SiPMBoxWidth=CenterBoxWidth/NLine;
  Double_t SiPMBoxHeight=CenterBoxHeight/NRow;
  for(int i = 0; i<NRow; i++){
	for(int j = 0; j<NLine; j++){
	  if(valid[i*NLine+j]==true){
		Double_t OneBoxXmax=CenterBoxOriginX-j*SiPMBoxWidth;
		Double_t OneBoxYmax=CenterBoxOriginY-i*SiPMBoxHeight;
		Double_t OneBoxXmin=OneBoxXmax-SiPMBoxWidth;
		Double_t OneBoxYmin=OneBoxYmax-SiPMBoxHeight;
		//TBox *b= new TBox(OneBoxXmin,OneBoxYmin,OneBoxXmax,OneBoxYmax);
		Double_t OneBoxXcenter=(OneBoxXmax+OneBoxXmin)/2;
		Double_t OneBoxYcenter=(OneBoxYmax+OneBoxYmin)/2;
		Double_t ratio = std::abs(PropertyAllSiPM[i*NLine+j]-Min)/std::abs(Max-Min);
		TArrow* ar = new TArrow(OneBoxXcenter,OneBoxYcenter,OneBoxXcenter+ratio*SiPMBoxWidth/2,OneBoxYcenter);
		ar->SetArrowSize(0.01);
		ar->Draw();
	  }
	  //TString chnum;
	  //chnum.Form("%d",i*NRow+j);
	  //TText *t=new TText(OneBoxXmin,OneBoxYmin,chnum);
	  /* Int_t Color=colordtm(PropertyAllSiPM[i*NLine+j],Min,Max); */
	  /* b->SetFillColor(Color); */
	  /* b->Draw(); */
	  // t->Draw();
	}
  }
}



