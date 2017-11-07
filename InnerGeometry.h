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

void InnerGeometry(TString Title,Double_t *PropertyAllSiPM,Bool_t *MeasuredSiPM,Bool_t *validSiPM,Double_t Min,Double_t Max){
  TPaveText *pttitle = new TPaveText(.2,.925,.8,.975);
  pttitle->AddText(Title);
  pttitle->Draw();
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
  Int_t Split=4;
  Double_t Num;
  TText* NumText;
  for (int i = 0; i < Split+1; i++) {
    Num=Min+i*(Max-Min)/Split;
    Double_t Position=0.085+i*CBHeight/Split;
    NumText=new TText(0.80,Position,Form("%.2lf",Num));
    NumText->SetTextSize(0.03);
    NumText->Draw();
  }
}
