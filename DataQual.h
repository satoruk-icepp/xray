Double_t PhiQCut[5]={0.05,0.2,0.2,8,0.4};
Double_t ZQCut[5]={0.5,2,2,100,0.5};

Bool_t FitQual(Double_t *FitErr, Bool_t PhiScan){
  Bool_t Quality=true;
  if(PhiScan==true){
	for(int i=0;i<5;i++){
	  if(FitErr[i]>PhiQCut[i]){
		Quality=false;
		break;
	  }
	}
  }else{
	for(int i=0;i<5;i++){
	  if(FitErr[i]>ZQCut[i]){
		Quality=false;
		break;
	  }
	}
  }
  return Quality;
}

Bool_t DevQual(Double_t Position, Double_t Design, Bool_t PhiScan){
  Bool_t Quality=true;
  Double_t Gap=std::abs(Position-Design);
  Double_t PhiGapMax=1;
  Double_t ZGapMax=20;
  if(PhiScan==true){
	if(Gap>PhiGapMax){
	  Quality=false;
	}
  }else{
	if(Gap>ZGapMax){
	  Quality=false;	  
	}
  }
  return Quality;
}

Bool_t Chi2Qual(Double_t ChiSq, Bool_t PhiScan){
  Bool_t Quality =true;
  Double_t ZChiSqMax = 5000;
  Double_t PhiChiSqMax = 500;
  if(PhiScan==true){
	if(ChiSq>PhiChiSqMax){
	  Quality=false;
	}
  }else{
	if(ChiSq>ZChiSqMax){
	  Quality = false;
	}
  }
  return Quality;
}
