Double_t PhiQCut[5]={0.05,0.2,0.2,8,0.4};
Double_t ZQCut[5]={0.5,2,2,100,0.5};

Bool_t DataQual(Double_t *FitErr, Bool_t PhiScan, Double_t Position, Double_t Design){
  Bool_t Quality=true;
  Double_t Gap=std::abs(Position-Design);
  Double_t PhiGapMax=1;
  Double_t ZGapMax=20;
  if(PhiScan==true){
	for(int i=0;i<5;i++){
	  if(FitErr[i]>PhiQCut[i]){
		Quality=false;
		break;
	  }
	  if(Gap>PhiGapMax){
		Quality=false;
	  }
	}
  }else{
	for(int i=0;i<5;i++){
	  if(FitErr[i]>ZQCut[i]){
		Quality=false;
		break;
	  }
	}
	if(Gap>ZGapMax){
	  Quality=false;	  
	}
  }
  return Quality;
}

