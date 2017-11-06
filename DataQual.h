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

Bool_t LargeZQual(Double_t ZPosition, Bool_t ZMeasured){
  Bool_t Quality=false;
  if(std::abs(ZPosition)<120&&ZMeasured==true){
    Quality= true;
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


Bool_t JudgeQual(Int_t *QualArray,Bool_t PhiScan,Bool_t ZMeasured,Bool_t PhiMeasured, Double_t ZPos,Double_t PhiPos,Double_t PosDesign ,Double_t *Err){
  Bool_t DataQual=false;
  Bool_t Measured=false;
  Double_t Pos;
  if(PhiScan==true&&PhiMeasured==true){
    Measured=true;
  }else if(PhiScan==false&&ZMeasured==true){
    Measured=true;
  }
  if(PhiScan==true){
    Pos=PhiPos;
  }else{
    Pos=ZPos;
  }
  if(Measured==true){
    ++QualArray[0];
    if(LargeZQual(ZPos,ZMeasured)==true){
      ++QualArray[1];
      if(FitQual(Err,PhiScan)==true){
        ++QualArray[2];
        if(DevQual(Pos,PosDesign,PhiScan)==true){
          ++QualArray[3];
          DataQual=true;
        }
      }
    }
  }
  return DataQual;
}
