Int_t arraysearch(std::vector<Double_t> array, Double_t value){
  Bool_t found=false;
  Int_t sizearr=array.size();
  int i=-1;
  if(sizearr!=0){
    for (i = 0; i < sizearr; i++) {
      if (array[i]==value) {
	found=true;
	break;
      }
    }
  }
  if(found==false){
    i=-1;  
  }
  return i;
}
