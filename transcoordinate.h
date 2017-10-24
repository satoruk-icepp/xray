Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */) {
  return x == 0.0 && y == 0.0 ? 0.0 : 180 + TMath::ATan2(-y, -x) *180/ TMath::Pi();
}
