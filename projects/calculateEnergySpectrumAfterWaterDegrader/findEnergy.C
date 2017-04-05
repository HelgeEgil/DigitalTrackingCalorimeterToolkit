#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

void Run(int mm) {
   TFile *f = new TFile(Form("energyspec_%dmm.root", mm));
   TH1F *h = (TH1F*) f->Get("energySpectrum");

   Float_t expectedEnergy = 255.3  - 0.58 * mm + 0.001387 * pow(mm,2) - 3.96e-6 * pow(mm,3);

   h->SetBinContent(1, 0);
   h->SetBinContent(2, 0);
   h->Draw();
   TF1 *fit = new TF1("fit", "gaus");
//   fit->SetParameter(0, h->GetMaximum());
   fit->SetParameter(1, expectedEnergy);
   fit->SetParameter(2, 2);
   h->Fit(fit, "B,Q");

   cout << mm << " " << fit->GetParameter(1) << " " << fit->GetParameter(2) << endl;
}
