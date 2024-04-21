// \author  Muhammad Alibordi 

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;
using namespace std; 


void DrawLabels(TCanvas* c){//, const TString& eraLabel) {
    TLatex* cmsLabel = new TLatex();
    cmsLabel->SetTextFont(42);
    cmsLabel->SetTextSize(0.04);
    cmsLabel->SetTextAlign(11); // Left-align
    cmsLabel->DrawLatexNDC(0.17, 0.92, "#bf{FUW} #it{Private}");

    TLatex* lumiLabel = new TLatex();
    lumiLabel->SetTextFont(42);
    lumiLabel->SetTextSize(0.04);
    lumiLabel->SetTextAlign(31); // Right-align
    TString lumiText =  "";//32fb^{-1}"DrellYan";// now just era we have lumi info though eraLabel;
    lumiLabel->DrawLatexNDC(0.94444, 0.92, lumiText);

    c->Update();
}
void toy2D_DecayTime_Mass() {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    Int_t nbins = 100;
    RooRealVar *svmass = new RooRealVar("svmass", "M_{B_{s}} GeV/c^{2}", 5.25, 5.49);
    RooRealVar *BsCt2DMC = new RooRealVar("BsCt2DMC", "c#tau_{B_{s}}", 0.007, 0.3, "cm");
    RooRealVar *BsCt2DMCErr = new RooRealVar("BsCt2DMCErr", "#Delta c#tau_{B_{s}}", 0.0002, 0.005, "cm");
    RooRealVar *tau = new RooRealVar("tau", "tau", 0.0441, 0.03, 0.055, "ps");

    RooRealVar *g1 = new RooRealVar("g1", "g1", 13.638, 1.0, 20.0);
    RooRealVar *g2 = new RooRealVar("g2", "g2", 6.082, 1.0, 15.00);
    RooRealVar *b1 = new RooRealVar("b1", "b1", 0.00009512, 0, 0.001);
    RooRealVar *b2 = new RooRealVar("b2", "b2", 0.0002825, 0, 0.001);
    RooRealVar *m1 = new RooRealVar("m1", "m1", 0.00019);
    RooRealVar *m2 = new RooRealVar("m2", "m2", 0.00019);
    RooGamma *G1 = new RooGamma("G1", "G1", *BsCt2DMCErr, *g1, *b1, *m1);
    RooGamma *G2 = new RooGamma("G2", "G2", *BsCt2DMCErr, *g2, *b2, *m2);
    RooRealVar *frac = new RooRealVar("frac", "frac", 0.5, 0, 1.0);
    RooRealVar norm("norm", "Normalization", 1.0);
    RooArgList pdfList("pdfList");
    pdfList.add(*G1);
    pdfList.add(*G2);
    RooArgList coefList("coefList");
    coefList.add(*frac);
    RooAddPdf *Di_Gamma = new RooAddPdf("Di_Gamma", "Di_Gamma", pdfList, coefList);
    BsCt2DMCErr->setRange("norm", BsCt2DMCErr->getMin(), BsCt2DMCErr->getMax());
    Di_Gamma->fixCoefRange("norm");

    RooRealVar *ctp0 = new RooRealVar("ctp0", "ctp0", -0.287321, -0.287321 - 0.0199639, -0.287321 + 0.0199639);
    RooRealVar *ctp1 = new RooRealVar("ctp1", "ctp1", -10.5223, -10.5223 - 0.0935714, -10.5223 + 0.0935714);
    RooRealVar *ctp2 = new RooRealVar("ctp2", "ctp2", 1.0, 1.0, 1.0); // Fixed value
    RooRealVar *ctp3 = new RooRealVar("ctp3", "ctp3", 13.8522, 13.8522 - 0.0517908, 13.8522 + 0.0517908);
    RooRealVar *ctp4 = new RooRealVar("ctp4", "ctp4", 1.03476, 1.03476 - 0.00213076, 1.03476 + 0.00213076);
    RooRealVar *ctp5 = new RooRealVar("ctp5", "ctp5", 4.29529, 4.29529 - 0.0169976, 4.29529 + 0.0169976);
    RooRealVar *ctp6 = new RooRealVar("ctp6", "ctp6", 0.118151, 0.118151 - 0.00212799, 0.118151 + 0.00212799);

    RooGenericPdf *cteffFunc = new RooGenericPdf("cteffFunc", "exp(ctp0+ctp1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,ctp2,ctp3,ctp4,ctp5,ctp6)",
        RooArgList(*ctp0, *ctp1, *ctp2, *ctp3, *ctp4, *ctp5, *ctp6, *BsCt2DMC));

    double kappa_val = 1.25029;
    RooRealVar *kappa = new RooRealVar("kappa", "kappa", kappa_val);
    RooFormulaVar *kapreso = new RooFormulaVar("kapreso", "kappa*BsCt2DMCErr", RooArgList(*kappa, *BsCt2DMCErr));
    RooRealVar *bias = new RooRealVar("bias", "bias", 0);
    RooRealVar *sigma = new RooRealVar("sigma", "per-event error scale factor", 1);
    RooGaussModel *gm = new RooGaussModel("gm", "gauss model scaled bt per-event error", *BsCt2DMC, *bias, *sigma, *kapreso);
    RooDecay *decay_gm = new RooDecay("decay_gm", "decay", *BsCt2DMC, *tau, *gm, RooDecay::SingleSided);
    RooEffProd *modelEff = new RooEffProd("modelEffbfrconv", "model with efficiency", *decay_gm, *cteffFunc);

    RooRealVar *mean = new RooRealVar("mean", "mean", 5.37, 5.25, 5.49);
    RooRealVar *sigma1 = new RooRealVar("sigma1", "sigma1", 0.03, 0., 0.5);
    RooRealVar *sigma2 = new RooRealVar("sigma2", "sigma2", 0.018, 0., 0.7);
    RooRealVar *sigma3 = new RooRealVar("sigma3", "sigma3", 0.022, 0., 0.9);
    RooGaussian *gauss1 = new RooGaussian("gauss1", "gauss1", *svmass, *mean, *sigma1);
    RooGaussian *gauss2 = new RooGaussian("gauss2", "gauss2", *svmass, *mean, *sigma2);
    RooGaussian *gauss3 = new RooGaussian("gauss3", "gauss3", *svmass, *mean, *sigma3);
    RooRealVar *frac1 = new RooRealVar("frac1", "frac1", 0.1681, 0.1, 0.4);
    RooRealVar *frac2 = new RooRealVar("frac2", "frac2", 0.45, 0.35, .60);
    RooAddPdf *T_gaus = new RooAddPdf("T_gaus", "T_gaus", RooArgList(*gauss1, *gauss2, *gauss3), RooArgList(*frac1, *frac2));

    RooProdPdf *sigmpdf = new RooProdPdf("sigmpdf", "mass*ctau*ctauerr", RooArgList(*T_gaus, *modelEff, *Di_Gamma));
    RooDataSet *data = sigmpdf->generate(RooArgList(*svmass, *BsCt2DMC, *BsCt2DMCErr), 100000);
    auto mbs_h = (TH1F *)data->createHistogram("mbs_h", *svmass);
    auto ctau_h = (TH1F *)data->createHistogram("ctau_h", *BsCt2DMC);
    auto ctauerr_h = (TH1F *)data->createHistogram("ctauerr_h", *BsCt2DMCErr);

    RooPlot *mbs_rp = svmass->frame(Title("m_{B_{s}} (GeV)"), Bins(nbins));
    data->plotOn(mbs_rp, DataError(RooAbsData::SumW2));
    RooPlot *ctau_rp = BsCt2DMC->frame(Title("B_{s} c#tau (cm)"), Bins(nbins));
    data->plotOn(ctau_rp, DataError(RooAbsData::SumW2));
    RooPlot *ctauerr_rp = BsCt2DMCErr->frame(Title("B_{s} c#tau error (cm)"), Bins(nbins));
    data->plotOn(ctauerr_rp, DataError(RooAbsData::SumW2));

    TCanvas *histCanvas = new TCanvas();
    histCanvas->Divide(3, 1);
    histCanvas->cd(1); 
    ctau_rp->Draw("HIST");
    DrawLabels(histCanvas);
    histCanvas->cd(2);
    mbs_rp->Draw("HIST");
    DrawLabels(histCanvas);
    histCanvas->cd(3);
    ctauerr_rp->Draw("HIST");
    DrawLabels(histCanvas);
    histCanvas->SaveAs("toy2D_DecayTime_Mass.png");
   
    RooRealVar *nSig_new = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC", 1400, 0., 100000);

    // Fit, instead copying the model, we use a new model with the same parameters

    RooRealVar *kappa_new = new RooRealVar("kappa_new", "kappa", kappa_val);
    RooFormulaVar *kapreso_new = new RooFormulaVar("kapreso_new", "kappa_new*BsCt2DMCErr", RooArgList(*kappa_new, *BsCt2DMCErr));
    RooRealVar *bias_new = new RooRealVar("bias_new", "bias", 0);
    RooRealVar *sigma_new = new RooRealVar("sigma_new", "per-event error scale factor", 1);
    RooGaussModel *gm_new = new RooGaussModel("gm_new", "gauss model scaled bt per-event error", *BsCt2DMC, *bias_new, *sigma_new, *kapreso_new);
    RooDecay *decay_gm_new = new RooDecay("decay_gm_new", "decay", *BsCt2DMC, *tau, *gm_new, RooDecay::SingleSided);
    RooEffProd *modelEff_new = new RooEffProd("modelEffbfrconv_new", "model with efficiency", *decay_gm, *cteffFunc);

    RooRealVar *mean_new = new RooRealVar("mean_new", "mean", 5.37, 5.25, 5.49);
    RooRealVar *sigma1_new = new RooRealVar("sigma1_new", "sigma1", 0.03, 0., 0.5);
    RooRealVar *sigma2_new = new RooRealVar("sigma2_new", "sigma2", 0.018, 0., 0.7);
    RooRealVar *sigma3_new = new RooRealVar("sigma3_new", "sigma3", 0.022, 0., 0.9);
    RooGaussian *gauss1_new = new RooGaussian("gauss1_new", "gauss1", *svmass, *mean_new, *sigma1_new);
    RooGaussian *gauss2_new = new RooGaussian("gauss2_new", "gauss2", *svmass, *mean_new, *sigma2_new);
    RooGaussian *gauss3_new = new RooGaussian("gauss3_new", "gauss3", *svmass, *mean_new, *sigma3_new);
    RooRealVar *frac1_new = new RooRealVar("frac1_new", "frac1", 0.1681, 0.1, 0.4);
    RooRealVar *frac2_new = new RooRealVar("frac2_new", "frac2", 0.45, 0.35, .60);
    RooAddPdf *T_gaus_new = new RooAddPdf("T_gaus", "T_gaus", RooArgList(*gauss1_new, *gauss2_new, *gauss3_new), RooArgList(*frac1_new, *frac2_new));

    RooProdPdf *sigmpdf_new = new RooProdPdf("sigmpdf", "mass*mistag", RooArgList(*T_gaus_new, *decay_gm_new));
    RooAddPdf *MTpdf_new = new RooAddPdf("MTpdf", "Total pdf", RooArgList(*sigmpdf_new), RooArgList(*nSig_new));

    RooFitResult *TwoDFitRes = MTpdf_new->fitTo(*data, Save());
    //RooFitResult* TwoDFitRes = MTpdf_new->fitTo(*data,Save(),Extended(1),ConditionalObservables(*BsCt2DMCErr));
    TwoDFitRes->Print("v");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    TH2 *hcorr = TwoDFitRes->correlationHist();
    TCanvas *cres = new TCanvas("Correlation Matrix", "Correlation Matrix", 1000, 600);
    hcorr->GetYaxis()->SetTitleOffset(1.4);
    hcorr->Draw("colz");
    DrawLabels(cres);

    data->plotOn(mbs_rp, DataError(RooAbsData::SumW2));
    MTpdf_new->plotOn(mbs_rp);
    MTpdf_new->paramOn(mbs_rp);
    MTpdf_new->plotOn(mbs_rp, LineColor(kBlue), LineWidth(2));
    Double_t chisquare_mass = mbs_rp->chiSquare();
    RooPlot *pullframe = svmass->frame(RooFit::Title("Mass pull"));
    RooHist *hpull1 = mbs_rp->pullHist();
    pullframe->addPlotable(hpull1, "P0");
    pullframe->SetMinimum(-3);
    pullframe->SetMaximum(+3);
    pullframe->SetYTitle("Mass-pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);

    MTpdf_new->plotOn(mbs_rp, Components(*gauss1_new), LineColor(3),RooFit::Name("signalone"), LineWidth(2), LineStyle(4));
    MTpdf_new->plotOn(mbs_rp, Components(*gauss2_new), LineColor(2),RooFit::Name("signaltwo"), LineWidth(2), LineStyle(2));
    MTpdf_new->plotOn(mbs_rp, Components(*gauss3_new), LineColor(6),RooFit::Name("signalthree"), LineWidth(2), LineStyle(2));
    
    TLegend *legm = new TLegend(0.7,0.7,0.9,0.9);
    legm->SetHeader("#bf{FUW} #it{Private}","C");
    legm->AddEntry(mbs_rp->findObject("signalone"),"Gauss(M_{B_{s}^{1}})","l");
    legm->AddEntry(mbs_rp->findObject("signaltwo"),"Gauss(M_{B_{s}^{2}})","l");
    legm->AddEntry(mbs_rp->findObject("signalthree"),"Gauss(M_{B_{s}^{3}})","l");
    

  

    data->plotOn(ctau_rp, DataError(RooAbsData::SumW2));
    MTpdf_new->plotOn(ctau_rp);
    MTpdf_new->paramOn(ctau_rp);
    MTpdf_new->plotOn(ctau_rp, LineColor(kBlue),RooFit::Name("signalct"), LineWidth(2), LineStyle(4));
    TLegend *legct = new TLegend(0.7,0.7,0.9,0.9);
    legct->SetHeader("#bf{FUW} #it{Private}","C");
    legct->AddEntry(ctau_rp->findObject("signalct"),"Expo(c#tau_{B_{s}^{1}})","l");


    RooPlot *pullframelt = BsCt2DMC->frame(100);
    RooHist *hpulllt = ctau_rp->pullHist();
    pullframelt->addPlotable(hpulllt, "P0");
    pullframelt->SetMinimum(-3);
    pullframelt->SetMaximum(+3);
    pullframelt->SetYTitle("ct-pull");
    pullframelt->SetMarkerStyle(20);
    pullframelt->SetNdivisions(10);

    TCanvas *c = new TCanvas("c", "c", 0, 0, 600, 600);
    c->SetFillColor(0);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05); // Set right margin to 5% of canvas width
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1);
    pad1->SetTopMargin(0.2); // Set top margin to 1% for title
    pad1->SetBottomMargin(0.1); // Set bottom margin for axis labels
    pad1->SetBorderMode(0); // Remove pad borders
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
    pad2->SetBottomMargin(0.1); // Set bottom margin for axis labels
    pad2->SetBorderMode(0); // Remove pad borders
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0); 
    mbs_rp->SetStats(0); 
    mbs_rp->Draw();
    legm->Draw();
    pad2->cd();
    pullframe->SetStats(0); 
    pullframe->Draw();
    DrawLabels(c);
    c->Update();

    TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 800, 600);
    TPad *pad11 = new TPad("pad11", "pad11", 0, 0.35, 1, 1);
    pad11->SetBottomMargin(0.00001);
    pad11->SetTopMargin(0.2); // Set top margin to 1% for title
    pad11->SetBottomMargin(0.1); // Set bottom margin for axis labels
    pad11->SetBorderMode(0); // Remove pad borders
    TPad *pad22 = new TPad("pad22", "pad22", 0, 0, 1, 0.33);
    pad22->SetTopMargin(0.00001);
    pad22->SetBottomMargin(0.1);
    pad22->SetBorderMode(0);
    pad11->Draw();
    pad22->Draw();
    pad11->cd();
    gStyle->SetOptTitle(0);
    cc->SetFillColor(0); 
    cc->SetBorderSize(2); 
    cc->SetLeftMargin(0.14); 
    cc->SetRightMargin(0.05); 
    ctau_rp->SetStats(0); 
    ctau_rp->Draw();
    legct->Draw();
    pad22->cd();
    pullframelt->Draw();
    cc->Update();

}
