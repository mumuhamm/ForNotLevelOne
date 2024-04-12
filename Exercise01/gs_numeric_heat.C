#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "TKey.h"
#include "RooExponential.h"
#include <map>
#include "TCut.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "TVirtualPad.h"
#include "RooDataHist.h"
#include <string>
#include "TEventList.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "TH2D.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooRealConstant.h"
#include "RooConstVar.h"
#include "Roo1DTable.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooRealSumPdf.h"
#include "Math/SpecFunc.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "Riostream.h"
#include "RooRandom.h"
#include "TMath.h"
#include "RooMCStudy.h"
#include "RooArgSet.h"
#include "RooLegendre.h"
#include "RooSpHarmonic.h"
#include "RooBifurGauss.h"
#include "complex.h"
#include <iostream>

using namespace RooFit; 
using namespace std;

const int M = 100;
const float dx =0.1;
const int kmax = 20;
const float a = 0.3;
class data{
    public:
    float kappa , theta0, thata1, thata2;
    float alpha[M]; //  alpha by hand = x*x +1 , x = a + (m-1)h
    void read();
}; 
void data::read(){
    kappa = 12.3;
    theta0 = 0.5;
    thata1 = 0.5;
    thata2 = 240;
    for(int i = 0; i < M; i++){
        alpha[i] = a + (i-1)*dx;
    }
}

void init (float *theta , ::data & d){
    for(int i = 0; i < M; i++){
        theta[i] =  d.thata1*(M-1-i)/(M-1.)+ d.thata2*i/(M-1.);
    }
}

void relax (float *theta, ::data & d){
   const float kappadx2 = d.kappa/dx/dx;
   for(int i =1; i< M-1; i++){
       theta[i] = (kappadx2*(theta[i+1] + theta[i-1]) + d.alpha[i]*d.theta0)/(kappadx2*2 + d.alpha[i]);
   }
}
void plot(float *theta ){
    TGraph *graph = new TGraph();
    for(int i =1 ; i < M-1; i++){
        graph->SetPoint(i, i, theta[i]);
    }
    TCanvas *canvas = new TCanvas("canvas", "Graph", 800, 600);
    graph->Draw("AL");  // "AL" option draws lines connecting the points
    graph->GetXaxis()->SetTitle("X Axis Title");
    graph->GetYaxis()->SetTitle("Y Axis Title");
    graph->SetTitle("Graph Title");
    canvas->Update();
    canvas->Draw();
}
void gs_numeric_heat(){

   ::data d;
    d.read();
    float theta[M];
    init(theta, d);
    for(int k = 0; k < kmax; k++)
        relax(theta, d);
    plot(theta);

   /*
    float kappa_values[] = {0.4, 13.5, 34.2, 53.2, 67.2};
    int num_kappa = sizeof(kappa_values) / sizeof(kappa_values[0]);

    TGraph *graphs[num_kappa];
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Define legend position

    for (int j = 0; j < num_kappa; j++) {
        ::data d;
        d.kappa = kappa_values[j];
        d.read();
        float theta[M];
        init(theta, d);
        relax(theta, d);
        graphs[j] = new TGraph(M);
        for (int i = 0; i < M; i++)graphs[j]->SetPoint(i, i, theta[i]);
        legend->AddEntry(graphs[j], Form("kappa=%.2f", kappa_values[j]), "l");
    }

    // Create a canvas to draw the graph
    TCanvas *canvas = new TCanvas("canvas", "Graph", 800, 600);

    // Draw each TGraph on the same canvas
    for (int j = 0; j < num_kappa; j++) {
        graphs[j]->SetLineColor(j + 1); // Set line color
        if (j == 0)
            graphs[j]->Draw("AL"); // Draw first graph with axis
        else
            graphs[j]->Draw("L SAME"); // Draw subsequent graphs on top
    }

    // Add legend to the canvas
    legend->Draw();

    // Update the canvas
    canvas->Update();
*/
}
