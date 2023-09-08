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


void stepOne() {

//some initial declaration 
Int_t nbins = 100;
RooRealVar *massJpsi = new RooRealVar("massJpsi", "massJpsi", 3.095, 3.097); // charmonium mu+mu- decya mode , not exact but functional

/*This rooreal var is a constaructor , even though you do not instatiate the base class you system 
wide installation of the roofit sofware ,and source to the roofit library takes care of the instatiation 
of the base class , and they a on fly library load callled gSystem which automatically loads the library, 
this constructor takes four inputs first two are const char ( consider static string, but you could also 
ask for the char name somewhere in your code, for instance when a roorealvar of histogram is filled you 
could just ask myhisto->GetName())and then you could use the same histogram for further decation and fit
*/

//The model which will generate the mass histogram 
RooRealVar * mean = new RooRealVar("mean","m_{J/#psi}", 3.096, 3.0, 3.2);
RooRealVar * sigma = new RooRealVar("sigma","sigma_{J/#psi}", 0, 0.0003);
RooGaussian *gaus = new RooGaussian("gaus","gaus", *massJpsi, *mean,*sigma);

/* the Gaussian model takes five arguments , so all these base class referneces can be found in the root site
Using pointer does not make any difference in physics or math, its the C++ thing, for revising C++ things please 
visit cplusplus or programize. Now the model is ready we can use some root functinality to create some random data */

RooDataSet* data = gaus->generate(*massJpsi,100000);

/*It takes the mass of the particle , but see how does it code , you fill the rooreal var , when people talk , they talk in the 
language of physics, which makes more sense, so in this case someone would say : we could use the Jpsi mass to generate 
some random events, and experienced people are more inclined toward talking in the language of physics, so the idea is translating
all lines of code in the language of physics so, we used the Jpsi mass within the scope of its uncertainty, though it is not exact,
then a natural width is assigned, even that is not exact too, 0.3MeV is way more biased than what PDG provides , jpsi width is in KeV, 
but due to the misreconstruction of the tracks the natural resolution decreases, thats why the leverage. 
How many events : 100000
exercise : please tanslate into physics language these lines ( from one of my muonselection code) :
bool isTkIsolated = false;
if (im->isIsolationValid() && im->isolationR03().sumPt/im->pt() < theConfig.getParameter<double>("cutTkIsoRel") ) isTkIsolated = true;
bool isPFIsolated = false;
if (im->isPFIsolationValid() && 
           (   im->pfIsolationR04().sumChargedHadronPt 
             + std::max(0., im->pfIsolationR04().sumNeutralHadronEt + im->pfIsolationR04().sumPhotonEt - 0.5*im->pfIsolationR04().sumPUPt)
            )/im->pt() < theConfig.getParameter<double>("cutPFIsoRel")
        ) isPFIsolated = true;
this one theConfig.getParameter<double>("cutTkIsoRel") is a predefined template 
 which returns some number used somewhere calles here ( basically in a python) 
 https://github.com/mumuhamm/UserCode/blob/master/OmtfAnalysis/python/omtfTree_cfi.py#L82
*/

// Now we plot the histogram 

auto mjpsi =(TH1F*)data->createHistogram("mjpsi", *massJpsi); // auto works because again the systemwide installation 

/*So the first arguments is the name you give : you could use Karim, Piotr, Lalu, Champarani, Sadulla, Fajli, Illish, 
VerZara, whatever. As the constructor says a const char : user defined: and your collaborator might ask : did you checked 
the base histogram before you fit, there must be some technique to translate the RooDataSet to construct histogram, and
most of the time the nottheLevel one answers : I will look into this, unfortunately the basic training is flawed, the student 
does not know where to look into, it has no relation with the students brilliance or knowledge, again prabachan but it is true. 
Now we have the histogram , we can plot it, this part is easy, we need a canvas and a draw option (try this), but we will use RooPlotFrame*/

RooPlot* mjpsi_rp = massJpsi->frame(Title("m_{J/#psi} (GeV)"),Bins(nbins)); // the mass jpsi is RooRealVar
data->plotOn(mjpsi_rp,DataError(RooAbsData::SumW2));// this part is not neccesary for the moment , we need this to plot data and fir function on the same canvas. 
TCanvas * massCanvas = new TCanvas();
mjpsi_rp->Draw();


 
// model to fit the data 

RooRealVar * my_mean = new RooRealVar("my_mean","m_{J/#psi}", 3.096, 3.0, 3.2);
RooRealVar * my_sigma = new RooRealVar("my_sigma","sigma_{J/#psi}", 0, 0.0003);
RooGaussian *my_gaus = new RooGaussian("my_gaus","my_gaus", *massJpsi, *my_mean,*my_sigma);

/* Here we are using the same roorealVar for the mass distribution since this would be  used as the mass parameter */

RooFitResult* fitRes = my_gaus->fitTo(*data,Save());//data_SigReg
fitRes->Print("v");


/*etokhon preproduction shooting , editing , r postproduction holo , kintu climax hoche asol jaiga jokhon movie release korbe, 
Kintu golpo hoche tiktok dekha ami ki r movie dekhbo, dil hai tumhara gaane malai kebab, ba arijit er gaane , yuvraj er chota choi 
Eheno topology jara kore tara baghawan rabidas, by the way toplogy kara kore ? jara cup r donought er difference bojhe na 
Any way lets see how the fit looks : 

*/

    




}