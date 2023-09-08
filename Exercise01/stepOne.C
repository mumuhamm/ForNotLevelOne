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

/*This rooreal var is a constaructor , even though you do not instantiate the base class : your system 
wide installation of the roofit sofware ,and source to the roofit library takes care of the instantiation 
of the base class , and they create a "onfly" library-load called gSystem which automatically loads the library. 
This constructor takes four inputs (arguments) : first two are const char( consider static string, but you could also 
ask for the character names somewhere in your code, for instance when a roorealvar or a histogram is filled (or initiated like Ben) you 
could just ask myhisto->GetName()and then you could use the same histogram (roorealvar) for further anallysis and/or fit) and the next two are double/float.
*/

//The model which will generate the mass histogram 
RooRealVar * mean = new RooRealVar("mean","m_{J/#psi}", 3.096, 3.0, 3.2);
RooRealVar * sigma = new RooRealVar("sigma","sigma_{J/#psi}", 0, 0.0003);
RooGaussian *gaus = new RooGaussian("gaus","gaus", *massJpsi, *mean,*sigma);

/* the Gaussian model takes five arguments , so all these base class references can be found in the root site. You could just search for 
some string name : TH1F class reference in Google. 
Using a pointer does not make any difference in physics or math, its the C++ thing, I always miss the point ( tiny, on the ... arrow is more visible).
for revisiting C++ functions please visit cplusplus or programize. Now, the model is ready, so, we can use some root functinality to create some random data */

RooDataSet* data = gaus->generate(*massJpsi,100000);

/*It takes the mass of the particle (the physics observable in question), but see how does it code , you fill the rooreal var, 
when people talk , they talk in the language of physics, which makes more sense, so in this case someone would say : we could use 
the Jpsi mass to generate some random events, and experienced people are more inclined toward talking in the language of physics, so the
idea is to translating all lines of code in the language of physics. So, we used the Jpsi mass/resonance within the scope of its uncertainty, 
though it is not exact, then a natural width is assigned, even that is not exact too, 0.3MeV is way more biased than what PDG provides , 
jpsi width is in KeV, but due to the misreconstruction of the tracks the natural resolution decreases, thats why the leverage. Even more expert
people sometimes say , it seems that some of the bins in the lower side band is overflowed and the pile-up correction/reweight factor seems have some bias. 
Is it integrated or differential in bins of some X variable ? What would be the answer to this question? And also, 
How many events : 100000
exercise : please tanslate the following c++ lines : into physics language ( from one of my muonselection code) :
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

/*So the first arguments is the const char name you give : you could use Karim, Piotr, Lalu, Champarani, Sadulla, Fajli, Illish, 
VerZara, whatever. As the constructor says a const char : "user defined": and your collaborator might ask : did you checked 
the base histogram before you fit, there must be some technique to translate the RooDataSet to construct histogram, and
most of the time the nottheLevelone answers : I will look into this, unfortunately the basic training is flawed, the student 
does not know where to look into, it has no relation with the students brilliance or knowledge, again prabachan but it is true. 
Now we have the histogram , we can plot it, this part is easy, we need a canvas and a draw option (try this), but we will use RooPlotFrame*/

RooPlot* mjpsi_rp = massJpsi->frame(Title("m_{J/#psi} (GeV)"),Bins(nbins)); // the mass jpsi is RooRealVar
data->plotOn(mjpsi_rp,DataError(RooAbsData::SumW2));// this part is not neccesary for the moment , we need this to plot data and fir function on the same canvas. 
TCanvas * massCanvas = new TCanvas();
mjpsi_rp->Draw();

/*Question what this function does :DataError(RooAbsData::SumW2), in case if we do not use , what changes do we expect in the results. */

 
/*model to fit the data : a Gaussin model , since we know the source , for real data , we creat our own function , looking at 
the distribution and from the knowledge of mathematics. Root provides a lot of basic functions, */ 

RooRealVar * my_mean = new RooRealVar("my_mean","m_{J/#psi}", 3.096, 3.0, 3.2);
RooRealVar * my_sigma = new RooRealVar("my_sigma","#sigma_{J/#psi}", 0, 0.0003);
RooGaussian *my_gaus = new RooGaussian("my_gaus","my_gaus", *massJpsi, *my_mean,*my_sigma);

/* Here we are using the same roorealVar for the mass distribution since this would be  used as the mass parameter */

RooFitResult* fitRes = my_gaus->fitTo(*data,Save());//data_SigReg
fitRes->Print("v");


/*etokhon preproduction shooting , editing , r postproduction holo , kintu climax hoche asol jaiga jokhon movie release korbe, 
Kintu golpo hoche tiktok dekha ami ki r movie dekhbo, dil hai tumhara gaane malai kebab, ba arijit er gaane , yuvraj er chota choi 
Eheno topology jara kore tara baghawan rabidas, by the way toplogy kara kore ? jara cup r donought er difference bojhe na 
Any way lets see how the fit looks : 
*/

//We already have a rooplot and data is already is plotted on it, so we are going to create another one. Line 135    
   my_gaus->plotOn(mjpsi_rp, LineColor(kBlue), LineWidth(1));//plotting fit function on top of the data
   my_gaus->paramOn(mjpsi_rp); // Fit parameters 
   
   
   RooPlot* pullmass = massJpsi->frame(RooFit::Title("mass pull")); //rooreal var mass is used 
   RooHist* hpull1  = mjpsi_rp->pullHist();
   pullmass->addPlotable(hpull1,"P0") ;
   pullmass->SetMinimum(-3) ;
   pullmass->SetMaximum(+3) ;
   pullmass->SetYTitle("pull");
   pullmass->SetMarkerStyle(20);
   pullmass->SetNdivisions(10);
   Double_t chisquare_mjpsi_rp = mjpsi_rp->chiSquare();
   cout<<"Chi square of mass fit is :"<< chisquare_mjpsi_rp<< endl;

/*Here I am going to divide one canvas in to two pads  and then fill the
upper pad with fitten mass distribution and then the lower one with pull distrubition, The notALevelOne 
student is always welcome to decorate the way the staudent wants. 
for instance for CMS size - I would use 
https://github.com/mumuhamm/RootAnalysis/blob/devel_ANYDATAFORMAT/GMT/Analysis/src/GMTHistograms.cc#L754
also for more descritive roofit example one would visit some more example : 
I have a repo with some small codes I needed to have for some test used in my thesis: 
for example a simple two dimensional momentum morphing. 
  */

   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   TPad *pad1 = new TPad("pad1","pad1",0,0.36,1,1);
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
   pad1->SetBottomMargin(0.00001);
   pad1->SetBorderMode(0);
   pad2->SetTopMargin(0.00001);
   pad2->SetBottomMargin(0.1);
   pad2->SetBorderMode(0);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();
   gStyle->SetOptTitle(0);
   c->SetFillColor(0);
   c->SetBorderSize(2);
   c->SetLeftMargin(0.1422222);
   c->SetRightMargin(0.04444445);
   mjpsi_rp->SetStats(0);
   mjpsi_rp->Draw();
   pad2->cd();
   pullmass->SetStats(0);
   pullmass->Draw();
   c->cd();


/*We will launch the likelihood, in math the region of convergence , or very near to the minima , where the derivative 
of the function exists, kind of analytic, that region would be the likelihood, analytically easy, take a function
take log to make it to the upper bound multivalued additive function, restricts the products to converge into zero by coverting those 
into addition, and I forgot the term some statibility. */

 Double_t amin,edm,errdef;
 Int_t nvpar,nparx;
 RooRealVar* mean_var = (RooRealVar*)fitRes->floatParsFinal().find("my_mean");
 RooRealVar* sigma_var = (RooRealVar*)fitRes->floatParsFinal().find("my_sigma");

 /*depending upon the availibity of the methods or functions in a class you can take the output you want*/
Float_t mean_central, sigma_central;
Float_t mean_err, sigma_err;
   
mean_central = mean_var->getVal();
mean_err = mean_var->getAsymErrorHi();
   
sigma_central = sigma_var->getVal();
sigma_err = sigma_var->getAsymErrorHi();
   
std::cout<< "--------------------------------------------------------------"<<"\n";   
std::cout<<"fitpar : "<<mean_central<<"\t"<<sigma_central<<"\n";
std::cout<<"fitpar_err : "<<mean_err<<"\t"<<sigma_err<<"\n";
std::cout<< "--------------------------------------------------------------" <<"\n";  

mean_err /= 4.;
sigma_err /= 4.;


//Negative log likelihood

RooAbsReal* nll = my_gaus->createNLL(*data,NumCPU(8));
nll->SetName("nll");
RooMinuit minuit(*nll);
minuit.migrad();




/* the LL calculation is simple , it just scanning grid of 20 bins , in two loop inwardPlease decorate the plot, this code is inpire by a code I wrote to calculate 
Wilson coefficient for some unpublished work, so I know it is the right procedure but still please send me your comments*/
const double nllMin = nll->getVal();
   TString tmp;
   tmp += nllMin;
   tmp += " )";
   RooFormulaVar * likl = new RooFormulaVar("likl", "L", "exp(-nll  + "+tmp, RooArgList(*nll));
   TH2D * manL = new TH2D("manL", "likelihood contour",40,mean_central-(20*mean_err),mean_central+(20*mean_err),40,sigma_central-(20*sigma_err),sigma_central+(20*sigma_err));
   manL->SetXTitle(my_mean->GetTitle());
   manL->SetYTitle(my_sigma->GetTitle());
   Float_t mean_m;
   Float_t sigma_m;
   mean_m = mean_central-(20*mean_err);
   for(int i = 1; i <= manL->GetNbinsX();++i,mean_m+=mean_err)
      {
      sigma_m = sigma_central-(20*sigma_err);
      for(int j = 1; j <= manL->GetNbinsY();++j,sigma_m+=sigma_err)
         {
         
         double x = manL->GetXaxis()->GetBinCenter(i);
         my_mean->setVal(x);
         double y = manL->GetYaxis()->GetBinCenter(j);
         my_sigma->setVal(y);
         
         TVirtualFitter *fitter = TVirtualFitter::Fitter(mjpsi);
         fitter->GetStats(amin,edm,errdef,nvpar,nparx);
         manL->SetBinContent(i, j, amin-likl->getVal());
         
         }
      }
   TCanvas *con = new TCanvas("con", "con", 0, 0, 800,600);
   manL->Draw("cont4z");


   /*Here we actually have done the exercise from 1-6, but this seems more of a roofit code than a tutorial of C++, so there is this question what sort C++ function
   can make this code more smarter, one that I wrote in the question, 7 which would make some of the hardcoded part of the code into function, so where we can put some arguments 
   and it will go like "onfly", but then what what would be the answer for the question 8 ? Please share your idea to change some of the hardcoded part into functions.*/

}