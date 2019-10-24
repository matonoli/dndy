// @(#)root/physics:$Id$
// Author: Adrian Bevan  2001

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * Copyright (C) 2001, Liverpool University.                             *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFeldmanCousins2
    \ingroup Physics

Class to calculate the CL upper limit using
the Feldman-Cousins method as described in PRD V57 #7, p3873-3889

The default confidence interval calculated using this method is 90%
This is set either by having a default the constructor, or using the
appropriate fraction when instantiating an object of this class (e.g. 0.9)

The simple extension to a gaussian resolution function bounded at zero
has not been addressed as yet -> `time is of the essence' as they write
on the wall of the maze in that classic game ...

#### VARIABLES THAT CAN BE ALTERED
=> depending on your desired precision: The initial values of fMuMin,
fMuMax, fMuStep and fNMax are those used in the PRD:
~~~ {.cpp}
  fMuMin = 0.0
  fMuMax = 50.0
  fMuStep= 0.005
~~~
but there is total flexibility in changing this should you desire.


see example of use in $ROOTSYS/tutorials/math/FeldmanCousins.C

see note about: "Should I use TRolke, TFeldmanCousins2, TLimit?"
 in the TRolke class description.

\author: Adrian Bevan, Liverpool University

Copyright Liverpool University 2001       bevan@slac.stanford.edu
*/

#include <iostream>
#include "TMath.h"
#include "TFeldmanCousins2.h"

ClassImp(TFeldmanCousins2);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TFeldmanCousins2::TFeldmanCousins2(Double_t newFC, TString options)
{
   fCL          = newFC;
   fUpperLimit  = 0.0;
   fLowerLimit  = 0.0;
   fNobserved   = 0.0;
   fNbackground = 0.0;
   options.ToLower();
   if (options.Contains("q")) fQUICK = 1;
   else                       fQUICK = 0;

   fNMax   = 50;
   fMuStep = 0.005;
   SetMuMin();
   SetMuMax();
   SetMuStep();
}


////////////////////////////////////////////////////////////////////////////////
// Destructor.

TFeldmanCousins2::~TFeldmanCousins2()
{
}

////////////////////////////////////////////////////////////////////////////////
/// given Nobserved and Nbackground, try different values of mu that give lower limits that
/// are consistent with Nobserved.  The closed interval (plus any stragglers) corresponds
/// to the F&C interval

Double_t TFeldmanCousins2::CalculateLowerLimit(Double_t Nobserved, Double_t Nbackground, Double_t Uncertainty)
{
   CalculateUpperLimit(Nobserved, Nbackground, Uncertainty);
   return fLowerLimit;
}


////////////////////////////////////////////////////////////////////////////////
/// given Nobserved and Nbackground, try different values of mu that give upper limits that
/// are consistent with Nobserved.  The closed interval (plus any stragglers) corresponds
/// to the F&C interval

Double_t TFeldmanCousins2::CalculateUpperLimit(Double_t Nobserved, Double_t Nbackground, Double_t Uncertainty)
{
   fNobserved   = Nobserved;
   fNbackground = Nbackground;
  fUnc = Uncertainty;

   Double_t mu = 0.0;

   // for each mu construct the ranked table of probabilities and test the
   // observed number of events with the upper limit
   Double_t min = -999.0;
   Double_t max = 0;
   Int_t iLower = 0;

   Int_t i;
   for(i = 0; i <= fNMuStep; i++) {
      mu = fMuMin + (Double_t)i*fMuStep;
      Int_t goodChoice = FindLimitsFromTable( mu );
      if( goodChoice ) {
         min = mu;
         iLower = i;
         break;
      }
   }

   //==================================================================
   // For quicker evaluation, assume that you get the same results when
   // you expect the uppper limit to be > Nobserved-Nbackground.
   // This is certainly true for all of the published tables in the PRD
   // and is a reasonable assumption in any case.
   //==================================================================

   Double_t quickJump = 0.0;
   if (fQUICK)          quickJump = Nobserved-Nbackground-fMuMin;
   if (quickJump < 0.0) quickJump = 0.0;

   for(i = iLower+1; i <= fNMuStep; i++) {
      mu = fMuMin + (Double_t)i*fMuStep + quickJump;
      Int_t goodChoice = FindLimitsFromTable( mu );
      if( !goodChoice ) {
         max = mu;
         break;
      }
   }

   fUpperLimit = max;
   fLowerLimit = min;

   return max;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the probability table for a given mu for n = 0, NMAX
/// and return 1 if the number of observed events is consistent
/// with the CL bad

Int_t TFeldmanCousins2::FindLimitsFromTable( Double_t mu )
{
   Double_t *p          = new Double_t[fNMax];   //the array of probabilities in the interval MUMIN-MUMAX
   Double_t *r          = new Double_t[fNMax];   //the ratio of likliehoods = P(Mu|Nobserved)/P(MuBest|Nobserved)
   Int_t    *rank       = new Int_t[fNMax];      //the ranked array corresponding to R (largest first)
   Double_t *muBest     = new Double_t[fNMax];
   Double_t *probMuBest = new Double_t[fNMax];

   //calculate P(i | mu) and P(i | mu)/P(i | mubest)
   Int_t i;
   for(i = 0; i < fNMax; i++) {
      muBest[i] = (Double_t)(i - fNbackground);
      if(muBest[i]<0.0) muBest[i] = 0.0;
      probMuBest[i] = Prob(i, muBest[i],  fNbackground);
      p[i]          = Prob(i, mu,  fNbackground, fUnc);
      if(probMuBest[i] == 0.0) r[i] = 0.0;
      else                     r[i] = p[i]/probMuBest[i];
   }

   //rank the likelihood ratio
   TMath::Sort(fNMax, r, rank);

   //search through the probability table and get the i for the CL
   Double_t sum = 0.0;
   Int_t iMax = rank[0];
   Int_t iMin = rank[0];
   for(i = 0; i < fNMax; i++) {
      sum += p[rank[i]];
      if(iMax < rank[i]) iMax = rank[i];
      if(iMin > rank[i]) iMin = rank[i];
      if(sum >= fCL) break;
   }

   delete [] p;
   delete [] r;
   delete [] rank;
   delete [] muBest;
   delete [] probMuBest;

   if((fNobserved <= iMax) && (fNobserved >= iMin)) return 1;
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the poissonian probability for a mean of mu+B events with a variance of N taking the systematic uncertainty into account.
Double_t TFeldmanCousins2::Prob(Int_t N, Double_t mu, Double_t B, Double_t Unc){

  if (Unc <= 0) { // belwo this value the result will be wrong
    return TMath::Poisson( N, mu+B);
  } else {
    TF1* fPDF = new TF1("fPDF","TMath::Poisson([0],x*([1]+[2]))*TMath::Gaus(x,1.,[3],1)",0.,2.);
    fPDF->SetParameters(N,mu,B,Unc);
    Double_t LowerBoundary = 1 - 5*Unc;
    if (LowerBoundary < 0) LowerBoundary = 0;
    return fPDF->Integral(LowerBoundary, 1 + 5*Unc);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Set maximum value of signal to use in calculating the tables.

void TFeldmanCousins2::SetMuMax(Double_t newMax)
{
   fMuMax   = newMax;
   fNMax    = (Int_t)newMax;
   SetMuStep(fMuStep);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the step in signal to use when generating tables.

void TFeldmanCousins2::SetMuStep(Double_t newMuStep)
{
   if(newMuStep == 0.0) {
      std::cout << "TFeldmanCousins2::SetMuStep ERROR New step size is zero - unable to change value"<< std::endl;
      return;
   } else {
      fMuStep = newMuStep;
      fNMuStep = (Int_t)((fMuMax - fMuMin)/fMuStep);
   }
}
