#ifndef _DELPHESBTAGGER_H_
#define _DELPHESBTAGGER_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <map>
#include <vector>

#include "TString.h"
#include "TTree.h"
#include "TDataType.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TRandom3.h"

#include "Report.h"

class
DelphesBtagger
{

private:

  TRandom3 r;

  TF1 *func0B;
  TF1 *func1B;
  TF1 *func2B;
  TF1 *func3B;
  TF1 *func4B;

  TF1 *func0C;
  TF1 *func1C;
  TF1 *func2C;
  TF1 *func3C;
  TF1 *func4C;

  TF1 *func0L;
  TF1 *func1L;
  TF1 *func2L;
  TF1 *func3L;
  TF1 *func4L;

public:

  DelphesBtagger( Int_t rseed = 8675309 )
    : r(8675309)
    , func0B( new TF1("func0B","0.85*TMath::TanH(0.0026*x[0])*(30.0/(1+0.063*x[0]))",0,500) )
    , func1B( new TF1("func1B","0.84*TMath::TanH(0.0025*x[0])*(28.0/(1+0.068*x[0]))",0,500) )
    , func2B( new TF1("func2B","0.82*TMath::TanH(0.0024*x[0])*(27.0/(1+0.07*x[0]))",0,500) )
    , func3B( new TF1("func3B","0.75*TMath::TanH(0.0023*x[0])*(25.0/(1+0.072*x[0]))",0,500) )
    , func4B( new TF1("func4B","0.7*TMath::TanH(0.0022*x[0])*(25.0/(1+0.077*x[0]))",0,500) )
    , func0C( new TF1("func0C","0.25*TMath::TanH(0.018*x[0])*(1/(1+0.0013*x[0]))",0,500) )
    , func1C( new TF1("func1C","0.24*TMath::TanH(0.016*x[0])*(1/(1+0.0012*x[0]))",0,500) )
    , func2C( new TF1("func2C","0.23*TMath::TanH(0.014*x[0])*(1/(1+0.0011*x[0]))",0,500) )
    , func3C( new TF1("func3C","0.22*TMath::TanH(0.011*x[0])*(1/(1+0.0010*x[0]))",0,500) )
    , func4C( new TF1("func4C","0.20*TMath::TanH(0.008*x[0])*(1/(1+0.0009*x[0]))",0,500) )
    , func0L( new TF1("func0L","0.01*0.00038*x[0]",0,500) )
    , func1L( new TF1("func1L","0.008*0.00036*x[0]",0,500) )
    , func2L( new TF1("func2L","0.006*0.0003*x[0]",0,500) )
    , func3L( new TF1("func3L","0.003*0.00025*x[0]",0,500) )
    , func4L( new TF1("func4L","0.001*0.0001*x[0]",0,500) )
  {}

  ~DelphesBtagger() {
    delete func0B;
    delete func1B;
    delete func2B;
    delete func3B;
    delete func4B;
    delete func0C;
    delete func1C;
    delete func2C;
    delete func3C;
    delete func4C;
    delete func0L;
    delete func1L;
    delete func2L;
    delete func3L;
    delete func4L;
  }

  Int_t getTagLevel( const float &pt , const int &flavor ) {

    Int_t tagLevel = 0;
    float efficiency[5];
    //P(tight | loose) = P(loose | tight) * P(tight)/ P (loose) = P(tight)/ P (loose)

    if (flavor == 5) {
      efficiency[0] = func0B->Eval(pt);//loosest
      efficiency[1] = func1B->Eval(pt);
      efficiency[2] = func2B->Eval(pt);
      efficiency[3] = func3B->Eval(pt);
      efficiency[4] = func4B->Eval(pt);//tightest
    }
    else if (flavor == 4) {
      efficiency[0] = func0C->Eval(pt);//loosest
      efficiency[1] = func1C->Eval(pt);
      efficiency[2] = func2C->Eval(pt);
      efficiency[3] = func3C->Eval(pt);
      efficiency[4] = func4C->Eval(pt);//tightest
    }
    else {
      efficiency[0] = func0L->Eval(pt);//loosest
      efficiency[1] = func1L->Eval(pt);
      efficiency[2] = func2L->Eval(pt);
      efficiency[3] = func3L->Eval(pt);
      efficiency[4] = func4L->Eval(pt);//tightest
    }
    float uniform = r.Uniform(1.0);

    if( uniform < efficiency[0] ) {
      tagLevel = 1;
      float looseGivenTight = efficiency[1] / efficiency[0];
      uniform = r.Uniform(1.0);
      if( uniform < looseGivenTight ) {
        tagLevel = 2;
        looseGivenTight = efficiency[2] / efficiency[1];
        uniform = r.Uniform(1.0);
        if( uniform < looseGivenTight ) {
          tagLevel = 3;
          looseGivenTight = efficiency[3] / efficiency[2];
          uniform = r.Uniform(1.0);
          if( uniform < looseGivenTight ) {
            tagLevel = 4;
            looseGivenTight = efficiency[4] / efficiency[3];
            uniform = r.Uniform(1.0);
            if( uniform < looseGivenTight ) {
              tagLevel = 5;
            }
          }
        }
      }
    }

    return tagLevel;
  }

};



#endif
