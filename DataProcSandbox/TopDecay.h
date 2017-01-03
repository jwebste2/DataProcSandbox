#ifndef _TOPDECAY_H_
#define _TOPDECAY_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <map>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"

namespace
topdecay
{

  // Top quark decays
  const int WB	 = 1;
  const int WLIGHT = 2;
  
  // W boson decays
  const int JETS  = 1;
  const int ELNU  = 2;
  const int MUNU  = 3;
  const int TAUNU = 4;

  // Undefined decay
  const int UNDEFINED = 0;

  //
  // Functions for initializing decay matrices
  //
  
  TH2D* getEmptyTopDecayMatrix() {
    TH2D *h2 = new TH2D( "h2_topDecayCategories" , ";t Decay;#bar{t} Decay;Rate" , 2 , 0.5 , 2.5 , 2 , 0.5 , 2.5 );
    h2->GetXaxis()->SetBinLabel( WB , "W^{+}+b" );
    h2->GetXaxis()->SetBinLabel( WLIGHT , "W^{+}+light q" );
    h2->GetYaxis()->SetBinLabel( WB , "W^{-}+#bar{b}" );
    h2->GetYaxis()->SetBinLabel( WLIGHT , "W^{-}+light #bar{q}" );
    return h2;
  }

  TH2D* getEmptyWDecayMatrix() {
    TH2D *h2 = new TH2D( "h2_WDecayCategories" , ";W^{+} Decay;W^{-} Decay;Rate" , 4 , 0.5 , 4.5 , 4 , 0.5 , 4.5 );
    h2->GetXaxis()->SetBinLabel( JETS , "jets" );
    h2->GetXaxis()->SetBinLabel( ELNU , "e#nu" );
    h2->GetXaxis()->SetBinLabel( MUNU , "#mu#nu" );
    h2->GetXaxis()->SetBinLabel( TAUNU , "#tau#nu" );
    h2->GetYaxis()->SetBinLabel( JETS , "jets" );
    h2->GetYaxis()->SetBinLabel( ELNU , "e#nu" );
    h2->GetYaxis()->SetBinLabel( MUNU , "#mu#nu" );
    h2->GetYaxis()->SetBinLabel( TAUNU , "#tau#nu" );
    return h2;
  }

};


#endif
