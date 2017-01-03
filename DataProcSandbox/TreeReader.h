#ifndef _TREEREADER_H_
#define _TREEREADER_H_

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

#include "Report.h"
#include "tth.h"

#include "TTHbbLeptonic/MVAVariables.h"
#include "TTHbbLeptonic/PairedSystem.h"
#include "TopEvent/Event.h"

class
TreeReader
{

protected:

  const Double_t Zmass = 91187.6;

  TTree *tree;

  std::vector<TTree*> outputTrees;

  // When we are scanning an input tree there's some extra output that we want to
  // produce only for signal. This bool is used to turn on/off this extra functionality.
  Bool_t signalMode;
  
  top::Event m_event;
  
  enum BranchType {
    UNDEFINED=0 ,
    UINT=1 ,
    INT=2 ,
    LONG64=3 ,
    FLOAT=4 ,
    CHAR=5 ,
    VECTOR_FLOAT=6 ,
    VECTOR_CHAR=7 ,
    VECTOR_INT=8 ,
    VECTOR_UINT=9
  };

  const std::vector<TString> vectortags = { "" , "_SIZE" , "_ENTRY0" , "_ENTRY1" , "_ENTRY2" , "_ENTRY3" };

  std::map< TString , BranchType > types;

  std::map< TString , UInt_t > map_uint;
  std::map< TString , Int_t > map_int;
  std::map< TString , Long64_t > map_long64;
  std::map< TString , Float_t > map_float;
  std::map< TString , Char_t > map_char;
  std::map< TString , std::vector<float>* > map_vector_float;
  std::map< TString , std::vector<char>* > map_vector_char;
  std::map< TString , std::vector<int>* > map_vector_int;
  std::map< TString , std::vector<uint>* > map_vector_uint;

  std::map< TString , Int_t > cutflow;
  
  //
  // Trigger config structure for taking care of trigger matching
  // for good leptons in each event. User defines number of matches
  // required for electrons/muons for each trigger.
  //
  struct Trig {

    TString name;
    TString elMatchName;
    TString muMatchName;
    Int_t nElectronsL1; // number of required electron matches
    Int_t nMuonsL1; // number of required muon matches

    Trig( TString _name , Int_t _nElectronsL1 , Int_t _nMuonsL1 )
      : name( _name )
      , elMatchName( "el_trigMatch_"+_name )
      , muMatchName( "mu_trigMatch_"+_name )
      , nElectronsL1( _nElectronsL1 )
      , nMuonsL1( _nMuonsL1 )
    {}
    
  };


  void setBranch_uint( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_uint[branchname] );
    }
    types[branchname] = UINT;
  }

  void setBranch_int( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_int[branchname] );
    }
    types[branchname] = INT;
  }

  void setBranch_long64( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_long64[branchname] );
    }
    types[branchname] = LONG64;
  }
  
  void setBranch_float( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_float[branchname] );
    }
    types[branchname] = FLOAT;
  }

  void setBranch_char( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_char[branchname] );
    }
    types[branchname] = CHAR;
  }

  void setBranch_vector_float( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_vector_float[branchname] );
    } else if( ! map_vector_float[branchname] ) {
      map_vector_float[branchname] = new std::vector<float>();
    }
    for( TString vtag : vectortags ) types[branchname+vtag] = VECTOR_FLOAT;
  }

  void setBranch_vector_char( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_vector_char[branchname] );
    } else if( ! map_vector_char[branchname] ) {
      map_vector_char[branchname] = new std::vector<char>();
    }
    for( TString vtag : vectortags ) types[branchname+vtag] = VECTOR_CHAR;
  }

  void setBranch_vector_int( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_vector_int[branchname] );
    } else if( ! map_vector_int[branchname] ) {
      map_vector_int[branchname] = new std::vector<int>();
    }
    for( TString vtag : vectortags ) types[branchname+vtag] = VECTOR_INT;
  }

  void setBranch_vector_uint( TString branchname , Bool_t intree = kTRUE ) {
    if( intree ) {
      tree->SetBranchStatus( branchname , 1 );
      tree->SetBranchAddress( branchname , &map_vector_uint[branchname] );
    } else if( ! map_vector_uint[branchname] ) {
      map_vector_uint[branchname] = new std::vector<uint>();
    }
    for( TString vtag : vectortags ) types[branchname+vtag] = VECTOR_UINT;
  }
  
  Double_t vectorElement( const TString &branchname , const std::size_t &i ) {
    if( types[branchname]==VECTOR_FLOAT ) return map_vector_float[branchname]->size() > i ? Double_t(map_vector_float[branchname]->at(i)) : -1;
    if( types[branchname]==VECTOR_INT ) return map_vector_int[branchname]->size() > i ? Double_t(map_vector_int[branchname]->at(i)) : -1;
    if( types[branchname]==VECTOR_UINT ) return map_vector_uint[branchname]->size() > i ? Double_t(map_vector_uint[branchname]->at(i)) : -1;
    if( types[branchname]==VECTOR_CHAR ) return map_vector_char[branchname]->size() > i ? Double_t( map_vector_char[branchname]->at(i)) : -1;
    report::error( "Failed to find branchname %s when looking for vector element %i" , branchname.Data() , int(i) );
    assert( false );
  }

  std::size_t vectorSize( const TString &branchname ) {
    if( types[branchname]==VECTOR_FLOAT ) return map_vector_float[branchname]->size();
    if( types[branchname]==VECTOR_INT ) return map_vector_int[branchname]->size();
    if( types[branchname]==VECTOR_UINT ) return map_vector_uint[branchname]->size();
    if( types[branchname]==VECTOR_CHAR ) return map_vector_char[branchname]->size();
    report::error( "Failed to find branchname %s when looking for vector->size" , branchname.Data() );
    assert( false );
  }

  Double_t read( const TString &branchname ) {
    Double_t result;
    fill( branchname , result );
    return result;
  }

  void setTree_Generic( TTree *_tree ) {
    // Declaration of leaf types
    tree = _tree;
    tree->SetBranchStatus( "*" , 0 );
    
    TIter next  = _tree->GetListOfLeaves();
    TLeaf *leaf = 0;
    while( (leaf = (TLeaf*)next()) ) {
      TString leafname  = leaf->GetName();
      TString leaftitle = leaf->GetTitle();
      TString leaftype  = leaf->GetTypeName();
      //Bool_t isarray    = leaftitle.Contains("[") && leaftitle.Contains("]");
      if( leafname.EndsWith("_") ) continue;
      if( leafname.BeginsWith("ljet_") ) continue;
      else if( leaftype == TString("Float_t") ) setBranch_float( leafname );
      else if( leaftype == TString("Int_t") ) setBranch_int( leafname );
      else if( leaftype == TString("UInt_t") ) setBranch_uint( leafname );
      else if( leaftype == TString("Long64_t") ) setBranch_long64( leafname );
      else if( leaftype == TString("Char_t") ) setBranch_char( leafname );
      else if( leaftype == TString("vector<float>") ) setBranch_vector_float( leafname );
      else if( leaftype == TString("vector<char>") ) setBranch_vector_char( leafname );
      else if( leaftype == TString("vector<int>") ) setBranch_vector_int( leafname );
      else if( leaftype == TString("vector<uint>") ) setBranch_vector_uint( leafname );
      else if( leaftype == TString("TRef") ) continue;
      else if( leaftype == TString("TLorentzVector") ) continue;
      else if( leaftype == TString("ULong64_t") ) continue;
      else {
	report::error( "Not sure how to handle unrecognized leaf type=%s (leaf name=%s)" , leaftype.Data() , leafname.Data() );
	report::blank( "If you want to simply ignore this leaftype, add the following line in TreeReader::setTree_Generic:" );
	report::blank( "  else if( leaftype == TString(\"%s\") ) continue;" , leaftype.Data() );
	assert( false );
      }
    }

  }

  typedef std::map<variable,TString>::const_iterator variableIter;
  typedef std::map<pairing,TString>::const_iterator pairingIter;
  typedef std::map<collection,TString>::const_iterator collectionIter;
  
  const std::map<variable,TString> mva_variable_names = {
    { variable::MaxPt , "MaxPt" } ,
    { variable::MinPt , "MinPt" } ,
    { variable::MaxdEta , "MaxdEta" } ,
    { variable::MindEta , "MindEta" } ,
    { variable::MaxdR , "MaxdR" } ,
    { variable::MindR , "MindR" } ,
    { variable::MaxM , "MaxM" } ,
    { variable::MinM , "MinM" } ,
    { variable::Avg , "Avg" } ,
    { variable::HiggsMass , "HiggsMass" } };

  const std::map<pairing,TString> mva_pairing_names = {
    { pairing::jj , "jj" } ,
    { pairing::bj , "bj" } ,
    { pairing::bb , "bb" } ,
    { pairing::uu , "uu" } ,
    { pairing::jl , "jl" } ,
    { pairing::lepbb , "lepbb" } ,
    { pairing::Hl , "Hl" } };
    
  const std::map<collection,TString> mva_collection_names = {
    { collection::all , "all" } ,
    { collection::jets , "jets" } ,
    { collection::bjets , "bjets" } ,
    { collection::ujets , "ujets" } };

  
public:

 
  TreeReader() {}
  ~TreeReader() {}

  Double_t operator[]( const TString &branchname ) { return read(branchname); }

  // Set input "val" to the value corresponding to "branchname"
  // Return false requesting a vector element that doesn't exist
  Bool_t fill( const TString  &branchname , Double_t &val ) {
    //report::debug( "filling %s" , branchname.Data() );
    BranchType type = types[branchname];
    if( type==UINT ) { val = Double_t(map_uint[branchname]); return kTRUE; }
    if( type==INT ) { val = Double_t(map_int[branchname]); return kTRUE; }
    if( type==LONG64 ) { val = Double_t(map_long64[branchname]); return kTRUE; }
    if( type==FLOAT ) { val = Double_t(map_float[branchname]); return kTRUE; }
    if( type==CHAR ) { val = Double_t(map_char[branchname]); return kTRUE; }
    if( type==VECTOR_FLOAT || type==VECTOR_INT || type==VECTOR_CHAR ) {
      TString origname( branchname );
      if( branchname.EndsWith("_SIZE") ) { val = Double_t( vectorSize(origname.ReplaceAll("_SIZE","")) ); return kTRUE; }
      if( branchname.EndsWith("_ENTRY0") ) { val = vectorElement( origname.ReplaceAll("_ENTRY0","") , 0 ); return vectorSize(origname) > 0; }
      if( branchname.EndsWith("_ENTRY1") ) { val = vectorElement( origname.ReplaceAll("_ENTRY1","") , 1 ); return vectorSize(origname) > 1; }
      if( branchname.EndsWith("_ENTRY2") ) { val = vectorElement( origname.ReplaceAll("_ENTRY2","") , 2 ); return vectorSize(origname) > 2; }
      if( branchname.EndsWith("_ENTRY3") ) { val = vectorElement( origname.ReplaceAll("_ENTRY3","") , 3 ); return vectorSize(origname) > 3; }
      report::printassert( false , "branch %s not identified" , branchname.Data() );
    }
    report::error( "TreeReader cant find entry %s" , branchname.Data() );
    assert( false );
  }

  
  void fillHist( TH1* h , const TString &branchname , const Double_t &units , const Double_t &w ) {
    Double_t val;
    TString origname( branchname );
    if( branchname.EndsWith("_ALL") ) {
      origname.ReplaceAll( "_ALL" , "" );
      for( size_t i = 0 ; i < vectorSize(origname) ; ++i ) {
	h->Fill( vectorElement(origname,i) * units , w );
      }
    } else {
      if( fill( origname.ReplaceAll("_INCL","") , val ) ) {
	if( branchname.EndsWith("_INCL") ) {
	  for( Int_t ibin = 1 ; ibin <= h->GetNbinsX() ; ++ibin ) {
	    if( h->GetBinLowEdge(ibin) <= (val * units) ) {
	      h->Fill( h->GetBinCenter(ibin) , w );
	    }
	  }
	} else {
	  h->Fill( val * units , w );
	}
      }
    }
  }

  virtual void setupOutputTree( const TString & ) {}
  virtual void fillOutputTree( const std::size_t& , const Long64_t& , const Double_t& ) {}
  virtual void saveOutputTrees( const TString & ) {}
  void setSignalMode( const Bool_t &v ) { signalMode = v; }
  
  virtual TString getSelectionTitle() { return "TitleNotSet"; }
  virtual TString getDressedSelectionTitle() { return getSelectionTitle(); }

  virtual void debugInfo() {}
  
  virtual void initCutflow() {}
  virtual void printCutflow() {
    std::vector< std::pair<TString,Int_t> > cnts;
    if( cutflow.size() == 0 ) return;
    for( auto itr = cutflow.begin() ; itr != cutflow.end() ; ++itr )
      cnts.push_back( *itr );

    std::sort( cnts.begin() , cnts.end() ,
	       [=](std::pair<TString,Int_t>& a,std::pair<TString,Int_t>& b) {return a.second > b.second;}
	       );

    Double_t maxcnt = Double_t( cnts.front().second );

    report::info( "Cutflow:" );
    for( std::size_t icut = 0 ; icut < cnts.size() ; ++icut ) {
      Double_t eff = Double_t(cnts[icut].second) / maxcnt;
      report::blank( "%30s %20i %20g" , cnts[icut].first.Data() , cnts[icut].second , eff );
    }
  }

  //
  // Functions that need to be overloaded
  //
  virtual void setTree( TTree *_tree ) { assert( false ); }
  virtual Bool_t getEntry( const Long64_t &ientry ) { assert( false ); return kFALSE; }
  virtual Bool_t passesSelection() { assert( false ); return kFALSE; }
  virtual Double_t leadingLeptonPt() { return -666; }

};

#endif
