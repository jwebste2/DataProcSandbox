#ifndef _PHYSICSPROCESS_H_
#define _PHYSICSPROCESS_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"

#include "Report.h"
#include "HistTools.h"
#include "Config.h"

static Bool_t youHaveBeenWarned_PUWeights = kFALSE;


class
PhysicsProcess
{

private:

  enum PhysicsProcessType { UNDEFINED , DATA , BACKGROUND , SIGNAL };
  
  TString name;
  TString roottitle;
  TString latextitle;
  TString treename;
  Int_t color;
  PhysicsProcessType type;

  struct SampleInfo {

    TString tag;
    Double_t xsec;
    Double_t kfactor;
    Double_t nevents_sumWeights; // just used to test the sumWeights tree
    Double_t nevents_cutflow;
    Double_t nevents_cutflow_mc; // used to test the cutflow_mc numbers
    Double_t nevents_cutflow_mc_pu; // used to test the cutflow_mc_pu numbers
    Double_t nevents_cutflow_mc_pu_zvtx; // used to test the cutflow_mc_pu_zvtx numbers    
    Double_t weight;
    TFile *file;

    SampleInfo()
      : tag( "null" )
      , xsec( -1.0 )
      , kfactor( -1.0 )
      , nevents_sumWeights( -1.0 )
      , nevents_cutflow( -1.0 )
      , nevents_cutflow_mc( -1.0 )
      , nevents_cutflow_mc_pu( -1.0 )
      , nevents_cutflow_mc_pu_zvtx( -1.0 )
      , weight( 0.0 )
      , file( 0 )
    {}

    // function to search for tag in a list of directories
    Int_t matchTag( std::vector<TString> &dirnames , TString &match ) {
      Int_t matches( 0 );
      for( size_t i = 0 ; i < dirnames.size() ; i++ ) {
        if( dirnames[i].Contains(tag) ) {
          matches++;
          match = dirnames[i];
        }
      }
      return matches;
    }

    // There are a bunch of different total yields loaded from different cutflow histograms and the sumWeights
    // tree... this function specifies which one to actually use when normalizing the MC
    Double_t getNevents() const {
      if( !youHaveBeenWarned_PUWeights ) {
	report::debug( "Physics process is using the variable \"%s\" to normalize" , (config::UsePileupWeights?"nevents_cutflow_mc_pu_zvtx":"nevents_cutflow_mc") );
	report::blank( "the MC. This can be changed (e.g. if you want to remove/add pileup weights" );
	report::blank( "by editing the SampleInfo::getNevents function in PhysicsProcess.h" );
	report::debug( "IMPORTANT - If pileup weights are included in the calculation of this number, then" );
	report::blank( "you should be careful that you are also including pileup weights in your event-by-" );
	report::blank( "event weights, defined in ExtraPlotTools.h. This should be handled automatically by" );
	report::blank( "the bool config::UsePileupWeights." );
	report::blank( "" );
	report::warn( "currently using sumWeights tree to get initial sample yields for cutflow studies!!" );
	report::blank( "" );
	youHaveBeenWarned_PUWeights = kTRUE;
      }
      return nevents_sumWeights;
      if( config::UsePileupWeights ) {
	return nevents_cutflow_mc_pu_zvtx;
      } else {
	return nevents_cutflow_mc;
      }
      assert( false );
    }

    void report( TString pname , Double_t lumi = -1 ) {
      report::info( "PhysicsProcess %-30s : Loaded file for tag %-15s lumi = %-10g xsec = %-10g nev = %-10g w = %-15g" , pname.Data() , tag.Data() , lumi , xsec , getNevents() , weight );
    }

    void print() const { report::blank( "tag = %s, xsec = %g, kfactor = %g" , tag.Data() , xsec , kfactor ); }

  };

  std::vector<SampleInfo> samples;

  // Load MC files - this function should be used when the sample tag represents the full path to the file
  void loadAsMC( Double_t lumi ) {
    for( std::size_t i = 0 ; i < samples.size() ; i++ ) {
      // FIXME - add Nevents once I figure out how to handle this with Delphes samples
      samples[i].file = TFile::Open( samples[i].tag );
      samples[i].weight = lumi * samples[i].xsec * samples[i].kfactor / samples[i].getNevents();
      samples[i].report( name , lumi );
    }
  }

  // Load MC files - this function should be used when the tag just represents a part of the full path
  // The full path is looked up automagically
  void loadAsMC( TString databasedir , std::vector<TString> &datadirs , Double_t lumi ) {
    TString tmpdir;
    for( std::size_t i = 0 ; i < samples.size() ; i++ ) {

      if( samples[i].matchTag(datadirs,tmpdir) != 1 ) {
        report::error( "nmatched for tag %s = %i" , samples[i].tag.Data() , samples[i].matchTag(datadirs,tmpdir) );
	report::blank( "You could tell the InputConfig to ignore this sample tag by blacklisting it in tth.h" );
        assert( false );
      }

      TString filepath = databasedir + "/" + tmpdir + "/output.root";
      samples[i].file = TFile::Open( filepath );
      report::printassert( samples[i].file , "Couldn't open file %s for tag %s" , filepath.Data() , samples[i].tag.Data() );

      // Search the file for a TDirectory that contains cutflow histograms. The cutflow histograms
      // can be used to extract the total yield prior to any selection provided there is no skimming.
      TString tdirname = "";
      TIter iter       = samples[i].file->GetListOfKeys();
      TKey *key	       = (TKey*)0;
      while( tdirname.Length()==0 && (key = (TKey*)iter()) ) {
	TObject *obj    = key->ReadObj();
	TString objname = obj->GetName();
	if( obj->IsA()->InheritsFrom("TDirectory") ) tdirname = objname;
      }

      report::printassert( tdirname.Length() > 0 , "Failed to find TDirectory in file when looking for cutflow histograms" );

      
      // Now that we have identified a TDirectory in the file (hopefully), load the bin contents from the first
      // bin of each cutflow histogram

      Bool_t success = kTRUE;
      TH1D *h = htools::quiet_load<TH1D>( samples[i].file , tdirname+"/cutflow" , success );
      report::printassert( success , "Failed to find cutflow histogram %s/cutflow" , tdirname.Data() );
      samples[i].nevents_cutflow = h->GetBinContent(1);
      delete h;

      h = htools::quiet_assert_load<TH1D>( samples[i].file , tdirname+"/cutflow_mc" );
      samples[i].nevents_cutflow_mc = h->GetBinContent(1);
      delete h;

      h = htools::quiet_assert_load<TH1D>( samples[i].file , tdirname+"/cutflow_mc_pu" );
      samples[i].nevents_cutflow_mc_pu = h->GetBinContent(1);
      delete h;

      h = htools::quiet_assert_load<TH1D>( samples[i].file , tdirname+"/cutflow_mc_pu_zvtx" );
      samples[i].nevents_cutflow_mc_pu_zvtx = h->GetBinContent(1);
      delete h;

      // Also load the yield from sumWeights tree
      TTree *t = htools::quiet_assert_load<TTree>( samples[i].file , "sumWeights" );
      samples[i].nevents_sumWeights  = 0.0;
      Float_t tmp_nev ( 0 );
      t->SetBranchAddress( "totalEventsWeighted" , &tmp_nev );
      for( Long64_t ientry = 0 ; ientry < t->GetEntries() ; ientry++ ) {
        t->GetEntry( ientry );
        samples[i].nevents_sumWeights += Double_t(tmp_nev);
      }
      delete t;

      // calculate the weight to scale the sample to luminosity
      samples[i].weight = lumi * samples[i].xsec * samples[i].kfactor / samples[i].getNevents();

      samples[i].report( name , lumi );

    }
  }

  
public:

  PhysicsProcess( TString _name , TString _roottitle , TString _latextitle , Int_t _color = -1 , TString _treename = "nominal" )
    : name( _name )
    , roottitle( _roottitle )
    , latextitle( _latextitle )
    , color( _color )
    , treename( _treename )
    , type( UNDEFINED )
  {}

  ~PhysicsProcess() {}

  TString getName() const { return name; }
  TString getRootTitle() const { return roottitle; }
  TString getLatexTitle() const { return latextitle; }
  Int_t getColor() const { return color; }
  Bool_t isData() const { return type==DATA; }
  Bool_t isBackground() const { return type==BACKGROUND; }
  Bool_t isSignal() const { return type==SIGNAL; }
  Bool_t isMC() const { return type==BACKGROUND || type==SIGNAL; }

  void setName( TString _name , TString _roottitle = "" , TString _latextitle = "" ) {
    name = _name;
    if( _roottitle.Length() ) roottitle = _roottitle;
    if( _latextitle.Length() ) latextitle = _latextitle;
  }

  std::size_t numSamples() const { return samples.size(); }
  TString getSampleTag( std::size_t i ) const { return samples[i].tag; }
  Double_t getSampleWeight( std::size_t i ) const { return samples[i].weight; }
  Double_t getSampleXsec( std::size_t i ) const { return samples[i].xsec; }
  Double_t getSampleNevents( std::size_t i ) const { return samples[i].getNevents(); }
  TFile* getSampleFile( std::size_t i ) const { return samples[i].file; }
  TTree* getSampleTree( std::size_t i ) const { return htools::quiet_assert_load<TTree>(samples[i].file,treename); }


  void cleanBlacklisted( std::vector<TString> &blacklist ) {
    for( TString blistTag : blacklist ) {
      for( std::vector<SampleInfo>::iterator is = samples.begin() , fs = samples.end() ; is != fs ; ++is ) {
	if( is->tag == blistTag ) {
	  report::warn( "removing blacklisted sample from %s" , name.Data() );
	  is->print();
	  samples.erase( is );
	  break;
	}
      }
    }
  }

  // Add a ROOT file to this process
  void addSample( TString tag , Double_t xsec = 1.0 , Double_t kfactor = 1.0 , Double_t nevents = 1.0 ) {
    SampleInfo s;
    s.tag			 = tag;
    s.xsec			 = xsec;
    s.kfactor			 = kfactor;
    s.nevents_sumWeights	 = nevents;
    s.nevents_cutflow		 = nevents;
    s.nevents_cutflow_mc	 = nevents;
    s.nevents_cutflow_mc_pu	 = nevents;
    s.nevents_cutflow_mc_pu_zvtx = nevents;
    s.nevents_cutflow_mc_pu_zvtx = nevents;
    samples.push_back( s );
  }

  // Add a ROOT file to this process, look up xsec and kfactor in text file
  void addSampleFindInfo( TString tag , Double_t x = 1.0 , Double_t y = 1.0 ) {
    SampleInfo s;
    s.tag = tag;

    std::ifstream xsecFile;
    xsecFile.open( config::XsecFilePath );
    Bool_t success = kFALSE;
    while( !xsecFile.eof() && !success ) {
      std::string cline;
      getline( xsecFile , cline );
      TString line( cline );
      TObjArray *tokens = line.Tokenize(" \t");
      if( tokens->GetEntries() < 3 ) {
        //report::debug( "line (%i) = %s" , tokens->GetEntries() , line.Data() );
        continue;
      }
      TString tokenTag     = ((TObjString*)tokens->At(0))->GetString().ReplaceAll(" ","").ReplaceAll("\t","");
      TString tokenXsec    = ((TObjString*)tokens->At(1))->GetString().ReplaceAll(" ","").ReplaceAll("\t","");
      TString tokenKfactor = ((TObjString*)tokens->At(2))->GetString().ReplaceAll(" ","").ReplaceAll("\t","");
      if( tag.BeginsWith(tokenTag) ) {
        s.xsec    = tokenXsec.Atof();
        s.kfactor = tokenKfactor.Atof();
        success   = kTRUE;
        break;
      } else if( tokenTag.Contains(tag) ) {
        report::warn( "Possible error: tag=[%s] not matched to token=[%s]" , tag.Data() , tokenTag.Data() );
      } else {
        //report::debug( "line (%i) = %s" , tokens->GetEntries() , line.Data() );
      }
    }
    xsecFile.close();

    if( !success ) {
      report::warn( "Failed to find cross-section / kfactor information for sample with tag = %s" , tag.Data() );
      report::blank( "Setting cross-section to 0!!" );
      s.xsec    = 0;
      s.kfactor = 0;
    }
      
    samples.push_back( s );
  }

  // Load Data file - this function should be used when the sample tag represents full path to the file
  void loadAsData() {
    type = DATA;
    for( SampleInfo s : samples ) {
      s.file = TFile::Open( s.tag );
      s.weight = 1.0; 
      s.report( name );
    }
  }

  // Load Data file - this function should be used when the sample tag just represents part of the full path
  // The full path is looked up automagically
  void loadAsData( TString databasedir , std::vector<TString> &datadirs ) {
    type = DATA;
    TString tmpdir;
    for( std::size_t i = 0 ; i < samples.size() ; i++ ) {
      report::printassert( samples[i].matchTag(datadirs,tmpdir) == 1 , "nmatched for tag %s = %i" , samples[i].tag.Data() , samples[i].matchTag(datadirs,tmpdir) );
      samples[i].file = TFile::Open( databasedir + "/" + tmpdir + "/output.root" );
      samples[i].weight = 1.0;
      samples[i].report( name );
    }
  }

 
  void loadAsBackground( Double_t lumi ) {
    type = BACKGROUND;
    loadAsMC(lumi);
  }

  void loadAsSignal( Double_t lumi ) {
    type = SIGNAL;
    loadAsMC(lumi);
  }
  
  void loadAsBackground( TString databasedir , std::vector<TString> &datadirs , Double_t lumi ) {
    type = BACKGROUND;
    loadAsMC(databasedir,datadirs,lumi);
  }

  void loadAsSignal( TString databasedir , std::vector<TString> &datadirs , Double_t lumi ) {
    type = SIGNAL;
    loadAsMC(databasedir,datadirs,lumi);
  }

  
  void dumpSamples( TString databasedir , Bool_t inclXsec = kTRUE ) {
    for( SampleInfo s : samples ) {
      if( inclXsec ) {
        report::blank( "%-80s %12g %12g %12g %12g" ,
		       s.file->GetName() ,
		       s.getNevents() ,
		       s.nevents_sumWeights ,
		       s.nevents_cutflow_mc ,
		       s.nevents_cutflow_mc_pu );
      } else {
        report::blank( "%-80s" , s.file->GetName() );
      }
    }
  }


};

#endif
