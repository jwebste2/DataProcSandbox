#ifndef _INPUTCONFIG_H_
#define _INPUTCONFIG_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <map>

#include "TString.h"
#include "TSystem.h"
#include "TLegend.h"

#include "Report.h"
#include "HistTools.h"
#include "PhysicsProcess.h"
#include "Config.h"


typedef std::map<TString,Double_t>::iterator it_lumiblocks;

class
InputConfig
{

private:

  TString databasedir;
  Double_t lumi;

  Bool_t sort_backgrounds;
  
  std::vector<PhysicsProcess> signals;
  std::vector<PhysicsProcess> backgrounds;
  std::vector<PhysicsProcess> data;

  std::vector<TString> blacklist;
  
  std::map<TString,Double_t> lumiblocks; // lumis in pb here taken from lumi calculator

public:

  InputConfig( TString _databasedir = "null" )
    : databasedir( _databasedir )
    , lumi( 0.0 )
    , sort_backgrounds( kTRUE )
  {

    // expand the base directory for input data
    gSystem->ExpandPathName(databasedir);

    // 50 ns
    lumiblocks["267073"] = 0.086042924;
    lumiblocks["267167"] = 0.2251302;
    lumiblocks["267358"] = 0.000081767;
    lumiblocks["267359"] = 0.000243973;
    lumiblocks["267360"] = 0.000324914;
    lumiblocks["267367"] = 0.00154844;
    lumiblocks["267385"] = 0.00600499;
    lumiblocks["267599"] = 0.005677128;
    lumiblocks["267638"] = 3.3730915;
    lumiblocks["267639"] = 2.9827723;
    lumiblocks["270806"] = 0.775558;
    lumiblocks["270953"] = 4.586361;
    lumiblocks["271048"] = 6.734014;
    lumiblocks["271298"] = 9.97045733;
    lumiblocks["271421"] = 13.6103267;
    lumiblocks["271516"] = 19.72344;
    lumiblocks["271595"] = 16.59377;
    lumiblocks["271744"] = 6.29282;

    // 25 ns
    lumiblocks["276262"] = 6.496017;
    lumiblocks["276329"] = 12.375005;
    lumiblocks["276336"] = 0.689043;
    lumiblocks["276416"] = 4.38556;
    lumiblocks["276511"] = 9.354875;
    lumiblocks["276689"] = 11.281437;
    lumiblocks["276731"] = 25.74843;
    lumiblocks["276778"] = 0.7012;
    lumiblocks["276790"] = 1.65671;
    lumiblocks["276952"] = 6.90003;
    lumiblocks["276954"] = 0.770813;
    lumiblocks["278880"] = 22.67796;
    lumiblocks["278912"] = 20.741551;
    lumiblocks["278968"] = 10.342965;
    lumiblocks["279169"] = 57.773407;
    lumiblocks["279259"] = 7.61594;
    lumiblocks["279279"] = 17.793527;
    lumiblocks["279284"] = 33.3266;
    lumiblocks["279345"] = 54.919944;
    lumiblocks["279515"] = 0.362112;
    lumiblocks["279598"] = 71.405581;
    lumiblocks["279685"] = 78.247129;
    lumiblocks["279764"] = 2.492831;
    lumiblocks["279813"] = 50.3832;
    lumiblocks["279867"] = 31.25825;
    lumiblocks["279928"] = 1.063;
    lumiblocks["279932"] = 46.91101;
    lumiblocks["279984"] = 68.746628;
    lumiblocks["280231"] = 93.249766;
    lumiblocks["280319"] = 98.60357;
    lumiblocks["280368"] = 8.739077;
    lumiblocks["280423"] = 72.32858;
    lumiblocks["280464"] = 61.505421;
    lumiblocks["280500"] = 7.7279;
    lumiblocks["280520"] = 12.880305;
    lumiblocks["280614"] = 25.77352;
    
  }
  ~InputConfig() {}

  TString getDataBaseDir() { return databasedir; }
  std::size_t numSignals() { return signals.size(); }
  std::size_t numBackgrounds() { return backgrounds.size(); }
  std::size_t numData() { return data.size(); }
  std::size_t numProcesses() { return signals.size() + backgrounds.size() + data.size(); }
  PhysicsProcess getSignal( std::size_t i ) const { return signals[i]; }
  PhysicsProcess getBackground( std::size_t i ) const { return backgrounds[i]; }
  PhysicsProcess getData( std::size_t i = 0 ) const { return data[i]; }
  PhysicsProcess getProcess( std::size_t i = 0 ) const {
    if( i < signals.size() ) return signals[i];
    if( i < signals.size() + backgrounds.size() ) return backgrounds[i-signals.size()];
    return data[i-signals.size()-backgrounds.size()];
  }

  Bool_t getSortBackgrounds() { return sort_backgrounds; }

  void setSortBackgrounds( Bool_t v ) { sort_backgrounds = v; }
  void addSignalProcess( PhysicsProcess p ) { signals.push_back(p); }
  void addBackgroundProcess( PhysicsProcess p ) { backgrounds.push_back(p); }
   void addDataProcess( PhysicsProcess p ) { data.push_back(p); }

  void addSignalProcess( TString _name , TString _roottitle , TString _latextitle , Int_t _color ) { signals.push_back( PhysicsProcess(_name,_roottitle,_latextitle,_color) ); }
  void addBackgroundProcess( TString _name , TString _roottitle , TString _latextitle , Int_t _color ) { backgrounds.push_back( PhysicsProcess(_name,_roottitle,_latextitle,_color) ); }
  void addDataProcess( TString _name , TString _roottitle , TString _latextitle ) { data.push_back( PhysicsProcess(_name,_roottitle,_latextitle) ); }

  void addSignalSample( TString _tag , Double_t _xsec , Double_t _kfactor ) { signals.back().addSample(_tag,_xsec,_kfactor); }
  void addBackgroundSample( TString _tag , Double_t _xsec , Double_t _kfactor ) { backgrounds.back().addSample(_tag,_xsec,_kfactor); }
  void addDataSample( TString _tag ) { data.back().addSample(_tag); }

  void addToBlacklist( TString _tag ) { blacklist.push_back(_tag); }
  void setBlacklist( std::vector<TString> v ) { blacklist = v; }
  
  void setLuminosity( Double_t _lumi ) { lumi = _lumi; }
  void setLuminosityFromDataTags() {
    assert( data.size() == 1 );
    lumi = 0.0;
    for( std::size_t i = 0 ; i < data[0].numSamples() ; i++ ) {
      // skip the blacklisted datasets
      if( std::find(blacklist.begin(),blacklist.end(),data[0].getSampleTag(i)) != blacklist.end() ) {
	//report::warn( "skipping blacklisted %s in lumi calculation" , data[0].getSampleTag(i).Data() );
	continue;
      }
      // loop over luminosity blocks until we find the one that matches this sample tag
      for( it_lumiblocks j = lumiblocks.begin() ; j != lumiblocks.end() ; j++ ) {
	if( data[0].getSampleTag(i).Contains(j->first) ) {
	  lumi += j->second;
	  break;
	}
      }
      //lumi += lumiblocks[data[0].getSampleTag(i)] / 1000.0;
    }
    report::info( "calculated Lumi = %g for InputConfig %s" , lumi , databasedir.Data() );
  }


  Double_t getLuminosity() { return lumi; }

  
  //
  // Fancy load function for tth analysis
  // Looks for directories with the correct tag containing "output.root" files
  //
  void load() {

    // get list of data directories in this base path
    std::vector<TString> datadirs;
    {
      void* dirp = gSystem->OpenDirectory(databasedir);
      
      const char* entry;
      Int_t n = 0;
      TString tmpdir;
      
      while( (entry = (char*)gSystem->GetDirEntry(dirp)) ) {
	tmpdir = entry;
	if( tmpdir.Contains("user.") ) datadirs.push_back( tmpdir );
      }

      delete entry;
    }

    for( std::size_t i = 0 ; i < signals.size() ; i++ ) { signals[i].cleanBlacklisted(blacklist); signals[i].loadAsSignal(databasedir,datadirs,lumi); }
    for( std::size_t i = 0 ; i < backgrounds.size() ; i++ ) { backgrounds[i].cleanBlacklisted(blacklist); backgrounds[i].loadAsBackground(databasedir,datadirs,lumi); }
    for( std::size_t i = 0 ; i < data.size() ; i++ ) { data[i].cleanBlacklisted(blacklist); data[i].loadAsData(databasedir,datadirs); }

  }

 

  //
  // Simple load function that just loads using filenames
  //
  void loadRootFiles() {
    for( std::size_t i = 0 ; i < signals.size() ; i++ ) signals[i].loadAsSignal(lumi);
    for( std::size_t i = 0 ; i < backgrounds.size() ; i++ ) backgrounds[i].loadAsBackground(lumi);
    for( std::size_t i = 0 ; i < data.size() ; i++ ) data[i].loadAsData();
  }

  
  void dumpSamples() {
    for( PhysicsProcess proc : signals ) proc.dumpSamples( databasedir );
    for( PhysicsProcess proc : backgrounds ) proc.dumpSamples( databasedir );
    for( PhysicsProcess proc : data ) proc.dumpSamples( databasedir , kFALSE );
  }

};

#endif

