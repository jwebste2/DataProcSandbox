#ifndef _ARGPARSER_H_
#define _ARGPARSER_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <map>

#include "TString.h"

#include "Report.h"

class
ArgParser
{ 
  
private:

  struct InputArgument {

    TString key;
    TString value;
    TString description;
    TString options;
    TString defaultval;
    Bool_t optional;

    InputArgument( const TString &_key , const TString &_description , const TString &_options , const TString &_defaultval = "" , const Bool_t &_optional = kFALSE ) {

      key = _key;
      value = "NULL";
      description = _description;
      options = _options;
      defaultval = _defaultval;
      optional = _optional;

      if( description.Length()==0 ) {
	description = "No description provided";
      }

    }

    TString getLabel() {
      if( defaultval.Length()==0 ) return key;
      return key+"[default="+defaultval+"]";
    }

    Bool_t required() { return !optional; }
    
  };


  std::vector<InputArgument> argvec;
  TString title;
  TString descr;

  TString tag; // argument keys separated by underscores, used for filename
  TString tagtail; // just the keys separated by underscores
  TString runcommand; // argument keys separated by spaces, used for example function call
  

public:

  ArgParser( const TString &_title , const TString &_descr = "" ) 
    : title( _title )
    , descr( _descr )
  { std::system( "hostname" ); }
  ~ArgParser() {}

  void addArg( const TString &key , const TString &description = "" , const TString &options = "" ) { argvec.push_back( InputArgument(key,description,options) ); }

  void addOptionalArg( const TString &key , const TString &description , const TString &defaultval ) { argvec.push_back( InputArgument(key,description,"DEFAULT="+defaultval,defaultval,kTRUE) ); }

  void parse( int argc , char* argv[] ) {

    // print description of function
    if( descr.Length() ) report::info( descr );
    
    // Define runcommand
    // Also count the number of arguments that need to be supplied by the user
    runcommand = title;
    int nargs_required = 0;
    bool last_required = true;
    for( std::size_t ik = 0 ; ik < argvec.size() ; ik++ ) {
      runcommand += (" "+argvec[ik].getLabel());
      bool this_required = argvec[ik].required();
      if( this_required ) {
	report::printassert( last_required , "Optional arguments must follow the required arguments" );
	nargs_required++;
      }
      last_required = this_required;
    }

    // assert correct number of arguments
    if( (argc-1) < nargs_required ) {
      report::error( "Command entered differs from what is expected" );
      report::info( "Please use the following format:" );
      report::spaceprint( "./run/%s" , runcommand.Data() );
      //report::blank( "Argument options..." );
      report::blank( "%-20s %-70s %-80s" , "ARG" , "DESCRIPTION" , "EXAMPLES" );
      report::blank( "%-20s %-70s %-80s" , "===" , "===========" , "========" );
      for( std::size_t ik = 0 ; ik < argvec.size() ; ik++ ) {
	report::blank( "%-20s %-70s %-80s" , argvec[ik].key.Data() , argvec[ik].description.Data() , argvec[ik].options.Data() );
      }
      assert( false );
    }

    tag = title;
    tagtail = "";
    std::size_t ik = 0;
    int argv_iter = 1;
    for( ; argv_iter < argc ; argv_iter++ ) {
      argvec[ik].value = TString( argv[argv_iter] );
      tagtail += ("_"+argvec[ik].value);
      report::info( "Parsed  Input #%2i : %-20s = %s" , argv_iter , argvec[ik].key.Data() , argvec[ik].value.Data() );
      ik++;
    }
    
    while( ik < argvec.size() ) {
      argvec[ik].value = argvec[ik].defaultval;
      tagtail += ("_"+argvec[ik].value);
      report::info( "Default Input #%2i : %-20s = %s" , argv_iter , argvec[ik].key.Data() , argvec[ik].value.Data() );      
      ik++;
    }

    tagtail.ReplaceAll( ".." , "-" );
    tagtail.ReplaceAll( "/" , "+" );
    tag = title+tagtail;
    report::info( "Parsed runtag = %s" , tag.Data() );

  }

  TString get( const TString &key ) { 
    for( std::size_t ik = 0 ; ik < argvec.size() ; ik++ ) {
      if( argvec[ik].key==key ) return argvec[ik].value;
    }
    report::printassert( "Could not locate key %s" , key.Data() );
    return "";
  }

  double getAtof( const TString &key ) { return get(key).Atof(); }
  int getAtoi( const TString &key ) { return get(key).Atoi(); }

  TString operator[]( const TString &key ) { return get(key); }

  TString getTitle() { return title; }
  TString getTag() { return tag; }
  TString getTagTail() { return tagtail; }
  
  

};

#endif
