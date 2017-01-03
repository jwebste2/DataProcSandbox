#ifndef _REPORT_H_
#define _REPORT_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>

#include "TString.h"
#include "TStopwatch.h"
#include "TRandom3.h"
#include "TMath.h"

namespace
report
{ 

  static bool suppress_output = false;
  const int dbg_lvl	      = 1;

  const int progressBarWidth    = 80;
  TString progressBarName       = "";
  TString progressBarBlankChar  = "░";
  TString progressBarFilledChar = "█";
  Double_t progressFraction     = 0;
  TStopwatch progressStopwatch;
  
  static void vreport( const TString tag , const char *fmt , va_list ap , bool incl_space = true ) {
    if( suppress_output ) return;
    std::cout << tag;
    vprintf( fmt , ap );
    if( incl_space ) std::cout << std::endl;
  }

  static void vreport( const TString tag , const TString str ) {
    if( suppress_output ) return;
    std::cout << tag << str << std::endl;
  }
  
  static void vreport_color( TString color , const TString tag , const char *fmt , va_list ap , bool incl_space = true ) {
    if( suppress_output ) return;
    // coloring doesn't work on local machine, so overriding
    //vreport( tag , fmt , ap , incl_space );
    //return;
    TString echocommand = TString::Format( "echo -e -n \"${%s}%s${FONTNONE}\"" , color.Data() , tag.Data() );
    std::system( echocommand.Data() );
    vprintf( fmt , ap );
    if( incl_space ) std::cout << std::endl;
  }

  static void vreport_color( TString color , const TString tag , const TString str ) {
    if( suppress_output ) return;
    // coloring doesn't work on local machine, so overriding
    vreport_color( color , tag , str );
  }

  static void vreport( TStopwatch &tw , const TString tag , const char *fmt , va_list ap ) {
    if( suppress_output ) return;
    vreport( tag , fmt , ap , kFALSE );
    std::cout << " -- ";
    tw.Print();
  }
  
  static void info( const char *fmt , ... )  { va_list ap; va_start(ap,fmt); vreport("INFO  : ",fmt,ap); }
  static void info( TStopwatch &tw , const char *fmt , ... )  { va_list ap; va_start(ap,fmt); vreport(tw,"INFO  : ",fmt,ap); }
  static void exam( const char *fmt , ... )  { va_list ap; va_start(ap,fmt); vreport_color("FONTPURPLE","EXAM  : ",fmt,ap); }
  static void warn( const char *fmt , ... )  { va_list ap; va_start(ap,fmt); vreport_color("FONTYELLOW","WARN  : ",fmt,ap); }
  static void error( const char *fmt , ... ) { va_list ap; va_start(ap,fmt); vreport_color("FONTRED","ERROR : ",fmt,ap); }
  static void blank( const char *fmt , ... ) { va_list ap; va_start(ap,fmt); vreport("  ",fmt,ap); }
  static void debug( const char *fmt , ... ) { va_list ap; va_start(ap,fmt); vreport_color("FONTCYAN","DEBUG : ",fmt,ap); }
  static void line() { std::cout << std::endl; }

  static void spaceprint( const char *fmt , ... ) { va_list ap; va_start(ap,fmt); vreport("      > ",fmt,ap); }
  static void spaceprint( TObject *o , TString option = "a" ) {
    std::cout << "      > ";
    o->Print(option.Data());
  }

  static void time( TStopwatch &tw ) { std::cout << "TIME  : "; tw.Print(); }
  
  static void debug( int lvl , const char *fmt , ... ) {
    if( lvl > dbg_lvl ) return;
    va_list ap; va_start(ap,fmt); vreport("DEBUG : ",fmt,ap); 
  }

  static void warnif( bool cond , const char *fmt , ... ) {
    if( ! cond ) return;
    va_list ap; va_start(ap,fmt); vreport_color("FONTYELLOW","WARN  : ",fmt,ap);
  }

  static void printassert( bool cond , const char* fmt , ... ) {
    if( cond ) return;
    va_list ap; va_start(ap,fmt); vreport("ASSERT FAILED : ",fmt,ap);
    assert( false );
  }

  static void title( const char *_info , const char *_exam ) {
    info( _info );
    exam( _exam );
  }

  static void startProgressBar( TString name = "" ) {
    progressFraction = 0;
    progressBarName = "  "+name;
    std::cout << progressBarName << " |";
    for( int i = 0 ; i < progressBarWidth ; ++i ) std::cout << progressBarBlankChar;
    std::cout << "| Starting\r";
    progressStopwatch.Start();
  }

  static void endProgressBar() {
    std::cout << std::endl;
  }

  static void updateProgressBar( Double_t finished , Double_t total = 1.0 , TString extrainfo = "" ) {

    if( finished==total ) {
      Double_t currTime = progressStopwatch.RealTime();
      std::cout << progressBarName << " |";
      for( int i = 0 ; i < progressBarWidth ; ++i ) std::cout << progressBarFilledChar;
      std::cout << "| nev=" << finished << "/" << total << TString::Format(" t=%0.0f/%0.0f r=%0.1fkHz",currTime,total*currTime/finished,finished/currTime/1000.0).Data() << extrainfo << "\r";
      endProgressBar();
      return;
    }

    if( (finished / total) - progressFraction < 0.001 ) return;
    progressFraction = finished / total;

    Double_t currTime = progressStopwatch.RealTime();
    progressStopwatch.Start(kFALSE);
    std::cout.flush();
    std::cout << progressBarName << " |";
    int pos = progressBarWidth * (finished / total);
    if( total == 0 ) pos = progressBarWidth;
    for( int i = 0 ; i < progressBarWidth ; ++i ) {
      if( i <= pos ) std::cout << progressBarFilledChar;
      //else if( i <= TMath::FloorNint(TMath::Power(pos,1.1)) ) std::cout << "▓";
      //else if( i <= TMath::FloorNint(TMath::Power(pos,1.5)) ) std::cout << "▒";
      //else if( i <= TMath::FloorNint(TMath::Power(pos,2.0)) ) std::cout << "░";
      //else if( i <= TMath::FloorNint(TMath::Power(pos,2.5)) ) std::cout << "ː";
      else std::cout << progressBarBlankChar;
    }
    std::cout << "| nev=" << finished << "/" << total << TString::Format(" t=%0.0f/%0.0f r=%0.1fkHz",currTime,total*currTime/finished,finished/currTime/1000.0).Data() << extrainfo << "\r";

  }

  
  /*
   * Function for keeping track of memory usage
   * Doesn't work without dump_mem.sh script
   *
  static void memory( TString jobname ) {
    TString script = tfstring::Format( ". scripts/dump_mem.sh %s" , jobname.Data() );
    system( script.Data() );
    std::cout << std::endl;
  }
  */

  
};

#endif
