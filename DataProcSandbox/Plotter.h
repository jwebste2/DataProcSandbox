#ifndef _PLOTTER_H_
#define _PLOTTER_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>
#include <TProfile.h>
#include <TArrow.h>
#include <TNamed.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooNDKeysPdf.h>
#include <RooBinning.h>

#include "Report.h"
#include "ArgParser.h"
#include "HistTools.h"


class
Plotter
{
 private:

  TCanvas *canvas;
  TFile *rootfile;

  TString title;
  TString psname;

  Int_t plotid;
  Int_t histid;

 public:

  TCanvas* getCanvas() { return canvas; }
  void updateCanvas() { canvas->Update(); }
  TString getPsName() { return psname; }

  TString getRootPath() { return TString::Format("output/%s.root",psname.Data()); }
  TString getPsPath() { return TString::Format("output/%s.ps",psname.Data()); }
  TString getEpsDir() { return TString::Format("output/%s_eps/",psname.Data()); }
  TString getPngDir() { return TString::Format("output/%s_png/",psname.Data()); }
  TString getPdfDir() { return TString::Format("output/%s_pdf/",psname.Data()); }  
  TString getCDir() { return TString::Format("output/%s_C/",psname.Data()); }

  // open a .ps file for saving plots
  TCanvas* openPs() { 
    canvas->Print( TString::Format("%s[",getPsPath().Data()) ); 
    return canvas;
  }

  void openDir() {
    // if directory already exists then clear it out
    std::system( TString::Format("rm -rf %s; mkdir -p %s",getEpsDir().Data(),getEpsDir().Data()).Data() );
    std::system( TString::Format("rm -rf %s; mkdir -p %s",getPngDir().Data(),getPngDir().Data()).Data() );
    std::system( TString::Format("rm -rf %s; mkdir -p %s",getPdfDir().Data(),getPdfDir().Data()).Data() );
    std::system( TString::Format("rm -rf %s; mkdir -p %s",getCDir().Data(),getCDir().Data()).Data() );
    // set the plot ID to 1 so we know to save plots
    plotid = 1;
  }
  void closeDir() { plotid = -1; }

  // open .root file for saving root objects
  TFile* openRoot() {
    rootfile = new TFile( getRootPath() , "RECREATE" );
    return rootfile;
  }

  // print canvas to .ps file
  Int_t print() { 
    canvas->Print( getPsPath() ); 
    if( plotid > 0 ) {
      canvas->SaveAs( TString::Format("%s/%.3i.eps",getEpsDir().Data(),plotid) );
      canvas->SaveAs( TString::Format("%s/%.3i.png",getPngDir().Data(),plotid) );
      canvas->SaveAs( TString::Format("%s/%.3i.pdf",getPdfDir().Data(),plotid) );
      canvas->SaveAs( TString::Format("%s/%.3i.C",getCDir().Data(),plotid) );
      return plotid++;
    }
    return plotid;
  }

  // print some text to .ps file
  void print( TString text1 , TString text2 , TString text3 = "" , TString text4 = "" , TString text5 = "" ) {
    canvas->Clear();
    TPaveText lt( 0. , 0. , 1. , 1. , "NDC" );
    lt.SetTextAlign( 22 );
    lt.AddText( text1 );
    if( text2.Length() ) lt.AddText( text2 );
    if( text3.Length() ) lt.AddText( text3 );
    if( text4.Length() ) lt.AddText( text4 );
    if( text5.Length() ) lt.AddText( text5 );
    lt.Draw();
    print();
  }

  // print string  to .ps file
  void print( std::string str ) {
    canvas->Clear();
    TPaveText lt( 0. , 0. , 1. , 1. , "NDC" );
    lt.SetTextAlign( 12 );
    std::istringstream iss( str );
    std::string line;
    while( std::getline( iss , line ) ) {
      std::cout << "  " << line << std::endl;
      lt.AddText( line.c_str() );
    }
    lt.Draw();
    print();
  }

  // close the .ps file and dump some useful commands for downloading and viewing
  void closePs() { 
    canvas->Print( TString::Format("%s]",getPsPath().Data()) ); 
    std::system( TString::Format( "echo \"  rm -rf %s_C; anl-get ${PWD}/%s\"" , psname.Data() , getCDir().Data() ).Data() );
    std::system( TString::Format( "echo \"  rm -rf %s_eps; anl-get ${PWD}/%s; open %s_eps\"" , psname.Data() , getEpsDir().Data() , psname.Data() ).Data() );
    std::system( TString::Format( "echo \"  rm -rf %s_png; anl-get ${PWD}/%s; open %s_png\"" , psname.Data() , getPngDir().Data() , psname.Data() ).Data() );
    std::system( TString::Format( "echo \"  rm -rf %s_pdf; anl-get ${PWD}/%s; open %s_pdf\"" , psname.Data() , getPdfDir().Data() , psname.Data() ).Data() );
    std::system( TString::Format( "echo \"  anl-get ${PWD}/%s; preview %s.ps\"" , getPsPath().Data() , psname.Data() ).Data() );
    //std::system( TString::Format( "open -a Preview %s" , getPsPath().Data() ) );
  }

  void closeRoot() {
    rootfile->Close();
    delete rootfile;
    std::system( TString::Format( "echo \"  rm -rf %s.root; anl-get %s .; root -l %s.root\"" , psname.Data() , getRootPath().Data() , psname.Data() ).Data() );
  }

  // save oss to txt file
  /*
  void saveTxt( std::ostringstream *oss ) {
    ofstream sf;
    sf.open( "output/" + psname + ".txt" );
    sf << oss->str();
    sf.close();
  }
  */
  
  template<typename valT> void divide( valT n ) { canvas->DivideSquare( Int_t(n) , 0.0 , 0.0 ); }
  template<typename valT> void cd( valT n ) { canvas->cd( Int_t(n) ); }



 public:
  
  Plotter( ArgParser ap ) 
    : title( ap.getTitle() )
    , psname( ap.getTag() )
    , plotid( -1 )
    , histid( 0 )
  {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    //canvas = new TCanvas( "canvas" , "canvas" , 800 , 600 );
    canvas = new TCanvas( "canvas" , "canvas" , 800 , 800 );
    htools::setup( canvas );
    htools::setup( gStyle );
    gErrorIgnoreLevel = kWarning; // turn off PostScript messages that appear every time we draw a plot
  }

  ~Plotter() {
    delete canvas;
  }

  void giveSimpleName( TNamed *n ) { n->SetName( TString::Format( "named_%i" , histid++ ) ); }
  
  void setCanvas1D() { htools::setCanvas1D( canvas ); }
  void setCanvas2D() { htools::setCanvas2D( canvas ); }
  void setCanvasPie() { htools::setCanvasPie( canvas ); }
  void setCanvasHorizontal() { htools::setCanvasHorizontal( canvas ); }
  void setCanvasLog( Bool_t logx=kFALSE , Bool_t logy=kFALSE ) { htools::setCanvasLog( canvas , logx , logy ); }

  template< typename RootObj > void write( RootObj *obj , const char* name = 0 ) { 
    assert( rootfile->IsOpen() );
    obj->SetDirectory( rootfile );
    obj->Write(name); 
  }

};

#endif



