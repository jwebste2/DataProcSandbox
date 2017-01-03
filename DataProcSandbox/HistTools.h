#ifndef _HISTTOOLS_H_
#define _HISTTOOLS_H_

#include <iostream>
#include <vector>

#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TEfficiency.h>
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
#include <TIterator.h>
#include <TRandom3.h>
#include <TPie.h>

#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooWorkspace.h>
#include <RooCmdArg.h>
#include <RooLinkedList.h>
#include <RooDataHist.h>

#include "Report.h"
#include "ArgParser.h"

namespace
htools
{

  const Double_t GenericTitleSize( 0.04 );
  const Double_t GenericLabelSize( 0.03 );
  const Double_t GenericTitleOffset( 1.4 );
  const Style_t GenericFont( 42 );

  // ATLAS Tags
  const TString AtlasInternal	 = "#bf{#it{ATLAS}} Internal";
  const TString AtlasPreliminary = "#bf{#it{ATLAS}} Preliminary";
  const TString Atlas		 = "#bf{#it{ATLAS}}";
  const TString AtlasDefault	 = AtlasInternal;

  std::vector<TObject*> garbage_collector_TObject;
  
  enum HistStyle {
    HOLLOW = 0 ,
    SOLID  = 1 ,
    DASHED = 3
  };

  enum HorizontalPosition {
    LEFT ,
    CENTER ,
    CENTERRIGHT ,
    RIGHTCENTER ,
    RIGHT
  };

  enum VerticalPosition {
    ABOVE ,
    TOP ,
    BOTTOM
  };
  

  static void clear() {
    for( TObject *o : garbage_collector_TObject ) if( o ) delete o;
    garbage_collector_TObject.clear();
  }
  

  // 
  // FUNCTIONS FOR LOADING ROOT OBJECTS
  //

  template< typename HistT > static HistT* load( TFile *f , TString name , Bool_t &success , Bool_t quiet , TString newname = "" ) {
    if( !f ) {
      if( !quiet ) report::warn( "file with hists not instantiated correctly" );
      success = kFALSE;
      return (HistT*) 0;
    }
    if( ! f->IsOpen() ) { 
      if( !quiet ) report::warn( "file with hists is not open" );
      success = kFALSE;
      return (HistT*) 0;
    }
    HistT *tmp = (HistT*) f->Get(name);
    if( !tmp ) {
      if( !quiet ) report::warn( TString::Format("missing hist %s in file %s",name.Data(),f->GetName()) );
      success = kFALSE;
      return (HistT*) 0;
    }
    //tmp->SetDirectory( 0 );
    if( newname.Length() ) tmp->SetName( newname );
    if( !quiet ) report::info( TString::Format("successfully loaded object %s from file %s",name.Data(),f->GetName()) );
    return tmp;
  }

  template< typename HistT > static HistT* load( TFile *f , TString name ) {
    Bool_t success = kTRUE;
    return load<HistT>( f , name , success , kFALSE , "" );
  }
  
  template< typename HistT > static HistT* load( TFile *f , TString name , Bool_t &success , TString newname = "" ) {
    return load<HistT>( f , name , success , kFALSE , newname );
  }

  template< typename HistT > static HistT* quiet_load( TFile *f , TString name , Bool_t &success , TString newname = "" ) {
    return load<HistT>( f , name , success , kTRUE , newname );
  }

  template< typename HistT > static HistT* assert_load( TFile *f , TString name ) {
    Bool_t success = kTRUE;
    HistT *h = load<HistT>( f , name , success );
    assert( success );
    return h;
  }

  template< typename HistT > static HistT* quiet_assert_load( TFile *f , TString name ) {
    Bool_t success = kTRUE;
    HistT *h = quiet_load<HistT>( f , name , success );
    assert( success );
    return h;
  }

  template< typename HistT > static HistT* createWeightHist( TString name , HistT *num , HistT *den ) {
    HistT *w = (HistT*) num->Clone( name );
    w->Divide( den );
    return w;
  }  

  // 
  // SAFE CLEANUP FUNCTIONS
  //
  
  template< typename HistT >
  static void safe_delete( HistT *h ) { if(h) delete h; }

  template< typename HistT >
  static void safe_delete( std::vector<HistT*> &v ) { for( HistT* h : v ) safe_delete(h); }

  //
  // SAFE HISTOGRAM ADD FUNCTIONS
  //

  template< typename HistT >
  static void safe_add( HistT *hsum , HistT *h , const TString &clonetitle ) {
    if( !hsum ) hsum = (HistT*) h->Clone(clonetitle);
    else hsum->Add( h );
  }

  
  //
  // LOADING WORKSPACES
  //

  static RooWorkspace* loadWorkspace( TString path = "ws_reco_factorized.root" , TString wsname = "ws_reco_factorized" ) {
    report::info( "attempting to load workspace %s from file %s" , wsname.Data() , path.Data() );
    TFile *wsfile = TFile::Open( path );
    RooWorkspace *ws = (RooWorkspace*) wsfile->Get( wsname );
    assert( ws );
    return ws;
  }

  // 
  // LOADING SPECIFIC DATASETS
  //


  // currently broken
  static RooLinkedList* loadFitCmds( RooWorkspace *ws , TString channel ) {
    RooLinkedList *fitcmds = new RooLinkedList();
    fitcmds->Add( new RooCmdArg( RooFit::CloneData(kFALSE) ) );
    fitcmds->Add( new RooCmdArg( RooFit::Extended(kTRUE) ) );
    fitcmds->Add( new RooCmdArg( RooFit::Constrain(*ws->set("constrained_"+channel)) ) );
    fitcmds->Add( new RooCmdArg( RooFit::Strategy(0) ) );
    fitcmds->Add( new RooCmdArg( RooFit::PrintLevel(-1) ) );
    fitcmds->Add( new RooCmdArg( RooFit::InitialHesse(kFALSE) ) );
    fitcmds->Add( new RooCmdArg( RooFit::Hesse(kFALSE) ) );
    fitcmds->Add( new RooCmdArg( RooFit::Minos(kFALSE) ) );
    fitcmds->Add( new RooCmdArg( RooFit::SumW2Error(kTRUE) ) );
    return fitcmds;
  }


  // 
  // HIST DRAWING/SETUP TOOLS
  //
  
  static short getColor( Int_t icolor ) {
    if( icolor==-1 ) return kBlack;
    if( icolor==0 ) return kBlue-6;
    if( icolor==1 ) return kRed-6;
    if( icolor==2 ) return kGreen-6;
    if( icolor==3 ) return kMagenta-6;
    if( icolor==4 ) return kCyan-6;
    if( icolor==5 ) return kYellow-6;
    if( icolor==6 ) return kBlue+3;
    if( icolor==7 ) return kRed+3;
    if( icolor==8 ) return kGreen+3;
    if( icolor==9 ) return kMagenta+3;
    if( icolor==10 ) return kCyan+3;
    if( icolor==11 ) return kYellow+3;
    if( icolor==68 ) return kGreen;
    if( icolor==95 ) return kYellow;
    if( icolor >= 100 ) return icolor-100;
    return kBlack;
  }

  static short getGray( std::size_t icolor ) {
    if( icolor==0 ) return 17;
    if( icolor==1 ) return 16;
    if( icolor==2 ) return 15;
    if( icolor==3 ) return 14;
    if( icolor==4 ) return 13;
    if( icolor==5 ) return 12;
    return kBlack;
  }

  template< typename HistT > static void normalize( HistT* h , Double_t n = 1.0 ) { h->Scale( n / h->Integral() ); }

  static void scaleMax( TH1 *h , Double_t n ) { h->SetMaximum( h->GetMaximum() * n ); }

  static void gaussianFlux( TH1 *h , TRandom3 &rand ) {
    // introduce Gaussian bin-by-bin fluctuations based on bin uncertainty
    Double_t entries( h->GetEntries() );
    Double_t sumweights( h->Integral() );
    assert( sumweights > 0 );
    for( Int_t ibx = 1 ; ibx <= h->GetNbinsX() ; ibx++ ) {
      Double_t old_content( h->GetBinContent(ibx) );
      Double_t old_error( h->GetBinError(ibx) );
      Double_t new_content = rand.Gaus( old_content , old_error );
      //report::debug( "%20s : rescaling bin %5i from %12g +/- %12g to %12g" , h->GetName() , ibx , old_content , old_error , new_content );
      h->SetBinContent( ibx , new_content );
    }
  }
  static void gaussianFlux( TH2 *h , TRandom3 &rand ) {
    Double_t entries( h->GetEntries() );
    Double_t sumweights( h->Integral() );
    assert( sumweights > 0 );
    for( Int_t ibx = 1 ; ibx <= h->GetNbinsX() ; ibx++ ) {
      for( Int_t iby = 1 ; iby <= h->GetNbinsY() ; iby++ ) {
	Double_t old_content( h->GetBinContent(ibx,iby) );
	Double_t old_error( h->GetBinError(ibx,iby) );
	Double_t new_content = rand.Gaus( old_content , old_error );
	//report::debug( "%20s : rescaling bin %5i %5i from %12g +/- %12g to %12g" , h->GetName() , ibx , iby , old_content , old_error , new_content );
	h->SetBinContent( ibx , iby , new_content );
      }
    }
  }
  static void gaussianFlux( TH3 *h , TRandom3 &rand ) {
    Double_t entries( h->GetEntries() );
    Double_t sumweights( h->Integral() );
    assert( sumweights > 0 );
    for( Int_t ibx = 1 ; ibx <= h->GetNbinsX() ; ibx++ ) {
      for( Int_t iby = 1 ; iby <= h->GetNbinsY() ; iby++ ) {
	for( Int_t ibz = 1 ; ibz <= h->GetNbinsZ() ; ibz++ ) {
	  Double_t old_content( h->GetBinContent(ibx,iby,ibz) );
	  Double_t old_error( h->GetBinError(ibx,iby,ibz) );
	  Double_t new_content = rand.Gaus( old_content , old_error );
	  //report::debug( "%20s : rescaling bin %5i %5i from %12g +/- %12g to %12g" , h->GetName() , ibx , iby , old_content , old_error , new_content );
	  h->SetBinContent( ibx , iby , ibz , new_content );
	}
      }
    }
  }


  static void setup( TLatex *l , Double_t textscale = 1.0 , std::size_t color = -1 ) {
    Double_t padheight = 1.;
    l->SetTextFont( 42 );
    l->SetTextColor( getColor(color) );
    Double_t textsize = textscale * 0.04 / padheight;
    l->SetTextSize( textsize );
  }

  template<typename TextObj> static void setupHorizontal( TextObj *p , HorizontalPosition hp , Double_t textscale ) {
    p->SetBorderSize( 0 );
    p->SetTextFont( 42 );
    p->SetFillStyle( 0 );
    p->SetTextColor( kBlack );
    p->SetTextSize( textscale * 0.04 );
    if( hp == LEFT ) {
      p->SetX1NDC( 0.15 );
      p->SetX2NDC( 0.4 );
      p->SetTextAlign( 12 );
    } else if( hp == CENTER ) {
      p->SetX1NDC( 0.4 );
      p->SetX2NDC( 0.5 );
      p->SetTextAlign( 32 );
    } else if( hp == CENTERRIGHT ) {
      p->SetX1NDC( 0.5 );
      p->SetX2NDC( 0.6 );
      p->SetTextAlign( 32);
    } else if( hp == RIGHTCENTER ) {
      p->SetX1NDC( 0.65 );
      p->SetX2NDC( 0.75 );
      p->SetTextAlign( 32);
    } else if( hp == RIGHT ) {
      p->SetX1NDC( 0.63 );
      p->SetX2NDC( 0.94 );
      p->SetTextAlign( 32 );
    } else assert( false );
  }

  static void setupVertical( TLegend *l , Double_t topmargin = 0.0 ) {
    Double_t textheight = Double_t(l->GetNRows()) * 1.125 * l->GetTextSize();
    l->SetY2NDC( 1.0 - 0.01 - 0.05 - topmargin );
    l->SetY1NDC( l->GetY2NDC() - textheight );
  }
  
  static void setupVertical( TPaveText *l , Double_t topmargin = 0.0 ) {
    Double_t textheight = Double_t(l->GetSize()) * 1.125 * l->GetTextSize();
    l->SetY2NDC( 1.0 - 0.01 - 0.05 - topmargin );
    l->SetY1NDC( l->GetY2NDC() - textheight );
  }

  template< typename TextObj > static void setup( TextObj *t , HorizontalPosition hp , VerticalPosition vp = TOP , Double_t textscale = 1.0 ) {
    setupHorizontal( t , hp , textscale );
    Double_t topmargin = 0.0;
    if( vp == ABOVE ) {
      topmargin = -0.015 - (0.04 * textscale);
    } else if( vp == TOP ) {
      // keep 0
    } else if( vp == BOTTOM ) {
      topmargin = 0.5;
    } else assert( false );
    setupVertical( t , topmargin );
    garbage_collector_TObject.push_back( t );
  }

  template< typename TextObj > static void setupBelow( TextObj *t , TPave *other , HorizontalPosition hp , Double_t textscale ) {
    setupHorizontal( t , hp , textscale );
    setupVertical( t , 0.05 + other->GetY2NDC() - other->GetY1NDC() );
    garbage_collector_TObject.push_back( t );
  }

  static void setup( TArrow *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetLineWidth( 2 );
  }

  static void setup( TLine *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetLineWidth( 2 );
  }

  static void setup( TAxis *h , Double_t offset_factor = 1.0 ) {
    h->SetTitleFont( GenericFont );
    h->SetLabelFont( GenericFont );
    h->SetTitleSize( GenericTitleSize );
    h->SetLabelSize( GenericLabelSize );
    h->SetTitleOffset( GenericTitleOffset * offset_factor );
  }
  
  static void setup( TProfile *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetMarkerColor( getColor(color) );
    h->SetLineWidth( 3 );
    h->SetMarkerStyle( 20+color );
  }

  static void setup( TGraph *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetMarkerColor( getColor(color) );
    h->SetLineWidth( 3 );
    h->SetMarkerStyle( 20+color );
  }

  static void setup( TGraphErrors *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetMarkerColor( getColor(color) );
    h->SetLineWidth( 3 );
    h->SetMarkerStyle( 20+color );
  }
  
  static void setup( TGraphAsymmErrors *h , std::size_t color = 0 , Int_t linewidth = 3 ) {
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetMarkerColor( getColor(color) );
    h->SetLineWidth( linewidth );
    h->SetMarkerStyle( 20+color );
  }

  static void setup( TF1 *h , std::size_t color = 0 ) {
    h->SetLineColor( getColor(color) );
    h->SetLineWidth( 3 );
  }

  static void setup( TGraph2D *h , std::size_t color = 0 ) {
    if( !h->GetHistogram() ) {
      std::cout << "trying to set up graph, but it has no histogram object yet" << std::endl;
      assert( false );
    }
    setup( h->GetXaxis() );
    setup( h->GetYaxis() );
    setup( h->GetZaxis() , 0.8 );
    h->SetLineColor( getColor(color) );
    h->SetMarkerColor( getColor(color) );
    // smooth pallete
    h->GetHistogram()->SetContour(160);
  }

  static void setupContours1D( TGraph2D *h , Double_t cl = 0.95 , std::size_t color = 0 ) {
    if( !h->GetHistogram() ) assert( false );

    h->SetLineWidth( 3 );
    h->SetLineStyle( 0 );
    h->SetLineColor( getColor(color) );

    Double_t cont[1] = { cl };
    h->GetHistogram()->SetContour( 1 , cont );
  }

  static void setupContours2D( TGraph2D *h ) {
    if( !h->GetHistogram() ) assert( false );

    h->SetLineWidth( 2 );
    h->SetLineColor( kBlack );
    h->SetLineStyle( 0 );

    Double_t cont[2] = { 2.3 , 6.0 };
    h->GetHistogram()->SetContour( 2 , cont );
  }

  static void setupContour2D( TGraph2D *h , Double_t contour , Int_t style , Int_t color ) {
    if( !h->GetHistogram() ) assert( false );
    
    h->SetLineWidth( 2 );
    h->SetLineColor( getColor(color) );
    h->SetLineStyle( style );

    Double_t cont[1] = { contour };
    h->GetHistogram()->SetContour( 1 , cont );
  }

  static void setupContour2D( TH2 *h , Double_t contour , Int_t style , Int_t color ) {
    h->SetLineWidth( 2 );
    h->SetLineColor( getColor(color) );
    h->SetLineStyle( style );
    Double_t cont[1] = { contour };
    h->SetContour( 1 , cont );
  }
  
  static void setup( TH1 *h , std::size_t color = 0 , int linestyle = 1 , int fillstyle = 0 ) {
    h->SetFillColor( getColor(color) );
    h->SetFillStyle( fillstyle );
    h->SetLineColor( fillstyle==1001 ? kBlack : getColor(color) );
    h->SetLineWidth( 3 );
    //h->SetLineWidth( 1 );
    h->SetLineStyle( linestyle );
    h->SetMarkerColor( getColor(color) );
    h->SetMarkerStyle( 0 );
    setup( h->GetXaxis() );
    setup( h->GetYaxis() );
    h->SetMinimum( 0. );
  }

  static void setup( TPie *h ) {
    h->SetTextFont( GenericFont );
    h->SetTextSize( GenericLabelSize );
    for( Int_t ientry = 0 ; ientry < h->GetEntries() ; ientry++ ) {
      h->SetEntryLineWidth( ientry , 2 );
      h->SetEntryLineColor( ientry , kBlack );
    }
    h->SetCircle( 0.5 , 0.5 , 0.35 );
  }
  
  static void setup( TEfficiency *h , std::size_t color = 0 , int linestyle = 1 , int fillstyle = 0 ) {
    h->SetFillColor( getColor(color) );
    h->SetFillStyle( fillstyle );
    h->SetLineColor( getColor(color) );
    h->SetLineWidth( 1 );
    h->SetLineStyle( linestyle );
    h->SetMarkerColor( getColor(color) );
    h->SetMarkerStyle( 0 );
  }

  
  static void setup( TH2 *h , std::size_t color = 0 ) {
    setup( h->GetXaxis() );
    setup( h->GetYaxis() );
    setup( h->GetZaxis() , 0.8 );
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetFillStyle( 0 );
    h->SetMarkerColor( getColor(color) );
    h->SetMinimum( 0. );
    h->SetContour(160);
  }

  static void setup( TH3 *h , std::size_t color = 0 ) {
    setup( h->GetXaxis() );
    setup( h->GetYaxis() );
    setup( h->GetZaxis() , 0.8 );
    h->SetLineColor( getColor(color) );
    h->SetFillColor( getColor(color) );
    h->SetFillStyle( 0 );
    h->SetMarkerColor( getColor(color) );
    h->SetMinimum( 0. );
  }

  static void setup( RooPlot *h ) {
    setup( h->GetXaxis() );
    setup( h->GetYaxis() );
  }

  static void setup( TStyle *s ) {
    s->SetPadTickY(1);
    s->SetPadTickX(1);
    s->SetOptStat(0);
    s->SetOptTitle(0);
    //s->SetPalette(510);
  }

  static void setup( TCanvas *c ) {
    c->SetTopMargin( 0.1 );
    c->SetRightMargin( 0.18 );
    c->SetBottomMargin( 0.18 );
    c->SetLeftMargin( 0.15 );    
    c->SetTicks( 0 , 1 );
  }

  static void setCanvas1D( TCanvas *c ) {
    setup(c);
    c->SetTopMargin( 0.05 );
    c->SetRightMargin( 0.05 );
    c->SetBottomMargin( 0.13 );
    c->SetLeftMargin( 0.13 );
  }

  static void setCanvas2D( TCanvas *c ) {
    setup(c);
    c->SetTopMargin( 0.05 );
    c->SetRightMargin( 0.18 );
    c->SetBottomMargin( 0.13 );
    c->SetLeftMargin( 0.13 );
  }

  static void setCanvasPie( TCanvas *c ) {
    c->SetTopMargin( 0.0 );
    c->SetRightMargin( 0.0 );
    c->SetBottomMargin( 0.0 );
    c->SetLeftMargin( 0.0 );
  }

  static void setCanvasHorizontal( TCanvas *c ) { 
    c->SetCanvasSize( 800 , 1100 );
    c->SetRightMargin( 0.05 ); 
    c->SetLeftMargin( 0.4 );
    c->SetBottomMargin( 0.05 );
  }

  static void setCanvasLog( TCanvas *c , Bool_t logx=kFALSE , Bool_t logy=kFALSE ) {
    c->SetLogx( Int_t(logx) );
    c->SetLogy( Int_t(logy) );
  }

  static void setup( TPad *p ) {
    //p->SetTopMargin( 0.1 );
    //p->SetRightMargin( 0.18 );
    //p->SetBottomMargin( 0.18 );
    //p->SetLeftMargin( 0.15 );    
    p->SetTopMargin( 0.05 );
    p->SetRightMargin( 0.05 );
    p->SetBottomMargin( 0.13 );
    p->SetLeftMargin( 0.13 );    
  }

  static void setupTop( TPad *p ) {
    setup( p );
    //p->SetRightMargin( 0.1 );
    //p->SetBottomMargin( 0.45 );
    p->SetBottomMargin( 0.4 );
    p->SetFillStyle( 4000 );
    p->SetFrameFillStyle( 4000 );
    p->SetTicks( 0 , 1 );
  }

  static void setupBottom( TPad *p ) {
    setup( p );
    //p->SetTopMargin( 0.6 );
    p->SetTopMargin( 0.65 );
    p->SetFillStyle( 4000 );
    p->SetFrameFillStyle( 4000 );
  }  
  
  static void setBottomHorizontalPad( TPad *p ) {
    setup( p );
    p->SetTopMargin( 0.18 );
    p->SetRightMargin( 0.1 );
    p->SetLeftMargin( 0.4 );
    p->SetTicks( 0 , 0 );
  }

  static void setTopHorizontalPad( TPad *p ) {
    setup( p );
    p->SetTopMargin( 0.18 );
    p->SetRightMargin( 0.1 );
    p->SetLeftMargin( 0.4 );
    p->SetFillStyle( 4000 );
    p->SetFillColor( 0 );
    p->SetFrameFillStyle( 4000 );
    p->SetTicks( 1 , 1 );
  }

  // Function to fill histogram but keep overflow/underflow bins
  template< typename ValT > static void fill( TH1* h , ValT val , Double_t w = 1.0 ) {
    assert( h );
    ValT LowerBound ( h->GetXaxis()->GetBinLowEdge( h->GetXaxis()->GetFirst() ) );
    ValT UpperBound ( h->GetXaxis()->GetBinUpEdge( h->GetXaxis()->GetLast() ) );
    if( val <= LowerBound )
      val = LowerBound + 0.0001;
    else if( val >= UpperBound )
      val = UpperBound - 0.0001;
    h->Fill( val , w );
  }


  //
  // ROOFIT MANIPULATION
  //

  void SetAllContaining( const RooAbsCollection& coll , Double_t val , TString tag = "" , Double_t otherval = 0.0 ) {
  
    TIterator *it = coll.createIterator();
    RooRealVar *v = (RooRealVar*)NULL;
  
    while( (v = (RooRealVar*)it->Next()) ) {
      TString xname( v->GetName() );
      if( xname.Contains(tag) ) {
	v->setVal( val );
      } else {
	v->setVal( otherval );
      }
    }
    
    coll.Print("v");
    v = 0;
    delete it;

  }


};



#endif

