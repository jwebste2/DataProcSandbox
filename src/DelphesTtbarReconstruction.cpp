#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TKey.h>
#include <TSystem.h>
#include <TClonesArray.h>

#include "Report.h"
#include "ArgParser.h"
#include "HistTools.h"
#include "Plotter.h"
#include "TopDecay.h"
#include "Particle.h"
#include "DelphesBtagger.h"
#include "DelphesTruthSelector.h"
#include "DelphesRecoSelector.h"
#include "TtbarFeatureExtractor.h"
#include "TtbarLjetFeatureExtractor.h"

// Delphes Includes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;


// Unit scale factors (set specifically for Delphes input)
const Double_t GeV = 1.0;
const Double_t MeV = 1000.0;


int
main( int argc , char* argv[] )
{

  ArgParser ap( "DelphesTtbar" , "Delphes ttbar reconstruction script, creates CSV file with output features" );
  ap.addOptionalArg( "recosel" , "Name for cut on reconstructed objects" , "ljet" );
  ap.addOptionalArg( "truthsel" , "Name for cut on truth decay" , "none" );
  ap.addOptionalArg( "minjets" , "Minimum number of reco jets required" , "4" );
  ap.addOptionalArg( "maxjets" , "Maximum number of reco jets required" , "4" );
  ap.addOptionalArg( "totalSplits" , "Number of splits for dividing job" , "1000" );
  ap.addOptionalArg( "splitId" , "The split to run in this job" , "0" );
  ap.parse( argc , argv );

  // Initialize a plotting object that we can use to save histograms, etc.
  Plotter p( ap );
  p.openPs();
  p.setCanvas1D();

  // Path to the input Delphes datasets
  TString datadir ( "/cnfs/data1/users/jwebster/ttbar_01p" );
  gSystem->ExpandPathName( datadir );

  // Look at the input data directory and find all of the Delphes
  // ROOT files that contain relevent data
  vector<TString> v_inputFilePaths;
  {
    void *dirp = gSystem->OpenDirectory(datadir);
    const char *entry;
    Int_t n = 0;
    TString tmpdir;
    while( (entry = (char*)gSystem->GetDirEntry(dirp)) ) {
      tmpdir = entry;
      if( tmpdir.BeginsWith("ttbar_01p_singlecore_") ) {
	v_inputFilePaths.push_back( datadir + "/" + tmpdir + "/delphes_output.root" );
      }
    }
    delete entry;
  }
  
  // Determine which files to process in this particular job
  // This is useful in case we need to submit a bunch of jobs to
  // a cluster.
  Double_t totalSplits = ap.getAtof("totalSplits");
  Double_t splitId     = ap.getAtof("splitId");
  Int_t iFile = TMath::FloorNint( splitId * Double_t(v_inputFilePaths.size()) / totalSplits );
  Int_t fFile = TMath::FloorNint( (splitId+1.0) * Double_t(v_inputFilePaths.size()) / totalSplits );
  report::info( "Processing files with IDs in range [ %i , %i )" , iFile , fFile );

  // Some book-keeping variables
  Double_t xsec	   = 336.354802117;
  Double_t totalev = 22300922.0;
  Double_t sumw    = 0.0;
  Int_t    nPassed = 0;
  Int_t    nTotal  = 0;
  Int_t    nSolved = 0;
  Int_t    nComb   = 0;

  //
  // Initialize the histograms that we want to fill and then save
  //

  map<TString,TH1D*> h1map;
  h1map["nJets"]	 = new TH1D( "h1_nJets" , ";# of Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nLJets"]	 = new TH1D( "h1_nLJets" , ";# of Light-Tagged Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nBJets"]	 = new TH1D( "h1_nBJets" , ";# of b-Tagged Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nEl"]		 = new TH1D( "h1_nEl" , ";# of Electrons;Events / Bin" , 5 , 0 , 5 );
  h1map["nMu"]		 = new TH1D( "h1_nMu" , ";# of Muons;Events / Bin" , 5 , 0 , 5 );
  h1map["nLep"]		 = new TH1D( "h1_nLep" , ";# of Leptons;Events / Bin" , 5 , 0 , 5 );
  h1map["nJetsMatched"]	 = new TH1D( "h1_nJetsMatched" , ";# of Matched Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nLJetsMatched"] = new TH1D( "h1_nLJetsMatched" , ";# of Matched Light-Tagged Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nBJetsMatched"] = new TH1D( "h1_nBJetsMatched" , ";# of Matched b-Tagged Jets;Events / Bin" , 10 , 0 , 10 );
  h1map["nElMatched"]	 = new TH1D( "h1_nElMatched" , ";# of Matched Electrons;Events / Bin" , 5 , 0 , 5 );
  h1map["nMuMatched"]	 = new TH1D( "h1_nMuMatched" , ";# of Matched Muons;Events / Bin" , 5 , 0 , 5 );
  h1map["nLepMatched"]	 = new TH1D( "h1_nLepMatched" , ";# of Matched Leptons;Events / Bin" , 5 , 0 , 5 );
  h1map["hadTopMass"]    = new TH1D( "h1_hadTopMass" , ";M_{jjj} [GeV];Truth-Matched Hadronic Tops / Bin" , 50 , 0 , 500 );

  map<TString,TH2D*> h2map;
  h2map["topDecayMatrix"] = topdecay::getEmptyTopDecayMatrix();
  h2map["WDecayMatrix"]   = topdecay::getEmptyWDecayMatrix();

  // Initialize an object to handle Delphes B-Tagging
  DelphesBtagger *bt = new DelphesBtagger();

  // Initialize an object for constructing all features
  TtbarLjetFeatureExtractor *fe = new TtbarLjetFeatureExtractor();
  fe->openCsv( TString::Format("output/%s.csv",ap.getTag().Data()) );

  TtbarLjetFeatureExtractor *fe_sandbox = new TtbarLjetFeatureExtractor();
  fe_sandbox->openCsv( TString::Format("output/%s_sandbox.csv",ap.getTag().Data()) );

  TtbarFeatureExtractor *fe_base = new TtbarFeatureExtractor();
  fe_base->openCsv( TString::Format("output/%s_base.csv",ap.getTag().Data()) );
  

  //
  // Loop over the files
  //
  for( ; iFile < fFile ; ++iFile ) {

    TFile *f_reco = TFile::Open( v_inputFilePaths[iFile] );
    TTree *t_reco = htools::quiet_assert_load<TTree>( f_reco , "Delphes" );

    ExRootTreeReader *ex = new ExRootTreeReader( t_reco );
    
    DelphesRecoSelector *rSel  = new DelphesRecoSelector( ex , bt , ap.getAtoi("minjets") , ap.getAtoi("maxjets") );
    DelphesTruthSelector *tSel = new DelphesTruthSelector( ex );

    fe->setRecoSelector( rSel );
    fe_sandbox->setRecoSelector( rSel );
    fe_base->setRecoSelector( rSel );
    fe_base->setTruthSelector( tSel );
    
    Long64_t nev = t_reco->GetEntries();
    //assert( nev == t_truth->GetEntries() );

    report::startProgressBar( TString::Format( "File %4i" , iFile ) );
    for( Int_t iev = 0 ; iev < nev ; iev++ ) {

      report::updateProgressBar( Double_t(iev+1) , Double_t(nev) ,
				 TString::Format(" SumW=%0.1f ε=%0.3f nuMomentumSolved=%0.4f",sumw,Double_t(nPassed)/Double_t(nTotal),Double_t(nSolved)/Double_t(nPassed)) );

      nTotal++;

      ex->ReadEntry( iev );
      
      // Process the truth and reconstruction records
      // RecoSelector returns false here if there's not at least one lepton and 2 jets
      // TruthSelector always returns true here
      if( ! rSel->processRecoRecord() ) continue;
      if( ! tSel->processTruthRecord() ) continue;

      // Ensure that the event passes user-defined event selection (cleaning data).
      // Otherwise drop the event.
      if( ! rSel->passesSelection(ap["recosel"]) ) continue;
      if( ! tSel->passesSelection(ap["truthsel"]) ) continue;

      // Fill the truth decay chain information into histograms
      h2map["topDecayMatrix"]->Fill( tSel->getDecayT() , tSel->getDecayTbar() );
      h2map["WDecayMatrix"]->Fill( tSel->getDecayWp() , tSel->getDecayWm() );
      
      // Match the reco particles to truth particles using
      // delta-R matching. Each reconstructed particle is matched
      // to the nearest truth particle of the same type, as long
      // as it is within dR < 0.4.
      RecoParticle::truthMatch( rSel->getAll() , tSel->getAll() );

      // Calculate some diagnostic information that we can
      // save to histograms for validation / debugging
      int nMatched_all	 = 0;
      int nMatched_jets	 = 0;
      int nMatched_ljets = 0;
      int nMatched_bjets = 0;
      int nMatched_el	 = 0;
      int nMatched_mu	 = 0;
      for( size_t ip = 0 ; ip < rSel->getN() ; ip++ ) {
	RecoParticle *p = rSel->getParticle(ip);
	if( p->isTruthMatched() ) {
	  nMatched_all++;
	  if( p->getType()==RecoParticle::JET ) {
	    nMatched_jets++;
	    if( p->getTagLevel()>0 )
	      nMatched_bjets++;
	    else
	      nMatched_ljets++;
	  }
	  else if( p->getType()==RecoParticle::EL ) nMatched_el++;
	  else if( p->getType()==RecoParticle::MU ) nMatched_mu++;
	  else assert( false );
	}
      }

      // Calculate the reconstructed mass of the particles matched to hadronic t decay
      if( tSel->getDecayWp() == topdecay::JETS ) {
	vector<Particle*> hadtop;
	for( size_t ij = 0 ; ij < rSel->getNjets() ; ij++ ) {
	  RecoParticle *j = rSel->getJet(ij);
	  if( j->isTruthMatched() ) {
	    if( j->getTruthParticle()->isFromWp() || j->getTruthParticle()->isFromTp() ) {
	      hadtop.push_back( j );
	    }
	  }
	}
	assert( hadtop.size() <= 3 );
	if( hadtop.size()==3 ) {
	  // All of the true top jets are matched
	  h1map["hadTopMass"]->Fill( Particle::getMergedMass(hadtop) );
	} else {
	  h1map["hadTopMass"]->Fill( 0 );
	}
      }

      // Calculate the reconstructed mass of particles matched to hadronic tbar decay
      if( tSel->getDecayWm() == topdecay::JETS ) {
	vector<Particle*> hadtop;
	for( size_t ij = 0 ; ij < rSel->getNjets() ; ij++ ) {
	  RecoParticle *j = rSel->getJet(ij);
	  if( j->isTruthMatched() ) {
	    if( j->getTruthParticle()->isFromWm() || j->getTruthParticle()->isFromTm() ) {
	      hadtop.push_back( j );
	    }
	  }
	}
	assert( hadtop.size() <= 3 );
	if( hadtop.size()==3 ) {
	  // All of the true top jets are matched
	  h1map["hadTopMass"]->Fill( Particle::getMergedMass(hadtop) );
	} else {
	  h1map["hadTopMass"]->Fill( 0 );
	}
      }

      // Fill additional histograms
      h1map["nJets"]->Fill( int(rSel->getNjets()) );
      h1map["nJetsMatched"]->Fill( nMatched_jets );
      h1map["nLJets"]->Fill( int(rSel->getNljets()) );
      h1map["nLJetsMatched"]->Fill( nMatched_ljets );
      h1map["nBJets"]->Fill( int(rSel->getNbjets()) );
      h1map["nBJetsMatched"]->Fill( nMatched_bjets );
      h1map["nEl"]->Fill( int(rSel->getNel()) );
      h1map["nElMatched"]->Fill( nMatched_el );
      h1map["nMu"]->Fill( int(rSel->getNmu()) );
      h1map["nMuMatched"]->Fill( nMatched_mu );
      h1map["nLep"]->Fill( int(rSel->getNlep()) );
      h1map["nLepMatched"]->Fill( nMatched_el + nMatched_mu );

      fe_base->fill( iev , -1 , 1 );
      fe_base->save();
      
      //
      // Now the very important bit...
      //
      // Loop over all top reconstruction combinations and calculate a bunch of
      // features for each combination.
      //
      // WARNING: I assume here that we are targetting the semi-leptonic decay. The
      //          dilepton logic would look quite a bit different.
      //
      int icombo = 0;
      size_t nj = size_t( TMath::Min( int(rSel->getNjets()) , 6 ) );

      // The first thing I do is loop over all possible jets that could come from
      // the leptonic top decay to W *b*.
      //      t --> W *b* --> l v *b* 
      // Note that more often then not this jet
      // should be b-tagged, but I don't explicitly require it to be b-tagged.
      for( size_t lepTopJet = 0 ; lepTopJet < nj ; ++lepTopJet ) {
	
	// Next, loop over the remaining jets that could come from the hadronic top
	// decay to W *b*. Again, this should be b-tagged, but I don't explicitly
	// require it.
	for( size_t hadTopJet = 0 ; hadTopJet < nj ; ++hadTopJet ) {
	  if( hadTopJet == lepTopJet ) continue;

	  // Now loop over the remaining pairs of jets that could come from the
	  // hadronic W decay.
	  //    t --> W b --> *j* *j* b
	  // More often then not these should be light jets, but it is not explicitly required here.
	  // The order doesn't matter now since the 2 jets come from the same W, hence the ranges
	  // in the following 2 for loops.
	  for( size_t hadWJet1 = 0 ; hadWJet1 < nj-1 ; ++hadWJet1 ) {
	    if( hadWJet1 == lepTopJet ) continue;
	    if( hadWJet1 == hadTopJet ) continue;
	    for(  size_t hadWJet2 = hadWJet1+1 ; hadWJet2 < nj ; ++hadWJet2 ) {
	      if( hadWJet2 == lepTopJet ) continue;
	      if( hadWJet2 == hadTopJet ) continue;

	      // Check if this combination is correctly truth-matched
	      //Bool_t hadTopMatched = rSel->getJet(hadTopJet)->fromCommonWofSameTop( rSel->getJet(hadWJet1) , rSel->getJet(hadWJet2) );
	      //Bool_t lepTopMatched = rSel->getJet(lepTopJet)->fromCommonTop( rSel->getLep(0) );
	      //int signal	   = hadTopMatched && lepTopMatched ? 1 : 0;
	      //report::debug( "combo=%i : lepTopJet=%i, hadTopJet=%i, hadWJet1=%i, hadWJet2=%i : hadTopMatched=%i, leptTopMatched=%i, signal=%i" ,
	      //	     icombo , int(lepTopJet) , int(hadTopJet) , int(hadWJet1) , int(hadWJet2) , int(hadTopMatched) , int(lepTopMatched) , signal );

	      if( true ) {
	      
		fe->fill( iev , icombo , rSel->getJet(lepTopJet) , rSel->getJet(hadTopJet) , rSel->getJet(hadWJet1) , rSel->getJet(hadWJet2) );
		//fe->dump(); assert( false );
		fe->save();

		if( icombo==0 ) {
		  fe_sandbox->fill( iev , icombo , rSel->getJet(lepTopJet) , rSel->getJet(hadTopJet) , rSel->getJet(hadWJet1) , rSel->getJet(hadWJet2) );
		  fe_sandbox->save();
		}

	      }
	      
	      nComb++;
	      icombo++;
	    }
	  }
	}
      }

      nPassed++;
      if( rSel->getNuMomentumSolved() ) nSolved++;
      sumw += xsec / totalev;

    }

    delete ex;
    delete rSel;
    delete tSel;
    delete t_reco;
    f_reco->Close();
    
  }

  report::debug( "Total events = %i, Passing = %i, Solved = %i, Combinations = %i" , nTotal , nPassed , nSolved , nComb );
  
  //
  // Save the 2-dim histograms
  //
  
  p.setCanvas2D();
  gStyle->SetPaintTextFormat( "0.4f" );
  for( map<TString,TH2D*>::iterator ibeg = h2map.begin() , iend = h2map.end() ; ibeg != iend ; ++ibeg ) {
    htools::normalize( ibeg->second );
    htools::setup( ibeg->second );
    ibeg->second->SetMarkerColor( kBlack );
    ibeg->second->Draw( "colz text" );
    p.print();
    delete ibeg->second;
  }
  
  //
  // Save the 1-dim histograms
  //
  
  p.setCanvas1D();
  for( map<TString,TH1D*>::iterator ibeg = h1map.begin() , iend = h1map.end() ; ibeg != iend ; ++ibeg ) {
    htools::setup( ibeg->second );
    ibeg->second->Draw( "hist e" );
    p.print();
    delete ibeg->second;
  }

  // Close output files and exit
  fe->closeAndZipCsv();
  fe_sandbox->closeAndZipCsv();
  fe_base->closeAndZipCsv();
  p.closePs();
  
  report::info( "done." );
  return 0;

}
