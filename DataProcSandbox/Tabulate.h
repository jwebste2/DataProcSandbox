#ifndef _TABULATE_H_
#define _TABULATE_H_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <TString.h>

#include "Report.h"
#include "ArgParser.h"

class 
Tabulate
{
  
private:

  struct Row {
    std::vector<TString> el;
    Double_t sort_stat;

    Row() : sort_stat(1.0e7) {}

    bool operator<( const Row &other ) const { return sort_stat < other.sort_stat; }
      
  };

  
  TString col_format;
  std::vector<Row> vec_rows;
  TString psname;


  std::vector<TString> produceLatexVector() {

    std::vector<TString> tab;
    tab.push_back( "\\documentclass[a4paper]{article}" );
    tab.push_back( "\\usepackage[usenames]{color}" );
    tab.push_back( "\\usepackage{amssymb}" );
    tab.push_back( "\\usepackage{amsmath}" );
    tab.push_back( "\\usepackage[utf8]{inputenc}" );
    //tab.push_back( "\\usepackage{feynmp}" );
    tab.push_back( "\\usepackage{multirow}" );
    tab.push_back( "\\usepackage{graphicx}" );
    tab.push_back( "\\usepackage{geometry}" );
    tab.push_back( "\\special{papersize=20in,15in}" );
    tab.push_back( "\\begin{document}" );
    tab.push_back( "\\section{FIXME}" );
    tab.push_back( "" );
    tab.push_back( "\\begin{table}[htbp]" );
    tab.push_back( "\\scriptsize" );
    tab.push_back( "\\centering" );
    tab.push_back( "\\begin{tabular}{"+col_format+"}" );

    // dump rows
    for( std::size_t irow = 0 ; irow < vec_rows.size() ; irow++ ) {

      Row r = vec_rows[irow];

      if( r.el.size()==0 ) {

	// if the row is empty, just do a horizontal line
	tab.push_back( " \\hline" );

      } else {

	// loop over elements in r.el
	TString l = "";
	for( std::size_t icol = 0 ; icol < r.el.size() ; icol++ ) {
	  if( icol > 0 ) l += " &";
	  l += " " + r.el[icol];
	}
	l += " \\\\";
	tab.push_back( l );

      }

    }

    // dump tail
    tab.push_back( "\\end{tabular}" );
    tab.push_back( "\\caption{FIXME}" );
    tab.push_back( "\\label{FIXME}" );
    tab.push_back( "\\end{table}" );
    tab.push_back( "" );
    tab.push_back( "\\end{document}" );

    return tab;
    
  }

  
  TString getCompileCommands() {
    TString c = "";
    c += "  cd output/; ";
    c += "rm *.dvi *.aux *.log; ";
    c += "latex " + psname + "; ";
    c += "dvipdf " + psname + "; ";
    c += "cd ..;";
    return c;
  }
    
  TString getDownloadCommands() {
    TString c = "";
    c += "  anl-get /users/jwebster/tth/tthAna/output/" + psname + ".pdf; ";
    c += "preview " + psname + ".pdf";
    return c;
  }
  
public:
  
  Tabulate( ArgParser ap , TString _col_format )
    : col_format( _col_format )
    , psname( ap.getTag() )
  {}
    
  ~Tabulate() {}

    void addRow( TString s1 = "" , TString s2 = "" , TString s3 = "" , TString s4 = "" , TString s5 = "" , TString s6 = "" , TString s7 = "" , TString s8 = "" ) {
    Row r;
    if( s1.Length() ) r.el.push_back( s1 );
    if( s2.Length() ) r.el.push_back( s2 );
    if( s3.Length() ) r.el.push_back( s3 );
    if( s4.Length() ) r.el.push_back( s4 );
    if( s5.Length() ) r.el.push_back( s5 );
    if( s6.Length() ) r.el.push_back( s6 );
    if( s7.Length() ) r.el.push_back( s7 );
    if( s8.Length() ) r.el.push_back( s8 );
    vec_rows.push_back( r );
  }

  void addRow( Double_t d_sort , TString s1 = "" , TString s2 = "" , TString s3 = "" , TString s4 = "" , TString s5 = "" , TString s6 = "" , TString s7 = "" , TString s8 = "" ) {
    Row r;
    r.sort_stat = d_sort;
    if( s1.Length() ) r.el.push_back( s1 );
    if( s2.Length() ) r.el.push_back( s2 );
    if( s3.Length() ) r.el.push_back( s3 );
    if( s4.Length() ) r.el.push_back( s4 );
    if( s5.Length() ) r.el.push_back( s5 );
    if( s6.Length() ) r.el.push_back( s6 );
    if( s7.Length() ) r.el.push_back( s7 );
    if( s8.Length() ) r.el.push_back( s8 );
    vec_rows.push_back( r );
  }

  void startRow() { Row newrow; vec_rows.push_back( newrow ); }
  void addPartOfRow( TString s ) { vec_rows.back().el.push_back( s ); }

  
  void sortRows() {

    // Figure out the first sortable element, and sort from there to the end of the list of Rows

    std::vector<Row>::iterator it = vec_rows.begin();
    while( it->sort_stat > 1.0e6 ) it++;

    if( it < vec_rows.end() ) {
      std::sort( it , vec_rows.end() );
    }

  }


  void dumpBash() {

    report::info( "dumping table to terminal" );

    std::cout << std::endl;
    
    for( std::size_t irow = 0 ; irow < vec_rows.size() ; irow++ ) {

      Row r = vec_rows[irow];
      
      if( r.el.size()==0 ) {

	// if the row is empty, just do a horizontal line
	for( Int_t i = 0 ; i < 80 ; i++ ) std::cout << "-";
	std::cout << std::endl;

      } else {
	
	// loop over elements in r.el
	for( std::size_t icol = 0 ; icol < r.el.size() ; icol++ ) {
	  std::cout << " " << r.el[icol] << std::flush;
	}
	
      }

    }
    
    std::cout << std::endl;
    
  }
  
  
  void dumpLatex() {
    
    report::info( "dumping table in latex format" );

    std::vector<TString> tab = produceLatexVector();
    std::cout << std::endl;
    for( TString l : tab ) {
      std::cout << l << std::endl;
    }
    std::cout << std::endl;

  }


  void save() {

    report::info( "writing table to file" );
    ofstream file( TString::Format("output/%s.tex",psname.Data()).Data() );

    std::vector<TString> tab = produceLatexVector();
    for( TString l : tab ) {
      file << l << std::endl;
    }

    file.close();

    std::cout << getCompileCommands() << std::endl;
    std::cout << getDownloadCommands() << std::endl;
    
  }

  void compile() {
    std::system( getCompileCommands().Data() );
    std::cout << getDownloadCommands() << std::endl;
  }
  
};

#endif
