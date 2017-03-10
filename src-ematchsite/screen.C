/*
===============================================================================
          ______                  _        _    _           
         |  ___ \       _        | |      | |  (_)_         
     ____| | _ | | ____| |_  ____| | _     \ \  _| |_  ____ 
    / _  ) || || |/ _  |  _)/ ___) || \     \ \| |  _)/ _  )
   ( (/ /| || || ( ( | | |_( (___| | | |_____) ) | |_( (/ / 
    \____)_||_||_|\_||_|\___)____)_| |_(______/|_|\___)____)

                                                  
   eMatchSite - sequence order independent binding site alignment

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org

   Copyright 2016 Michal Brylinski

   This file is part of eMatchSite.

   eMatchSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eMatchSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eMatchSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#include "screen.h"

using namespace std;

Screen::Screen( void )
{
 _kendall_tau = 0.0;
 
 _cmp_library.clear();
}

Screen::~Screen() {}


// ==================================================================================   getKendallTau

double Screen::getKendallTau( void )
{
 return _kendall_tau;
}


// ==================================================================================   cleanCompounds

void Screen::cleanCompounds( void )
{
 _cmp_library.clear();
}


// ==================================================================================   loadLibrary

bool Screen::loadLibrary( std::string l1_name )
{
 std::string line1;
 
 igzstream library_file( l1_name.c_str() );
 
 if ( ! library_file.good() )
 {
  std::cerr << "ERROR: Opening file `" << l1_name << "' failed\n";
  return EXIT_FAILURE;
 }
 
 std::string lib1[8];
 
 int lib2 = 0;
 
 while ( getline(library_file,line1) )
  if ( lib2++ < MAXSCR )
  {
   int lib3 = 0;
   
   istringstream lib4(line1);
   
   while (lib4)
    lib4 >> lib1[lib3++];
   
   double tmw   = atof(lib1[1].c_str());
   double tlogp = atof(lib1[2].c_str());
   double tpsa  = atof(lib1[3].c_str());
   
   int thbd = atoi(lib1[4].c_str());
   int thba = atoi(lib1[5].c_str());
   
   bitset<MAXSMI> tfpt1;
   bitset<MAXMAC> tfpt2;
   
   string::iterator lib5;
   
   int lib6;
   
   for ( lib5 = lib1[6].begin(), lib6 = 0; lib5 < lib1[6].end(); lib5++, lib6++ )
    switch(*lib5)
    {
     case '0': tfpt1[MAXSMI-1-lib6] = 0; break;
     case '1': tfpt1[MAXSMI-1-lib6] = 1; break;
    }
   
   for ( lib5 = lib1[7].begin(), lib6 = 0; lib5 < lib1[7].end(); lib5++, lib6++ )
    switch(*lib5)
    {
     case '0': tfpt2[MAXMAC-1-lib6] = 0; break;
     case '1': tfpt2[MAXMAC-1-lib6] = 1; break;
    }
   
   compound tcmp;
   
   tcmp.mw = tmw;
   tcmp.logp = tlogp;
   tcmp.psa = tpsa;
   tcmp.hbd = thbd;
   tcmp.hba = thba;
   tcmp.fpt1 = tfpt1;
   tcmp.fpt2 = tfpt2;
   
   for ( int lib7 = 0; lib7 < 2; lib7++ )
   {
    tcmp.sco_tst[lib7] = 0.0;
    tcmp.sco_tsa[lib7] = 0.0;
    tcmp.sco_tsc[lib7] = 0.0;
    
    tcmp.sco_tmt[lib7] = 0.0;
    tcmp.sco_tma[lib7] = 0.0;
    tcmp.sco_tmc[lib7] = 0.0;
    
    tcmp.sco_mw[lib7] = 0.0;
    tcmp.sco_logp[lib7] = 0.0;
    tcmp.sco_psa[lib7] = 0.0;
    tcmp.sco_hbd[lib7] = 0.0;
    tcmp.sco_hba[lib7] = 0.0;
    
    tcmp.ran_tst[lib7] = 0.0;
    tcmp.ran_tsa[lib7] = 0.0;
    tcmp.ran_tsc[lib7] = 0.0;
    
    tcmp.ran_tmt[lib7] = 0.0;
    tcmp.ran_tma[lib7] = 0.0;
    tcmp.ran_tmc[lib7] = 0.0;
    
    tcmp.sco_svm[lib7] = 0.0;
    
    tcmp.rank[lib7] = 0;
   }
   
   _cmp_library.insert( pair< std::string , compound >(lib1[0],tcmp) );
  }
 
 library_file.close();
 
 if ( _cmp_library.size() < 1 )
  return EXIT_FAILURE;
 
 else
  return EXIT_SUCCESS;
}


// ==================================================================================   rankLibrary

void Screen::rankLibrary( ModelSVM * l1_model, Protein * l1_protein, int l1_number, std::string l1_method )
{
 map< std::string , compound >::iterator i1;
 
 multimap< double , std::string > score1;
 
 multimap< double , std::string > score2;
 multimap< double , std::string > score3;
 multimap< double , std::string > score4;
 multimap< double , std::string > score5;
 multimap< double , std::string > score6;
 multimap< double , std::string > score7;
 
 for ( i1 = _cmp_library.begin() ; i1 != _cmp_library.end(); i1++ )
 {
  if ( l1_method == "tst" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tst[l1_number] = l1_protein->getTanimotoSMILES( ((*i1).second).fpt1, "T" );
  
  if ( l1_method == "tsa" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tsa[l1_number] = l1_protein->getTanimotoSMILES( ((*i1).second).fpt1, "A" );
  
  if ( l1_method == "tsc" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tsc[l1_number] = l1_protein->getTanimotoSMILES( ((*i1).second).fpt1, "C" );
  
  if ( l1_method == "tmt" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tmt[l1_number] = l1_protein->getTanimotoMACCS( ((*i1).second).fpt2, "T" );
  
  if ( l1_method == "tma" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tma[l1_number] = l1_protein->getTanimotoMACCS( ((*i1).second).fpt2, "A" );
  
  if ( l1_method == "tmc" || l1_method == "svm" || l1_method == "sum" || l1_method == "max" || l1_method == "min" )
   ((*i1).second).sco_tmc[l1_number] = l1_protein->getTanimotoMACCS( ((*i1).second).fpt2, "C" );
  
  if ( l1_method == "svm" )
  {
   ((*i1).second).sco_mw[l1_number]   = 0.5 * pow( ( l1_protein->getPocketProperty(1) - ((*i1).second).mw   ) / l1_protein->getPocketProperty(2),  2.0) - log( 1.0 / ( l1_protein->getPocketProperty(2)  * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_logp[l1_number] = 0.5 * pow( ( l1_protein->getPocketProperty(3) - ((*i1).second).logp ) / l1_protein->getPocketProperty(4),  2.0) - log( 1.0 / ( l1_protein->getPocketProperty(4)  * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_psa[l1_number]  = 0.5 * pow( ( l1_protein->getPocketProperty(5) - ((*i1).second).psa  ) / l1_protein->getPocketProperty(6),  2.0) - log( 1.0 / ( l1_protein->getPocketProperty(6)  * sqrt( 2.0 * PI ) ) );
   
   ((*i1).second).sco_hbd[l1_number]  = 0.5 * pow( ( l1_protein->getPocketProperty(7) - ((*i1).second).hbd  ) / l1_protein->getPocketProperty(8),  2.0) - log( 1.0 / ( l1_protein->getPocketProperty(8)  * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_hba[l1_number]  = 0.5 * pow( ( l1_protein->getPocketProperty(9) - ((*i1).second).hba  ) / l1_protein->getPocketProperty(10), 2.0) - log( 1.0 / ( l1_protein->getPocketProperty(10) * sqrt( 2.0 * PI ) ) );
   
   double ssvm1[MAXSV2];
   
   ssvm1[0]  = ((*i1).second).sco_tst[l1_number];
   ssvm1[1]  = ((*i1).second).sco_tsa[l1_number];
   ssvm1[2]  = ((*i1).second).sco_tsc[l1_number];
   ssvm1[3]  = ((*i1).second).sco_tmt[l1_number];
   ssvm1[4]  = ((*i1).second).sco_tma[l1_number];
   ssvm1[5]  = ((*i1).second).sco_tmc[l1_number];
   ssvm1[6]  = ((*i1).second).sco_mw[l1_number];
   ssvm1[7]  = ((*i1).second).sco_logp[l1_number];
   ssvm1[8]  = ((*i1).second).sco_psa[l1_number];
   ssvm1[9]  = ((*i1).second).sco_hbd[l1_number];
   ssvm1[10] = ((*i1).second).sco_hba[l1_number];
   
   ((*i1).second).sco_svm[l1_number] = l1_model->SVMpredict( 3, ssvm1 );
  }
  
       if ( l1_method == "tst" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tst[l1_number],(*i1).first) );
  }
  else if ( l1_method == "tsa" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsa[l1_number],(*i1).first) );
  }
  else if ( l1_method == "tsc" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsc[l1_number],(*i1).first) );
  }
  else if ( l1_method == "tmt" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmt[l1_number],(*i1).first) );
  }
  else if ( l1_method == "tma" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tma[l1_number],(*i1).first) );
  }
  else if ( l1_method == "tmc" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmc[l1_number],(*i1).first) );
  }
  else if ( l1_method == "svm" )
  {
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_svm[l1_number],(*i1).first) );
  }
  else if ( l1_method == "sum" || l1_method == "max" || l1_method == "min" )
  {
   score2.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tst[l1_number],(*i1).first) );
   score3.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsa[l1_number],(*i1).first) );
   score4.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsc[l1_number],(*i1).first) );
   score5.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmt[l1_number],(*i1).first) );
   score6.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tma[l1_number],(*i1).first) );
   score7.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmc[l1_number],(*i1).first) );
  }
 }
 
 if ( l1_method == "sum" || l1_method == "max" || l1_method == "min" )
 {
  multimap< double , std::string >::iterator i3;
  
  int i4;
  
  for ( i3 = score2.begin(), i4 = 1; i3 != score2.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tst[l1_number] = i4 / (double) _cmp_library.size();
  
  for ( i3 = score3.begin(), i4 = 1; i3 != score3.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tsa[l1_number] = i4 / (double) _cmp_library.size();
  
  for ( i3 = score4.begin(), i4 = 1; i3 != score4.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tsc[l1_number] = i4 / (double) _cmp_library.size();
  
  for ( i3 = score5.begin(), i4 = 1; i3 != score5.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tmt[l1_number] = i4 / (double) _cmp_library.size();
  
  for ( i3 = score6.begin(), i4 = 1; i3 != score6.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tma[l1_number] = i4 / (double) _cmp_library.size();
  
  for ( i3 = score7.begin(), i4 = 1; i3 != score7.end(); i3++, i4++ )
   (_cmp_library.find((*i3).second)->second).ran_tmc[l1_number] = i4 / (double) _cmp_library.size();
  
  map< std::string , compound >::iterator i5;
  
  for ( i5 = _cmp_library.begin() ; i5 != _cmp_library.end(); i5++ )
  {
   double tsc1 = 0.0;
   
   if ( l1_method == "sum" )
    tsc1 = ( ((*i5).second).ran_tst[l1_number] + ((*i5).second).ran_tsa[l1_number] + ((*i5).second).ran_tsc[l1_number] + ((*i5).second).ran_tmt[l1_number] + ((*i5).second).ran_tma[l1_number] + ((*i5).second).ran_tmc[l1_number] ) / 6.0;
   else if ( l1_method == "max" )
   {
    tsc1 = ((*i5).second).ran_tst[l1_number];
    
    if ( ((*i5).second).ran_tsa[l1_number] < tsc1 ) tsc1 = ((*i5).second).ran_tsa[l1_number];
    if ( ((*i5).second).ran_tsc[l1_number] < tsc1 ) tsc1 = ((*i5).second).ran_tsc[l1_number];
    if ( ((*i5).second).ran_tmt[l1_number] < tsc1 ) tsc1 = ((*i5).second).ran_tmt[l1_number];
    if ( ((*i5).second).ran_tma[l1_number] < tsc1 ) tsc1 = ((*i5).second).ran_tma[l1_number];
    if ( ((*i5).second).ran_tmc[l1_number] < tsc1 ) tsc1 = ((*i5).second).ran_tmc[l1_number];
   }
   else if ( l1_method == "min" )
   {
    tsc1 = ((*i5).second).ran_tst[l1_number];
    
    if ( ((*i5).second).ran_tsa[l1_number] > tsc1 ) tsc1 = ((*i5).second).ran_tsa[l1_number];
    if ( ((*i5).second).ran_tsc[l1_number] > tsc1 ) tsc1 = ((*i5).second).ran_tsc[l1_number];
    if ( ((*i5).second).ran_tmt[l1_number] > tsc1 ) tsc1 = ((*i5).second).ran_tmt[l1_number];
    if ( ((*i5).second).ran_tma[l1_number] > tsc1 ) tsc1 = ((*i5).second).ran_tma[l1_number];
    if ( ((*i5).second).ran_tmc[l1_number] > tsc1 ) tsc1 = ((*i5).second).ran_tmc[l1_number];
   }
   
   score1.insert( pair< double , std::string >(tsc1,(*i5).first) );
  }
 }
 
 multimap< double , std::string >::iterator i6;
 
 int i7;
 
 for ( i6 = score1.begin(), i7 = 1; i6 != score1.end(); i6++, i7++ )
  (_cmp_library.find((*i6).second)->second).rank[l1_number] = i7;
}


// ==================================================================================   calculateKendallTau

void Screen::calculateKendallTau( void )
{
 double ranks1[MAXSCR];
 double ranks2[MAXSCR];
 
 int ranks3 = 0;
 
 map< std::string , compound >::iterator i1;
 
 for ( i1 = _cmp_library.begin() ; i1 != _cmp_library.end(); i1++ )
 {
  ranks1[ranks3] = ((*i1).second).rank[0];
  ranks2[ranks3] = ((*i1).second).rank[1];
  
  ranks3++;
 }
 
 _kendall_tau = 0.0;
 
 for ( int i = 0; i < ranks3; i++ )
 {
  for ( int j = i + 1; j < ranks3; j++ )
  {
   double prod = ( ranks1[i] - ranks1[j] ) * ( ranks2[i] - ranks2[j] );
   
   _kendall_tau += ( prod > 0 ) - ( prod < 0 );
  }
 }
 
 _kendall_tau /= ( 0.5 * ( (double) ranks3 * ( (double) ranks3 - 1.0 ) ) );
}
