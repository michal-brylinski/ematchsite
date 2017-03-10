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


#include "data.h"

using namespace std;


std::string three2oneS( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return "A"; }
 else if ( resnam1 == "CYS" ) { return "C"; }
 else if ( resnam1 == "ASP" ) { return "D"; }
 else if ( resnam1 == "GLU" ) { return "E"; }
 else if ( resnam1 == "PHE" ) { return "F"; }
 else if ( resnam1 == "GLY" ) { return "G"; }
 else if ( resnam1 == "HIS" ) { return "H"; }
 else if ( resnam1 == "ILE" ) { return "I"; }
 else if ( resnam1 == "LYS" ) { return "K"; }
 else if ( resnam1 == "LEU" ) { return "L"; }
 else if ( resnam1 == "MET" ) { return "M"; }
 else if ( resnam1 == "ASN" ) { return "N"; }
 else if ( resnam1 == "PRO" ) { return "P"; }
 else if ( resnam1 == "GLN" ) { return "Q"; }
 else if ( resnam1 == "ARG" ) { return "R"; }
 else if ( resnam1 == "SER" ) { return "S"; }
 else if ( resnam1 == "THR" ) { return "T"; }
 else if ( resnam1 == "VAL" ) { return "V"; }
 else if ( resnam1 == "TRP" ) { return "W"; }
 else if ( resnam1 == "TYR" ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

char three2oneC( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return 'A'; }
 else if ( resnam1 == "CYS" ) { return 'C'; }
 else if ( resnam1 == "ASP" ) { return 'D'; }
 else if ( resnam1 == "GLU" ) { return 'E'; }
 else if ( resnam1 == "PHE" ) { return 'F'; }
 else if ( resnam1 == "GLY" ) { return 'G'; }
 else if ( resnam1 == "HIS" ) { return 'H'; }
 else if ( resnam1 == "ILE" ) { return 'I'; }
 else if ( resnam1 == "LYS" ) { return 'K'; }
 else if ( resnam1 == "LEU" ) { return 'L'; }
 else if ( resnam1 == "MET" ) { return 'M'; }
 else if ( resnam1 == "ASN" ) { return 'N'; }
 else if ( resnam1 == "PRO" ) { return 'P'; }
 else if ( resnam1 == "GLN" ) { return 'Q'; }
 else if ( resnam1 == "ARG" ) { return 'R'; }
 else if ( resnam1 == "SER" ) { return 'S'; }
 else if ( resnam1 == "THR" ) { return 'T'; }
 else if ( resnam1 == "VAL" ) { return 'V'; }
 else if ( resnam1 == "TRP" ) { return 'W'; }
 else if ( resnam1 == "TYR" ) { return 'Y'; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string one2three( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return "ALA"; }
 else if ( resnam1 == "C" ) { return "CYS"; }
 else if ( resnam1 == "D" ) { return "ASP"; }
 else if ( resnam1 == "E" ) { return "GLU"; }
 else if ( resnam1 == "F" ) { return "PHE"; }
 else if ( resnam1 == "G" ) { return "GLY"; }
 else if ( resnam1 == "H" ) { return "HIS"; }
 else if ( resnam1 == "I" ) { return "ILE"; }
 else if ( resnam1 == "K" ) { return "LYS"; }
 else if ( resnam1 == "L" ) { return "LEU"; }
 else if ( resnam1 == "M" ) { return "MET"; }
 else if ( resnam1 == "N" ) { return "ASN"; }
 else if ( resnam1 == "P" ) { return "PRO"; }
 else if ( resnam1 == "Q" ) { return "GLN"; }
 else if ( resnam1 == "R" ) { return "ARG"; }
 else if ( resnam1 == "S" ) { return "SER"; }
 else if ( resnam1 == "T" ) { return "THR"; }
 else if ( resnam1 == "V" ) { return "VAL"; }
 else if ( resnam1 == "W" ) { return "TRP"; }
 else if ( resnam1 == "Y" ) { return "TYR"; }
 
 else
 {
  cout << "Unknown residue passed to one2three: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int one2num( char resnam1 )
{
      if ( resnam1 == 'A' ) { return  0; }
 else if ( resnam1 == 'C' ) { return  1; }
 else if ( resnam1 == 'D' ) { return  2; }
 else if ( resnam1 == 'E' ) { return  3; }
 else if ( resnam1 == 'F' ) { return  4; }
 else if ( resnam1 == 'G' ) { return  5; }
 else if ( resnam1 == 'H' ) { return  6; }
 else if ( resnam1 == 'I' ) { return  7; }
 else if ( resnam1 == 'K' ) { return  8; }
 else if ( resnam1 == 'L' ) { return  9; }
 else if ( resnam1 == 'M' ) { return 10; }
 else if ( resnam1 == 'N' ) { return 11; }
 else if ( resnam1 == 'P' ) { return 12; }
 else if ( resnam1 == 'Q' ) { return 13; }
 else if ( resnam1 == 'R' ) { return 14; }
 else if ( resnam1 == 'S' ) { return 15; }
 else if ( resnam1 == 'T' ) { return 16; }
 else if ( resnam1 == 'V' ) { return 17; }
 else if ( resnam1 == 'W' ) { return 18; }
 else if ( resnam1 == 'Y' ) { return 19; }
 
 else
 {
  cout << "Unknown residue passed to one2num: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int one2group( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return 0; }
 else if ( resnam1 == "C" ) { return 4; }
 else if ( resnam1 == "D" ) { return 2; }
 else if ( resnam1 == "E" ) { return 2; }
 else if ( resnam1 == "F" ) { return 3; }
 else if ( resnam1 == "G" ) { return 0; }
 else if ( resnam1 == "H" ) { return 1; }
 else if ( resnam1 == "I" ) { return 0; }
 else if ( resnam1 == "K" ) { return 1; }
 else if ( resnam1 == "L" ) { return 0; }
 else if ( resnam1 == "M" ) { return 0; }
 else if ( resnam1 == "N" ) { return 2; }
 else if ( resnam1 == "P" ) { return 0; }
 else if ( resnam1 == "Q" ) { return 2; }
 else if ( resnam1 == "R" ) { return 1; }
 else if ( resnam1 == "S" ) { return 4; }
 else if ( resnam1 == "T" ) { return 4; }
 else if ( resnam1 == "V" ) { return 0; }
 else if ( resnam1 == "W" ) { return 3; }
 else if ( resnam1 == "Y" ) { return 3; }
 
 else
 {
  cout << "Unknown residue passed to one2group: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string num2one( int resnum1 )
{
      if ( resnum1 ==  0 ) { return "A"; }
 else if ( resnum1 ==  1 ) { return "C"; }
 else if ( resnum1 ==  2 ) { return "D"; }
 else if ( resnum1 ==  3 ) { return "E"; }
 else if ( resnum1 ==  4 ) { return "F"; }
 else if ( resnum1 ==  5 ) { return "G"; }
 else if ( resnum1 ==  6 ) { return "H"; }
 else if ( resnum1 ==  7 ) { return "I"; }
 else if ( resnum1 ==  8 ) { return "K"; }
 else if ( resnum1 ==  9 ) { return "L"; }
 else if ( resnum1 == 10 ) { return "M"; }
 else if ( resnum1 == 11 ) { return "N"; }
 else if ( resnum1 == 12 ) { return "P"; }
 else if ( resnum1 == 13 ) { return "Q"; }
 else if ( resnum1 == 14 ) { return "R"; }
 else if ( resnum1 == 15 ) { return "S"; }
 else if ( resnum1 == 16 ) { return "T"; }
 else if ( resnum1 == 17 ) { return "V"; }
 else if ( resnum1 == 18 ) { return "W"; }
 else if ( resnum1 == 19 ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to num2one: " << resnum1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

void three2hpp( std::string resnam1, double reshpp1[] )
{
      if ( resnam1 == "ALA" ) { double reshpp2[20] = { 0.440, 0.610, 0.100,5.330, 0.390,  1.940,-0.300, 0.620,-0.500,12.970,0.616, 0.310, 0.300,1.360, 0.620, 1.150,  2.100, 0.420, 0.350, 1.800 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "CYS" ) { double reshpp2[20] = { 0.580, 0.360,-1.420,7.930, 0.250, -1.240, 6.300, 0.290,-1.000,14.630,0.680, 1.540, 0.900,1.270, 0.290,-1.200,  1.400, 0.840, 0.760, 2.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "ASP" ) { double reshpp2[20] = {-0.310, 0.610, 0.780,3.590,-3.810,-10.950,-1.400,-0.900, 3.000,10.850,0.028,-0.770,-0.600,0.110,-0.090, 0.650, 10.000,-0.510,-2.150,-3.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "GLU" ) { double reshpp2[20] = {-0.340, 0.510, 0.830,3.650,-2.910,-10.200, 0.000,-0.740, 3.000,11.890,0.043,-0.640,-0.700,0.250,-0.740,-0.710,  7.800,-0.370,-1.950,-3.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "PHE" ) { double reshpp2[20] = { 2.540,-1.520,-2.120,9.030, 2.270, -0.760, 7.500, 1.190,-2.500,14.000,1.000, 1.790, 0.500,1.570, 1.190,-1.410, -9.200, 1.740, 1.690, 2.800 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "GLY" ) { double reshpp2[20] = { 0.000, 0.810, 0.330,4.480, 0.000,  2.390, 1.200, 0.480, 0.000,12.430,0.501, 0.000, 0.300,1.090, 0.480,-1.840,  5.700, 0.000, 0.000,-0.400 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "HIS" ) { double reshpp2[20] = {-0.010, 0.690,-0.500,5.100,-0.640,-10.270,-1.300,-0.400,-0.500,12.160,0.165, 0.130,-0.100,0.680,-0.400, 3.120,  2.100,-2.280,-0.650,-3.200 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "ILE" ) { double reshpp2[20] = { 2.460,-1.450,-1.130,8.830, 1.820,  2.150, 4.300, 1.380,-1.800,15.670,0.943, 1.800, 0.700,1.440, 1.380,-2.920, -8.000, 1.810, 1.830, 4.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "LYS" ) { double reshpp2[20] = {-2.450, 0.460, 1.400,2.950,-2.770, -9.520,-3.600,-1.500, 3.000,11.360,0.283,-0.990,-1.800,0.090,-1.500, 2.060,  5.700,-2.030,-1.540,-3.900 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "LEU" ) { double reshpp2[20] = { 2.460,-1.650,-1.180,8.470, 1.820,  2.280, 6.600, 1.060,-1.800,14.900,0.943, 1.700, 0.500,1.470, 1.530, 0.750, -9.200, 1.800, 1.800, 3.800 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "MET" ) { double reshpp2[20] = { 1.100,-0.660,-1.590,8.950, 0.960, -1.480, 2.500, 0.640,-1.300,14.390,0.738, 1.230, 0.400,1.420, 0.640,-3.850, -4.200, 1.180, 1.100, 1.900 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "ASN" ) { double reshpp2[20] = {-1.320, 0.890, 0.480,3.710,-1.910, -9.680,-0.200,-0.780, 0.200,11.420,0.236,-0.600,-0.500,0.330,-0.780,-0.770,  7.000,-1.030,-0.990,-3.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "PRO" ) { double reshpp2[20] = { 1.290,-0.170, 0.730,3.870, 0.990,  0.000, 2.200, 0.120, 0.000,11.370,0.711, 0.720,-0.300,0.540, 0.120,-0.530,  2.100, 0.860, 0.840,-1.600 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "GLN" ) { double reshpp2[20] = {-0.710, 0.970, 0.950,3.870,-1.300, -9.380,-0.200,-0.850, 0.200,11.760,0.251,-0.220,-0.700,0.330,-0.850,-0.110,  6.000,-0.960,-0.930,-3.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "ARG" ) { double reshpp2[20] = {-2.420, 0.690, 1.910,4.180,-3.950,-19.920,-1.100,-2.530, 3.000,11.720,0.000,-1.010,-1.400,0.150,-2.530, 0.580,  4.200,-1.560,-1.500,-4.500 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "SER" ) { double reshpp2[20] = {-0.840, 0.420, 0.520,4.090,-1.240, -5.060,-0.600,-0.180, 0.300,11.230,0.359,-0.040,-0.100,0.970,-0.180,-0.260,  6.500,-0.640,-0.630,-0.800 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "THR" ) { double reshpp2[20] = {-0.410, 0.290, 0.070,4.490,-1.000, -4.880,-2.200,-0.050,-0.400,11.690,0.450, 0.260,-0.200,1.080,-0.050,-0.450,  5.200,-0.260,-0.270,-0.700 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "VAL" ) { double reshpp2[20] = { 1.730,-0.750,-1.270,7.630, 1.300,  1.990, 5.900, 1.080,-1.500,15.710,0.825, 1.220, 0.600,1.370, 1.800,-0.130, -3.700, 1.340, 1.320, 4.200 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "TRP" ) { double reshpp2[20] = { 2.560,-1.200,-0.510,7.660, 2.130, -5.880, 7.900, 0.810,-3.400,13.930,0.878, 2.250, 0.300,1.000, 0.810,-1.140,-10.000, 1.460, 1.350,-0.900 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 else if ( resnam1 == "TYR" ) { double reshpp2[20] = { 1.630,-1.430,-0.210,5.890, 1.470, -6.110, 7.100, 0.260,-2.300,13.420,0.880, 0.960,-0.400,0.830, 0.260, 0.130, -1.900, 0.510, 0.390,-1.300 }; for ( int i = 0; i < 20; i++ ) reshpp1[i] = reshpp2[i]; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
}
