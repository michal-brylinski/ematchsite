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


#include "fisherpitman.h"

using namespace std;

double FisherPitman( int fp_rounds, vector<double> &fp_dst1, vector<double> &fp_dst2 )
{
 srand(1);
 
 double fp_T = 0.0;
 
 double fp_sum = 0.0;
 
 vector<double> fp_perm;
 
 fp_perm.clear();
 
 vector<double>::iterator it1;
 
 for ( it1 = fp_dst1.begin(); it1 < fp_dst1.end(); it1++ )
 {
  fp_sum += (*it1);
  
  fp_perm.push_back(*it1);
 }
 
 fp_T += ( (double) fp_dst1.size() * pow( 1.0 / (double) fp_dst1.size() * fp_sum, 2.0) );
 
 fp_sum = 0.0;
 
 for ( it1 = fp_dst2.begin(); it1 < fp_dst2.end(); it1++ )
 {
  fp_sum += (*it1);
  
  fp_perm.push_back(*it1);
 }
 
 fp_T += ( (double) fp_dst2.size() * pow( 1.0 / (double) fp_dst2.size() * fp_sum, 2.0) );
 
 int fp_wge = 0;
 
 for ( int it2 = 0; it2 < fp_rounds; it2++ )
 {
  std::random_shuffle(fp_perm.begin(), fp_perm.end());
  
  double fp_sum1 = 0.0;
  double fp_sum2 = 0.0;
  
  unsigned int it3 = 0;
  
  for ( it1 = fp_perm.begin(); it1 < fp_perm.end(); it1++ )
  {
   if ( it3++ < fp_dst1.size() )
    fp_sum1 += (*it1);
   else
    fp_sum2 += (*it1);
  }
  
  double fp_P = ( (double) fp_dst1.size() * pow( 1.0 / (double) fp_dst1.size() * fp_sum1, 2.0) ) + 
                ( (double) fp_dst2.size() * pow( 1.0 / (double) fp_dst2.size() * fp_sum2, 2.0) );
  
  if ( fp_P >= fp_T )
   fp_wge++;
 }
 
 return (double) fp_wge / (double) fp_rounds;
}
