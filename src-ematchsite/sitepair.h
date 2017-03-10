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


#ifndef __SITEPAIR_H_
#define __SITEPAIR_H_

#include<map>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<stdio.h>

#include "size.h"
#include "protein.h"
#include "rmsd.h"
#include "runsvm.h"
#include "fisherpitman.h"
#include "munkres.h"

using namespace std;

struct pair_score
{
 int res1_index;
 int res2_index;
 
 double seq_prf;
 double sec_prf;
 double hpp_prf;
 double bin_prb;
 double nbr_dis;
 double seq_ent;
 double mag_dst;
 double sco_svc;
 double sco_svr;
};

struct sup_data
{
 double sup_rms;
 double sup_t[3];
 double sup_u[3][3];
};

void SitePairMake( map< pair<int,int>, pair_score > &, Protein *, Protein * );

void SitePairSVM( map< pair<int,int>, pair_score > &, ModelSVM * );

double SitePairMunkres( map< pair<int,int>, pair_score > &, map< pair<int,int>, pair_aligned > &, int, int );

double SitePairScoreProb( map< pair<int,int>, pair_aligned > & );

double SitePairTLscore( Protein *, Protein * );

double SitePairPMscore( Protein *, Protein * );

#endif
