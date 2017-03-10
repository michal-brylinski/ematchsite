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


#ifndef __SCREEN_H_
#define __SCREEN_H_

#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<bitset>
#include<map>
#include<sstream>
#include<gzstream.h>

#include "size.h"
#include "protein.h"
#include "runsvm.h"

using namespace std;

struct compound
{
 double mw;
 double logp;
 double psa;
 
 int hbd;
 int hba;
 
 bitset<MAXSMI> fpt1;
 bitset<MAXMAC> fpt2;
 
 double sco_tst[2];
 double sco_tsa[2];
 double sco_tsc[2];
 
 double sco_tmt[2];
 double sco_tma[2];
 double sco_tmc[2];
 
 double sco_mw[2];
 double sco_logp[2];
 double sco_psa[2];
 
 double sco_hbd[2];
 double sco_hba[2];
 
 double ran_tst[2];
 double ran_tsa[2];
 double ran_tsc[2];
 
 double ran_tmt[2];
 double ran_tma[2];
 double ran_tmc[2];
 
 double sco_svm[2];
 
 int rank[2];
};

class Screen {
        
  private:
    
    map< std::string , compound > _cmp_library; // screening library
    
    double                        _kendall_tau; // Kendall Tau
    
  public:
    
    Screen( void );
    
    ~Screen();
    
    double getKendallTau( void );
    
    void cleanCompounds( void );
    
    bool loadLibrary( std::string );
    
    void rankLibrary( ModelSVM *, Protein *, int, std::string );
    
    void calculateKendallTau( void );
};

#endif
