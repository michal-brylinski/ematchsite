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


#ifndef __COORDS_H_
#define __COORDS_H_

#include<string>
#include<deque>

using namespace std;

// ==================================================================================   CoordsProtein

class CoordsProtein {
        
  private:
    
    int    n; // atom number
    int    r; // residue number
    double x; // x coord
    double y; // y coord
    double z; // z coord
    string t; // residue name
    string a; // atom name
    string c; // chain id

  public:

    CoordsProtein( int, int, double, double, double, string, string, string );
    
    CoordsProtein( void );

    ~CoordsProtein();

    int getAtomNumber();
    
    int getResidueNumber();
    
    double getCoords( int );
    
    void setCoords( double, double, double );
    
    string getResidueName();
    
    string getAtomName();
    
    string getChainID();
};

// ==================================================================================   CoordsLigand

class CoordsLigand {
        
  private:
    
    int    n; // atom number
    double x; // x coordinate
    double y; // y coordinate
    double z; // z coordinate
    string a; // atom name

  public:

    CoordsLigand( int, double, double, double, string );
    
    CoordsLigand( void );

    ~CoordsLigand();

    int getAtomNumber();
    
    double getCoords( int );
    
    void setCoords( double, double, double );
    
    string getAtomName();
};

#endif
