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


#ifndef __RUNSVM_H_
#define __RUNSVM_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cstdio>
#include<cctype>
#include<string>

#include "size.h"
#include "svm.h"

using namespace std;

class ModelSVM {
        
  private:
    
    bool _msvc_loaded; /* residue matching svc */
    bool _msvr_loaded; /* residue matching svr */
    bool _vssc_loaded; /* virtual screening scoring */
    bool _prob_loaded; /* probability of being similar */
    
    int _ems_attr;
    
    struct svm_model * _ems_model_msvc;
    struct svm_node * _ems_node_msvc;
    double _ems_scale_msvc[MAXSV1][2];
    
    struct svm_model * _ems_model_msvr;
    struct svm_node * _ems_node_msvr;
    double _ems_scale_msvr[MAXSV1][2];
    
    struct svm_model * _ems_model_vssc;
    struct svm_node * _ems_node_vssc;
    double _ems_scale_vssc[MAXSV2][2];
    
    struct svm_model * _ems_model_prob;
    struct svm_node * _ems_node_prob;
    double _ems_scale_prob[MAXSV4][2];
    
  public:
    
    ModelSVM( bool, bool, bool, bool );
    
    ModelSVM( void );
    
    ~ModelSVM();
    
    void loadModel( int, std::string );
    
    void loadScale( int, std::string );
    
    double SVMpredict( int, double [] );
};

#endif
