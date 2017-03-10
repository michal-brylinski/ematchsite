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


#include "runsvm.h"

using namespace std;

ModelSVM::ModelSVM( bool aa, bool ab, bool ac, bool ad )
{
 _msvc_loaded = aa;
 _msvr_loaded = ab;
 _vssc_loaded = ac;
 _prob_loaded = ad;
 
 _ems_attr = 64;
 
 _ems_node_msvc = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_msvr = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_vssc = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_prob = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
}

ModelSVM::ModelSVM( void )
{
 _msvc_loaded = false;
 _msvr_loaded = false;
 _vssc_loaded = false;
 _prob_loaded = false;
 
 _ems_attr = 64;
 
 _ems_node_msvc = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_msvr = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_vssc = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
 _ems_node_prob = ( struct svm_node * ) malloc( _ems_attr * sizeof( struct svm_node ) );
}

ModelSVM::~ModelSVM() {}


// ==================================================================================   loadModel

void ModelSVM::loadModel( int m1_num, std::string m1_name )
{
 switch ( m1_num )
 {
  case 1: _ems_model_msvc = svm_load_model( m1_name.c_str() ); break;
  case 2: _ems_model_msvr = svm_load_model( m1_name.c_str() ); break;
  case 3: _ems_model_vssc = svm_load_model( m1_name.c_str() ); break;
  case 4: _ems_model_prob = svm_load_model( m1_name.c_str() ); break;
 }
}


// ==================================================================================   loadScale

void ModelSVM::loadScale( int s1_num, std::string s1_name )
{
 string line1, line2;
 
 ifstream s1_file( s1_name.c_str() );
 
 int i1 = 0;
 
 while (getline(s1_file,line1))
  if ( i1++ > 1 )
  {
   istringstream iss(line1, istringstream::in);
   
   int i2 = 0;
   int i3 = 0;
   double i4 = 0.0;
   double i5 = 0.0;
   
   while( iss >> line2 )     
   {
         if ( i2 == 0 ) { i3 = atoi(line2.c_str()) - 1; }
    else if ( i2 == 1 ) { i4 = atof(line2.c_str()); }
    else if ( i2 == 2 ) { i5 = atof(line2.c_str()); }
    
    i2++;
   }
   
   switch ( s1_num )
   {
    case 1: _ems_scale_msvc[i3][0] = i4; _ems_scale_msvc[i3][1] = i5; break;
    case 2: _ems_scale_msvr[i3][0] = i4; _ems_scale_msvr[i3][1] = i5; break;
    case 3: _ems_scale_vssc[i3][0] = i4; _ems_scale_vssc[i3][1] = i5; break;
    case 4: _ems_scale_prob[i3][0] = i4; _ems_scale_prob[i3][1] = i5; break;
   }
  }
 
 s1_file.close();
}


// ==================================================================================   SVMpredict

double ModelSVM::SVMpredict( int m1_num, double m1_fet [] )
{
 int m2; m2 = 0;
 
 switch ( m1_num )
 {
  case 1: m2 = MAXSV1; break;
  case 2: m2 = MAXSV1; break;
  case 3: m2 = MAXSV2; break;
  case 4: m2 = MAXSV4; break;
 }
 
 for ( int m3 = 0; m3 < m2; m3++ )
 {
  double lb; lb = 0.0;
  double ub; ub = 0.0;
  
  switch ( m1_num )
  {
   case 1: lb = _ems_scale_msvc[m3][0]; ub = _ems_scale_msvc[m3][1]; break;
   case 2: lb = _ems_scale_msvr[m3][0]; ub = _ems_scale_msvr[m3][1]; break;
   case 3: lb = _ems_scale_vssc[m3][0]; ub = _ems_scale_vssc[m3][1]; break;
   case 4: lb = _ems_scale_prob[m3][0]; ub = _ems_scale_prob[m3][1]; break;
  }
  
  m1_fet[m3] = ( ( m1_fet[m3] - lb ) / ( ub - lb ) ) * 2.0 - 1.0;
  
  if ( m1_fet[m3] < -1.0 )
   m1_fet[m3] = -1.0;
  if ( m1_fet[m3] > 1.0 )
   m1_fet[m3] = 1.0;
 }
 
 int nr_class = 0;
 
 double *prob_estimates = NULL;
 
 int *labels = ( int * ) malloc( nr_class * sizeof( int ) );
 
 switch ( m1_num )
 {
  case 1: nr_class = svm_get_nr_class(_ems_model_msvc); svm_get_labels(_ems_model_msvc,labels); break;
  case 2: nr_class = svm_get_nr_class(_ems_model_msvr); svm_get_labels(_ems_model_msvr,labels); break;
  case 3: nr_class = svm_get_nr_class(_ems_model_vssc); svm_get_labels(_ems_model_vssc,labels); break;
  case 4: nr_class = svm_get_nr_class(_ems_model_prob); svm_get_labels(_ems_model_prob,labels); break;
 }
 
 prob_estimates = ( double * ) malloc( nr_class * sizeof( double ) );
 
 int r1 = 0;
 
 if ( labels[1] )
  r1 = 1;
 
 free(labels);
 
 for ( int i = 0; i < m2; i++ )
 {
  switch ( m1_num )
  {
   case 1: _ems_node_msvc[i].index = i + 1; _ems_node_msvc[i].value = m1_fet[i]; break;
   case 2: _ems_node_msvr[i].index = i + 1; _ems_node_msvr[i].value = m1_fet[i]; break;
   case 3: _ems_node_vssc[i].index = i + 1; _ems_node_vssc[i].value = m1_fet[i]; break;
   case 4: _ems_node_prob[i].index = i + 1; _ems_node_prob[i].value = m1_fet[i]; break;
  }
 }
 
 double p1 = 0.0;
 double p2 = 0.0;
 
 switch ( m1_num )
 {
  case 1: _ems_node_msvc[m2].index = -1; p1 = svm_predict_probability(_ems_model_msvc,_ems_node_msvc,prob_estimates); p2 = prob_estimates[r1]; break;
  case 2: _ems_node_msvr[m2].index = -1; p1 = svm_predict(_ems_model_msvr,_ems_node_msvr); break;
  case 3: _ems_node_vssc[m2].index = -1; p1 = svm_predict_probability(_ems_model_vssc,_ems_node_vssc,prob_estimates); p2 = prob_estimates[r1]; break;
  case 4: _ems_node_prob[m2].index = -1; p1 = svm_predict_probability(_ems_model_prob,_ems_node_prob,prob_estimates); p2 = prob_estimates[r1]; break;
 }
 
 free(prob_estimates);
 
 if ( m1_num == 2 )
 {
  p2 = p1;
 }
 else
 {
  if ( p2 < 0.0 )
   p2 = 0.0;
  if ( p2 > 1.0 )
  p2 = 1.0;
 }
 
 return p2;
}
