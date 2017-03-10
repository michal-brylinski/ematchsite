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


#include "sitepair.h"

using namespace std;


// ==================================================================================   SitePairMake

void SitePairMake( map< pair<int,int>, pair_score > &aln_pair, Protein * target1, Protein * target2 )
{
 aln_pair.clear();
 
 map< pair<int,int>, sup_data > sup_matrix;
 
 for ( int ip1 = 0; ip1 < target1->getBindingResiduesTotal(); ip1++ )
 {
  double res1_xyz[3];
  
  target1->getBindingResCoords(ip1, res1_xyz, 0);
  
  for ( int ip2 = 0; ip2 < target2->getBindingResiduesTotal(); ip2++ )
  {
   double res2_xyz[3];
   
   target2->getBindingResCoords(ip2, res2_xyz, 0);
   
   int res1 = target1->getBindingResNum(ip1);
   int res2 = target2->getBindingResNum(ip2);
   
   double seq_prof1[20];
   double seq_prof2[20];
   
   target1->getSeqProfile(res1-1, seq_prof1);
   target2->getSeqProfile(res2-1, seq_prof2);
   
   pair_score scores = { ip1, ip2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   
   /* Sequence profile score */
   
   for ( int i1 = 0; i1 < 20; i1++ )
    scores.seq_prf += ( seq_prof1[i1] * seq_prof2[i1] );
   
   /* Secondary structure score */
   
   double sec_prof1[3];
   double sec_prof2[3];
   
   target1->getSecProfile(res1-1, sec_prof1);
   target2->getSecProfile(res2-1, sec_prof2);
   
   scores.sec_prf = sqrt( pow( sec_prof1[0] - sec_prof2[0], 2.0 ) + pow( sec_prof1[1] - sec_prof2[1], 2.0 ) + pow( sec_prof1[2] - sec_prof2[2], 2.0 ) );
   
   /* Hydrophobicity score */
   
   double hpp_prof1[20];
   double hpp_prof2[20];
   
   target1->getHPProfile(res1-1, hpp_prof1);
   target2->getHPProfile(res2-1, hpp_prof2);
   
   double cc_sx  = 0.0;
   double cc_sy  = 0.0;
   double cc_sxy = 0.0;
   double cc_sxx = 0.0;
   double cc_syy = 0.0;
   
   for ( int i1 = 0; i1 < 20; i1++ )
   {
    cc_sx  += hpp_prof1[i1];
    cc_sxx += ( hpp_prof1[i1] * hpp_prof1[i1] );
    cc_sy  += hpp_prof2[i1];
    cc_syy += ( hpp_prof2[i1] * hpp_prof2[i1] );
    cc_sxy += ( hpp_prof1[i1] * hpp_prof2[i1] );
   }
   
   double cc_mx = cc_sx / 20.0;
   double cc_my = cc_sy / 20.0;
   double cc_vx = ( cc_sxx / 20.0 ) - ( cc_mx * cc_mx );
   double cc_vy = ( cc_syy / 20.0 ) - ( cc_my * cc_my );
   double cc_sdx = sqrt(cc_vx);
   double cc_sdy = sqrt(cc_vy);
   double cc_cxy = ( cc_sxy / 20.0 ) - ( cc_mx * cc_my );
   
   scores.hpp_prf = cc_cxy / ( cc_sdx * cc_sdy );
   
   /* Binding probability */
   
   scores.bin_prb = pow( target1->getBindingProb(ip1) - target2->getBindingProb(ip2), 2.0);
   
   /* Neighbor distribution */
   
   vector<double> neighbor1;
   vector<double> neighbor2;
   
   target1->getNeighborDist(ip1, neighbor1);
   target2->getNeighborDist(ip2, neighbor2);
   
   scores.nbr_dis = FisherPitman(1000, neighbor1, neighbor2);
   
   /* Sequence entropy */
   
   double entropy1 = target1->getResEntropy(res1);
   double entropy2 = target2->getResEntropy(res2);
   
   scores.seq_ent = pow( entropy1 - entropy2, 2.0 );
   
   /* Magic distance */
   
   double magic_dst = 0.0;
   double magic_tot = 0.0;
   
   string ligA_id[MAXCMP];
   string ligB_id[MAXCMP];
   
   int ligA_nc = target1->getCmpsList( ligA_id );
   int ligB_nc = target2->getCmpsList( ligB_id );
   
   for ( int i1 = 0; i1 < ligA_nc; i1++ )
   {
    double ligA_xyz[MAXLIG][3];
    
    int ligA_ix = target1->getLigandNum( ligA_id[i1] );
    
    if ( ligA_ix > -1 )
    {
     int ligA_na = target1->getLigandCoords( ligA_ix, ligA_xyz ); ligA_na += 0;
     
     for ( int i2 = 0; i2 < ligB_nc; i2++ )
     {
      double ligB_xyz[MAXLIG][3];
      
      int ligB_ix = target2->getLigandNum( ligB_id[i2] );
      
      if ( ligB_ix > -1 )
      {
       int ligB_na = target2->getLigandCoords( ligB_ix, ligB_xyz ); ligB_na += 0;
       
       pair<int,int> k_match[MAXLIG];
       
       double s_match;
       
       int n_match = target1->getLigandMatch( ligA_ix, ligB_id[i2], k_match, s_match );
       
       if ( n_match > 0 )
       {
        double weights[MAXLIG];
        
        double mob_xyz[MAXLIG][3];
        double ref_xyz[MAXLIG][3];
        
        for ( int i3 = 0; i3 < n_match; i3++ )
        {
         weights[i3] = 1.0;
         
         for ( int i4 = 0; i4 < 3; i4++ )
         {
          ref_xyz[i3][i4] = ligA_xyz[k_match[i3].first-1][i4];
          mob_xyz[i3][i4] = ligB_xyz[k_match[i3].second-1][i4];
         }
        }
        
        double u[3][3];
        double t[3];
        double rms1 = 0.0;
        
        map< pair<int,int>, sup_data >::iterator it7;
        
        it7 = sup_matrix.find( std::make_pair(i1,i2) );
        
        if ( it7 != sup_matrix.end() )
        {
         rms1 = (*it7).second.sup_rms;
         
         for ( int i5 = 0; i5 < 3; i5++ )
         {
          t[i5] = (*it7).second.sup_t[i5];
          
          for ( int i6 = 0; i6 < 3; i6++ )
           u[i5][i6] = (*it7).second.sup_u[i5][i6];
         }
        }
        else
        {
         int mode = 1;
         int ier = 0;
         
         u3b_(&weights, &mob_xyz, &ref_xyz, &n_match, &mode, &rms1, &u, &t, &ier);
         
         rms1 = sqrt( rms1 / (double) n_match );
         
         sup_data tmp_sup;
         
         tmp_sup.sup_rms = rms1;
         
         for ( int i5 = 0; i5 < 3; i5++ )
         {
          tmp_sup.sup_t[i5] = t[i5];
          
          for ( int i6 = 0; i6 < 3; i6++ )
           tmp_sup.sup_u[i5][i6] = u[i5][i6];
         }
         
         sup_matrix.insert( pair< pair<int,int> ,sup_data >( std::make_pair(i1,i2), tmp_sup ) );
        }
        
        double res3_xyz[3];
        
        res3_xyz[0] = t[0] + u[0][0] * res2_xyz[0] + u[1][0] * res2_xyz[1] + u[2][0] * res2_xyz[2];
        res3_xyz[1] = t[1] + u[0][1] * res2_xyz[0] + u[1][1] * res2_xyz[1] + u[2][1] * res2_xyz[2];
        res3_xyz[2] = t[2] + u[0][2] * res2_xyz[0] + u[1][2] * res2_xyz[1] + u[2][2] * res2_xyz[2];
        
        magic_dst += ( s_match * s_match ) * sqrt( pow(res1_xyz[0]-res3_xyz[0], 2.0) + pow(res1_xyz[1]-res3_xyz[1], 2.0) + pow(res1_xyz[2]-res3_xyz[2], 2.0) );
        
        magic_tot += ( s_match * s_match );
       }
      }
     }
    }
   }
   
   magic_dst /= magic_tot;
   
   if ( magic_dst > 0.0 )
    scores.mag_dst = log( 1.0 / magic_dst );
   else
    scores.mag_dst = 0.0;
   
   aln_pair.insert( pair< pair<int,int>, pair_score > ( std::make_pair(res1,res2), scores ) );
  }
 }
}


// ==================================================================================   SitePairSVM

void SitePairSVM( map< pair<int,int>, pair_score > &a1_pair, ModelSVM * a1_model )
{
 map< pair<int,int>, pair_score >::iterator it1;
 
 for ( it1 = a1_pair.begin(); it1 != a1_pair.end(); it1++ )
 {
  double svm1[MAXSV1];
  
  svm1[0] = (it1->second).seq_prf;
  svm1[1] = (it1->second).sec_prf;
  svm1[2] = (it1->second).hpp_prf;
  svm1[3] = (it1->second).bin_prb;
  svm1[4] = (it1->second).nbr_dis;
  svm1[5] = (it1->second).seq_ent;
  svm1[6] = (it1->second).mag_dst;
  
  (it1->second).sco_svc = a1_model->SVMpredict( 1, svm1 );
  
  svm1[0] = (it1->second).seq_prf;
  svm1[1] = (it1->second).sec_prf;
  svm1[2] = (it1->second).hpp_prf;
  svm1[3] = (it1->second).bin_prb;
  svm1[4] = (it1->second).nbr_dis;
  svm1[5] = (it1->second).seq_ent;
  svm1[6] = (it1->second).mag_dst;
  
  (it1->second).sco_svr = a1_model->SVMpredict( 2, svm1 );
 }
}


// ==================================================================================   SitePairMunkres

double SitePairMunkres( map< pair<int,int>, pair_score > &a1_pair, map< pair<int,int>, pair_aligned > &a1_aln, int r1_tot, int r2_tot )
{
 a1_aln.clear();
 
 Matrix<double> a1_matrix(r1_tot, r2_tot);
 
 map< pair<int,int>, pair_score >::iterator it3;
 
 for (int it1 = 0; it1 < r1_tot; it1++)
  for (int it2 = 0; it2 < r2_tot; it2++)
   for ( it3 = a1_pair.begin(); it3 != a1_pair.end(); it3++ )
    if ( (it3->second).res1_index == it1 && (it3->second).res2_index == it2 )
     a1_matrix(it1,it2) = (it3->second).sco_svr;
 
 Munkres<double> a1_munkres;
 
 a1_munkres.solve(a1_matrix);
 
 double a1_tot = 0.0;
 
 for (int it1 = 0; it1 < r1_tot; it1++)
  for (int it2 = 0; it2 < r2_tot; it2++)
   if ( a1_matrix(it1,it2) > -0.9 )
    for ( it3 = a1_pair.begin(); it3 != a1_pair.end(); it3++ )
     if ( (it3->second).res1_index == it1 && (it3->second).res2_index == it2 )
     {
      pair_aligned tmp1 = { (it3->second).res1_index, (it3->second).res2_index, (it3->second).sco_svc, (it3->second).sco_svr, 0.0 };
      
      a1_aln.insert( pair< pair<int,int>, pair_aligned > ( std::make_pair( (it3->first).first, (it3->first).second ), tmp1 ) );
      
      a1_tot += (it3->second).sco_svr;
     }
 
 if ( a1_aln.size() > 0 )
  a1_tot /= (double) a1_aln.size();
 
 return a1_tot;
}


// ==================================================================================   SitePairScore

double SitePairScoreProb( map< pair<int,int>, pair_aligned > &a1_aln )
{
 double a1_prob = 0.0;
 int    a1_tot = 0.0;
 
 map< pair<int,int>, pair_aligned >::iterator it1;
 
 for ( it1 = a1_aln.begin(); it1 != a1_aln.end(); it1++ )
 {
  a1_prob += (it1->second).sco_svc;
  
  a1_tot++;
 }
 
 return sqrt( a1_prob / (double) a1_tot );
}


// ==================================================================================   SitePairTLscore

double SitePairTLscore( Protein * target1, Protein * target2 )
{
 double tls1 = 0.0;
 
 for ( int it1 = 0; it1 < 5; it1++ )
  if ( target1->getPocketProperty(it1*2+2) > 0.0 || target2->getPocketProperty(it1*2+2) > 0.0 )
   tls1 += fabs( target1->getPocketProperty(it1*2+1) - target2->getPocketProperty(it1*2+1) ) / sqrt( pow( target1->getPocketProperty(it1*2+2), 2.0 ) + pow( target2->getPocketProperty(it1*2+2), 2.0 ) );
 
 tls1 /= 5.0;
 
 return tls1;
}


// ==================================================================================   SitePairPMscore

double SitePairPMscore( Protein * target1, Protein * target2 )
{
 int sco1 = 0;
 int sco2 = 0;
 int sco3 = 0;
 
 for ( int it1 = 0; it1 < 5; it1++ )
  for ( int it2 = it1; it2 < 5; it2++ )
   for ( int it3 = 0; it3 < 3; it3++ )
    for ( int it4 = it3; it4 < 3; it4++ )
    {
     vector<double> binA;
     vector<double> binB;
     
     int nbinA = target1->getPMbin(binA, it1, it2, it3, it4);
     int nbinB = target2->getPMbin(binB, it1, it2, it3, it4);
     
     sco1 += nbinA;
     sco2 += nbinB;
     
     int i1 = 0;
     int i2 = 0;
     
     if ( nbinA > 0 && nbinB > 0 )
      while ( i1 < nbinA && i2 < nbinB )
      {
       if ( fabs( binA[i1] - binB[i2] ) <= 0.5 )
       {
        i1++;
        i2++;
        
        sco3++;
       }
       else
       {
        if ( binA[i1] < binB[i2] )
         i1++;
        else
         i2++;
       }
      }
    }
 
 double pms1 = 0.0;
 
 if ( sco1 > sco2 )
 {
  if ( sco1 > 0 )
   pms1 = ( (double) sco3 ) / ( (double) sco1 );
 }
 else
 {
  if ( sco2 > 0 )
   pms1 = ( (double) sco3 ) / ( (double) sco2 );
 }
 
 return pms1;
}
