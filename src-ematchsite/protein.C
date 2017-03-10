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


#include "protein.h"

using namespace std;

Protein::Protein( void )
{
 _protein_na = 0;
 _protein_nr = 0;
 
 _protein_xyz.clear();
 
 _protein_seq1.clear();
 
 _protein_aln1.clear();
 
 for ( int i1 = 0; i1 < MAXPRO; i1++ )
 {
  _protein_seq3[i1] = 0;
  
  for ( int i2 = 0; i2 < 20; i2++ )
  {
   _protein_prf1[i1][i2] = 0.0;
   _protein_prf3[i1][i2] = 0.0;
  }
  
  for ( int i2 = 0; i2 < 3; i2++ )
   _protein_prf2[i1][i2] = 0.0;
 }
 
 _pocket_number = 0;
 _ligclust_tot_smi = 0;
 _ligclust_tot_mac = 0;
 
 _protein_ntpl = 0;
 
 for ( int i1 = 0; i1 < 4; i1++ )
  _pocket_center[i1] = 0.0;
 
 _binres_tot = 0;
 
 _binding_res.clear();
 
 _ligand_fpts_smi.clear();
 _ligand_fpts_mac.clear();
 
 for ( int i1 = 0; i1 < MAXSMI; i1++ )
  _profile_fpts_smi[i1] = 0.0;
  
 for ( int i1 = 0; i1 < MAXMAC; i1++ )
  _profile_fpts_mac[i1] = 0.0;
 
 _ave_mw = 0.0;
 _ave_logp = 0.0;
 _ave_psa = 0.0;
 _ave_mr = 0.0;
 _ave_hba = 0.0;
 _ave_hbd = 0.0;
 
 _std_mw = 0.0;
 _std_logp = 0.0;
 _std_psa = 0.0;
 _std_mr = 0.0;
 _std_hba = 0.0;
 _std_hbd = 0.0;
 
 _pocket_plb = 0.0;
 
 _svm_pocket_confidence = 0.0;
 
 _cmps_nc = 0;
 
 for ( int i1 = 0; i1 < 3; i1++ )
 {
  _sup_t[i1] = 0.0;
  
  for ( int i2 = 0; i2 < 3; i2++ )
   _sup_u[i1][i2] = 0.0;
 }
}

Protein::~Protein() {}


// ==================================================================================   loadStructure

bool Protein::loadStructure( std::string p1_name )
{
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  if ( line1.size() > 53 )
   if ( line1.substr(0,6) == "ATOM  ")
   {
    string residue1 = line1.substr(17,3);
    string atom1    = line1.substr(12,4);
    string chain1   = line1.substr(21,1);
    
    int residue2 = atoi(line1.substr(22,4).c_str());
    int atom2    = atoi(line1.substr(6,5).c_str());
    
    double x1 = atof(line1.substr(30,8).c_str());
    double y1 = atof(line1.substr(38,8).c_str());
    double z1 = atof(line1.substr(46,8).c_str());
    
    _protein_xyz.push_back( CoordsProtein( atom2, residue2, x1, y1, z1, residue1, atom1, chain1 ) );
    
    _protein_na++;
    
    if ( atom1 == " CA " )
    {
     _protein_seq1.append( three2oneS( line1.substr(17,3) ) );
     
     double hp1[20];
     
     three2hpp( line1.substr(17,3), hp1 );
     
     for ( int i1 = 0; i1 < 20; i1++ )
      _protein_prf3[_protein_nr][i1] = hp1[i1];
     
     _protein_seq3[_protein_nr] = residue2;
     
     _protein_nr++;
    }
   }
 
 p1_file.close();
 
 strcpy(_protein_seq2, _protein_seq1.c_str());
 
 return EXIT_SUCCESS;
}


// ==================================================================================   loadPsipred

bool Protein::loadPsipred( std::string p1_name )
{
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  if ( line1.size() > 29 )
   if ( line1.find("PSIPRED") == string::npos )
   {
    int res1 = atoi(line1.substr(0,4).c_str());
    
    _protein_prf2[res1-1][0] = atof(line1.substr(9,7).c_str());
    _protein_prf2[res1-1][1] = atof(line1.substr(16,7).c_str());
    _protein_prf2[res1-1][2] = atof(line1.substr(23,7).c_str());
   }
 
 p1_file.close();
 
 return EXIT_SUCCESS;
}


// ==================================================================================   loadSequence

bool Protein::loadSequence( std::string p1_name )
{
 for ( int ii1 = 0; ii1 < _protein_nr; ii1++ )
  for ( int ii2 = 0; ii2 < 20; ii2++ )
   _protein_prf1[ii1][ii2] = 0.0;
 
 string line1;
 
 int ii3 = -2;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() ) { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  if ( ++ii3 > -1 )
  {
   istringstream iss1(line1);
   
   std::string iss2;
   
   int iss3 = -2;
   
   while ( iss1 && iss3++ < 19 )
   {
    iss1 >> iss2;
    
    if ( iss3 > -1 )
     _protein_prf1[ii3][iss3] = atof(iss2.c_str());
   }
  }
 
 p1_file.close();
 
 return EXIT_SUCCESS;
}


// ==================================================================================   loadPocket

bool Protein::loadPocket( std::string p1_name, int p1_number )
{
 string line1;
 
 _binres_tot = 0;
 
 bool pkt1 = false;
 bool pkt2 = false;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
 {
  if ( line1.size() > 31 )
   if ( line1.substr(0,6) == "POCKET" )
    if ( atoi(line1.substr(7,4).c_str()) == p1_number )
     pkt1 = true;
  
  if ( line1 == "TER" )
   pkt1 = false;
  
  if ( pkt1 )
  {
   if ( line1.size() > 31 )
    if ( line1.substr(0,6) == "POCKET" )
    {
     _pocket_number = atoi(line1.substr(7,4).c_str());
     
     _svm_pocket_confidence = atof(line1.substr(24,8).c_str());
    }
   
   if ( line1.size() > 41 )
    if ( line1.substr(0,6) == "CENTER" )
    {
     _pocket_center[0] = atof(line1.substr(6,9).c_str());
     _pocket_center[1] = atof(line1.substr(15,9).c_str());
     _pocket_center[2] = atof(line1.substr(24,9).c_str());
     _pocket_center[3] = atof(line1.substr(33,9).c_str());
    }
   
   if ( line1.size() > 1051 )
    if ( line1.substr(0,6) == "SMILES" )
    {
     std::string lig1 = line1.substr(28,MAXSMI).c_str();
     
     bitset<MAXSMI> lig2;
     
     string::iterator pkt3;
     int pkt4;
     
     for ( pkt3 = lig1.begin(), pkt4 = 0; pkt3 < lig1.end(); pkt3++, pkt4++ )
      switch(*pkt3)
      {
       case '0': lig2[MAXSMI-1-pkt4] = 0; break;
       case '1': lig2[MAXSMI-1-pkt4] = 1; break;
      }
     
     ligand_fpt_smi lig3;
     
     lig3.ftp_number = atoi(line1.substr(6,5).c_str());
     lig3.fpt_fraction = atof(line1.substr(19,8).c_str());
     lig3.fpt_pdbid = line1.substr(12,7);
     lig3.fpt_fingerprint = lig2;
     
     _ligand_fpts_smi.push_back( lig3 );
     
     _ligclust_tot_smi++;
    }
   
   if ( line1.size() > 195 )
    if ( line1.substr(0,5) == "MACCS" )
    {
     std::string lig1 = line1.substr(28,MAXMAC).c_str();
     
     bitset<MAXMAC> lig2;
     
     string::iterator pkt3;
     int pkt4;
     
     for ( pkt3 = lig1.begin(), pkt4 = 0; pkt3 < lig1.end(); pkt3++, pkt4++ )
      switch(*pkt3)
      {
       case '0': lig2[MAXMAC-1-pkt4] = 0; break;
       case '1': lig2[MAXMAC-1-pkt4] = 1; break;
      }
     
     ligand_fpt_mac lig3;
     
     lig3.ftp_number = atoi(line1.substr(6,5).c_str());
     lig3.fpt_fraction = atof(line1.substr(19,8).c_str());
     lig3.fpt_pdbid = line1.substr(12,7);
     lig3.fpt_fingerprint = lig2;
     
     _ligand_fpts_mac.push_back( lig3 );
     
     _ligclust_tot_mac++;
    }
   
   if ( line1.size() > 5126 )
    if ( line1.substr(0,7) == "PROFSMI" )
     for ( int prf1 = 0; prf1 < MAXSMI; prf1++ )
      _profile_fpts_smi[prf1] = atof(line1.substr(7+prf1*5,5).c_str());
   
   if ( line1.size() > 846 )
    if ( line1.substr(0,7) == "PROFMAC" )
     for ( int prf1 = 0; prf1 < MAXMAC; prf1++ )
      _profile_fpts_mac[prf1] = atof(line1.substr(7+prf1*5,5).c_str());
   
   if ( line1.size() > 29 )
   {
         if ( line1.substr(0,9) == "PROP_MW  " )
    {
     _ave_mw = atof(line1.substr(10,10).c_str());
     _std_mw = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_LOGP" )
    {
     _ave_logp = atof(line1.substr(10,10).c_str());
     _std_logp = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_PSA " )
    {
     _ave_psa = atof(line1.substr(10,10).c_str());
     _std_psa = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_MR  " )
    {
     _ave_mr = atof(line1.substr(10,10).c_str());
     _std_mr = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_HBD " )
    {
     _ave_hbd = atof(line1.substr(10,10).c_str());
     _std_hbd = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_HBA " )
    {
     _ave_hba = atof(line1.substr(10,10).c_str());
     _std_hba = atof(line1.substr(20,10).c_str());
    }
   }
   
   if ( line1.size() > 13 )
    if ( line1.substr(0,8) == "PLBINDEX" )
     _pocket_plb = atof(line1.substr(8,6).c_str());
   
   if ( line1.size() > 46 )
    if ( line1.substr(0,7) == "RESIDUE" )
    {
     binding_residue res1;
     
     int res2 = atoi(line1.substr(7,6).c_str());
     
     bool res3 = false;
     
     if ( line1.substr(16,1) == "*" )
     {
      res3 = true;
      
      res1.residue_number = res2;
      res1.residue_code = line1.substr(14,1);
      res1.residue_binding = res3;
      res1.residue_score = atof(line1.substr(39,8).c_str());
      res1.residue_distance = atof(line1.substr(22,9).c_str());
      res1.residue_fraction = atof(line1.substr(31,8).c_str());
      res1.residue_group = one2group(res1.residue_code);
      
      vector<CoordsProtein>::iterator ipca1;
      
      for ( int ipca2 = 0; ipca2 < 3; ipca2++ )
       res1.residue_coordsSC[ipca2] = 0.0;
      
      int tmp_nsc = 0;
      
      for ( ipca1 = _protein_xyz.begin(); ipca1 < _protein_xyz.end(); ipca1++ )
       if ( (*ipca1).getResidueNumber() == res2 )
       {
        if ( (*ipca1).getAtomName() == " CA " )
         for ( int ipca2 = 0; ipca2 < 3; ipca2++ )
          res1.residue_coordsCA[ipca2] = (*ipca1).getCoords(ipca2+1);
        
        if ( (*ipca1).getAtomName() == " CB " || ( (*ipca1).getAtomName() == " CA " && res1.residue_code == "G" ) )
         for ( int ipca2 = 0; ipca2 < 3; ipca2++ )
          res1.residue_coordsCB[ipca2] = (*ipca1).getCoords(ipca2+1);
        
        for ( int ipca2 = 0; ipca2 < 3; ipca2++ )
         res1.residue_coordsSC[ipca2] += (*ipca1).getCoords(ipca2+1);
        
        tmp_nsc++;
       }
      
      for ( int ipca2 = 0; ipca2 < 3; ipca2++ )
       res1.residue_coordsSC[ipca2] /= ( (double) tmp_nsc );
      
      _binding_res.insert( pair<int,binding_residue>(_binres_tot,res1) );
      
      if ( res3 )
       _binres_tot++;
     }
    }
   
   pkt2 = true;
  }
 }
 
 p1_file.close();
 
 if ( pkt2 )
  return EXIT_SUCCESS;
 else
  return EXIT_FAILURE;
}


// ==================================================================================   loadLigands

bool Protein::loadLigands( std::string c1_name, int c1_number )
{
 string line1;
 
 string tmp_sdf_t[MAXSDF];
 int    tmp_sdf_n = 0;
 
 ifstream c1_file( c1_name.c_str() );
 
 if ( !c1_file.is_open() ) { cout << "Cannot open " << c1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(c1_file,line1))
 {
  tmp_sdf_t[tmp_sdf_n++] = line1;
  
  if ( line1.compare("$$$$") == 0 )
  {
   if ( _cmps_nc < MAXCMP )
   {
    bool pkt1 = false;
    
    _cmps_nla[_cmps_nc] = atoi(tmp_sdf_t[3].substr(0,3).c_str());
    _cmps_nlb[_cmps_nc] = atoi(tmp_sdf_t[3].substr(3,3).c_str());
    
    for ( int i1 = 4 + _cmps_nla[_cmps_nc] + _cmps_nlb[_cmps_nc]; i1 < tmp_sdf_n - 1; i1++ )
     if ( tmp_sdf_t[i1].find("EFINDSITE_POCKET") != string::npos )
      if ( atoi(tmp_sdf_t[i1+1].c_str()) == c1_number )
       pkt1 = true;
     
    if ( pkt1 )
    {
     bool w1 = false;
     bool w2 = false;
     bool w3 = false;
     bool w4 = false;
     bool w5 = false;
     bool w6 = false;
     bool w7 = false;
     bool w8 = false;
     bool w9 = false;
     
     for ( int i1 = 4 + _cmps_nla[_cmps_nc] + _cmps_nlb[_cmps_nc]; i1 < tmp_sdf_n - 1; i1++ )
     {
      if ( tmp_sdf_t[i1].find("FINGERPRINT") != string::npos )
      {
       string::iterator i2;
       int i3;
       
       for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
        switch(*i2)
        {
         case '0': _cmps_fpt_smi[_cmps_nc][MAXSMI-1-i3] = 0; break;
         case '1': _cmps_fpt_smi[_cmps_nc][MAXSMI-1-i3] = 1; break;
        }
       
       w1 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("MACCS166") != string::npos )
      {
       string::iterator i2;
       int i3;
       
       for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
        switch(*i2)
        {
         case '0': _cmps_fpt_mac[_cmps_nc][MAXMAC-1-i3] = 0; break;
         case '1': _cmps_fpt_mac[_cmps_nc][MAXMAC-1-i3] = 1; break;
        }
       
       w2 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("OB_MW") != string::npos )
      {
       _cmps_mw[_cmps_nc] = atof(tmp_sdf_t[i1+1].c_str());
       
       w3 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("OB_logP") != string::npos )
      {
       _cmps_logp[_cmps_nc] = atof(tmp_sdf_t[i1+1].c_str());
       
       w4 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("OB_PSA") != string::npos )
      {
       _cmps_psa[_cmps_nc] = atof(tmp_sdf_t[i1+1].c_str());
       
       w5 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("OB_MR") != string::npos )
      {
       _cmps_mr[_cmps_nc] = atof(tmp_sdf_t[i1+1].c_str());
       
       w6 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("MCT_HBD") != string::npos )
      {
       _cmps_hbd[_cmps_nc] = atoi(tmp_sdf_t[i1+1].c_str());
       
       w7 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("MCT_HBA") != string::npos )
      {
       _cmps_hba[_cmps_nc] = atoi(tmp_sdf_t[i1+1].c_str());
       
       w8 = true;
      }
      
      else if ( tmp_sdf_t[i1].find("MOLID") != string::npos )
      {
       _cmps_id[_cmps_nc] = tmp_sdf_t[i1+1].c_str();
       
       w9 = true;
      }
     }
     
     for ( int i1 = 0; i1 < 3; i1++ )
      _cmps_cen[_cmps_nc][i1] = 0.0;
     
     _cmps_xyz[_cmps_nc].clear();
     
     for ( int i1 = 4; i1 < _cmps_nla[_cmps_nc] + 4; i1++ )
     {
      double l1_x = atof(tmp_sdf_t[i1].substr(0,10).c_str());
      double l1_y = atof(tmp_sdf_t[i1].substr(10,10).c_str());
      double l1_z = atof(tmp_sdf_t[i1].substr(20,10).c_str());
      
      _cmps_cen[_cmps_nc][0] += l1_x;
      _cmps_cen[_cmps_nc][1] += l1_y;
      _cmps_cen[_cmps_nc][2] += l1_z;
      
      _cmps_xyz[_cmps_nc].push_back( CoordsLigand( i1-3, l1_x, l1_y, l1_z, tmp_sdf_t[i1].substr(30,3) ) );
     }
     
     for ( int i1 = _cmps_nla[_cmps_nc] + 4; i1 < _cmps_nla[_cmps_nc] + 4 + _cmps_nlb[_cmps_nc]; i1++ )
     {
      ligand_bond t_bnd;
      
      t_bnd.atom1 = atoi(tmp_sdf_t[i1].substr(0,3).c_str());
      t_bnd.atom2 = atoi(tmp_sdf_t[i1].substr(3,3).c_str());
      t_bnd.bond_type = atoi(tmp_sdf_t[i1].substr(6,3).c_str());
      
      _cmps_bnd[_cmps_nc].push_back(t_bnd);
     }
     
     for ( int i1 = 0; i1 < 3; i1++ )
      _cmps_cen[_cmps_nc][i1] /= _cmps_nla[_cmps_nc];
     
     _cmps_matching[_cmps_nc].clear();
     
     if ( w1 && w2 && w3 && w4 && w5 && w6 && w7 && w8 && w9 )
      _cmps_nc++;
    }
   }
   
   tmp_sdf_n = 0;
  }
 }
 
 c1_file.close();
 
 if ( _cmps_nc > 0 )
  return EXIT_SUCCESS;
 else
  return EXIT_FAILURE;
}


// ==================================================================================   loadAlignments

bool Protein::loadAlignments( std::string a1_name )
{
 _protein_ntpl = 0;
 
 string line1;
 
 string ln_inf;
 string ln_tar;
 string ln_dot;
 string ln_tpl;
 
 ifstream a1_file( a1_name.c_str() );
 
 if ( !a1_file.is_open() )  { cout << "Cannot open " << a1_name << endl; exit(EXIT_FAILURE); }
 
 while ( getline( a1_file, line1 ) )
 {
       if ( ln_inf.size() < 1 )
  {
   ln_inf = line1;
  }
  else if ( ln_tar.size() < 1 )
  {
   ln_tar = line1;
  }
  else if ( ln_dot.size() < 1 )
  {
   ln_dot = line1;
  }
  else if ( ln_tpl.size() < 1 )
  {
   ln_tpl = line1;
  }
  
  if ( line1 == "*" )
  {
   std::string tpl_nam;
   
   tpl_alignment tmp_al;
   
   istringstream iss1(ln_inf);
   
   std::string iss2;
   
   int iss3 = 0;
   
   while (iss1)
   {
    iss1 >> iss2;
    
    switch(iss3)
    {
     case 0 : tpl_nam = iss2.substr(1, iss2.size()-1); break;
     case 1 : tmp_al.tpl_res = atoi(iss2.c_str()); break;
     case 2 : tmp_al.tpl_len = atoi(iss2.c_str()); break;
     case 3 : tmp_al.tpl_tms = atof(iss2.c_str()); break;
     case 4 : tmp_al.tpl_rms = atof(iss2.c_str()); break;
     case 5 : tmp_al.tpl_sid = atof(iss2.c_str()); break;
    }
    
    iss3++;
   }
   
   tmp_al.tpl_pairs.clear();
   
   int rn1 = 0;
   int rn2 = 0;
   
   for ( int i1 = 0; i1 < (int) ln_tar.size(); i1++ )
   {
    if ( ln_tar.substr(i1, 1) != "-" )
     rn1++;
    if ( ln_tpl.substr(i1, 1) != "-" )
     rn2++;
    
    if ( ln_tar.substr(i1, 1) != "-" && ln_tpl.substr(i1, 1) != "-" )
     tmp_al.tpl_pairs.push_back( std::make_pair(rn1,ln_tpl.substr(i1, 1)) );
   }
   
   _protein_aln1.insert( pair<string,tpl_alignment>(tpl_nam, tmp_al) );
   
   _protein_ntpl++;
   
   ln_inf.clear();
   ln_tar.clear();
   ln_dot.clear();
   ln_tpl.clear();
  }
 }
 
 a1_file.close();
 
 if ( _protein_ntpl > 0 )
  return EXIT_SUCCESS;
 else
  return EXIT_FAILURE;
}


// ==================================================================================   loadLibrary

bool Protein::loadLibrary( std::string m1_lib, int m1_tot, std::string m1_ids[], int & m1_pairs )
{
 map<string,bool> m1_check;
 
 for ( int it1 = 0; it1 < m1_tot; it1++ )
  m1_check.insert( pair<string,bool>(m1_ids[it1], true) );
 
 map<string,bool> m1_list;
 
 string line1;
 
 ifstream m1_file( (m1_lib+"/cmps.lst").c_str() );
 
 if ( !m1_file.is_open() )  { cout << "Cannot open " << m1_lib << "/cmps.lst" << endl; exit(EXIT_FAILURE); }
 
 while (getline(m1_file,line1))
  m1_list.insert( pair<string,bool>(line1, true) );
 
 m1_file.close();
 
 for ( int it1 = 0; it1 < _cmps_nc; it1++ )
 {
  if ( m1_list.count(_cmps_id[it1]) )
  {
   igzstream m1_data( (m1_lib+"/data/"+_cmps_id[it1].substr(1,2)+"/"+_cmps_id[it1]).c_str() );
   
   if ( ! m1_data.good() )
   {
    std::cerr << "ERROR: Opening file `" << m1_lib << "/data/" << _cmps_id[it1].substr(1,2) << "/" << _cmps_id[it1] << "' failed\n";
    
    return EXIT_FAILURE;
   }
   
   while ( getline(m1_data,line1) )
    if ( line1.length() > 3 )
     if ( line1.compare(0, 3, "ALN") == 0 )
     {
      string tmp_id = line1.substr(4,7);
      double tmp_tc = atof(line1.substr(11,6).c_str());
      int    tmp_np = atoi(line1.substr(17,4).c_str());
      
      if ( tmp_np >= 3 && m1_list.count(tmp_id) && m1_check.count(tmp_id) )
      {
       compound_match tmp_pr;
       
       tmp_pr.tanimoto = tmp_tc;
       tmp_pr.total_pairs = tmp_np;
       
       istringstream iss1(line1.substr(22,line1.length()-22));
       
       std::string iss2;
       
       int iss3 = 0;
       int iss4 = 0;
       
       while (iss1)
       {
        iss1 >> iss2;
        
        if ( iss3 == 0 )
         iss3 = atoi(iss2.c_str());
        else
        {
         iss4 = atoi(iss2.c_str());
         
         tmp_pr.atom_pairs.push_back( std::make_pair(iss3,iss4) );
         
         iss3 = 0;
         iss4 = 0;
        }
       }
       
       (_cmps_matching[it1]).insert( pair<string,compound_match>(tmp_id, tmp_pr) );
       
       m1_pairs++;
      }
     }
   
   m1_data.close();
  }
 }
 
 return EXIT_SUCCESS;
}


// ==================================================================================   getTanimotoSMILES

double Protein::getTanimotoSMILES( std::bitset<MAXSMI> cmp_fpt, std::string cmp_method )
{
 double score1 = 0.0;
 
 list<ligand_fpt_smi>::iterator i1;
 
 if ( cmp_method == "T" )
  for ( i1 = _ligand_fpts_smi.begin() ; i1 != _ligand_fpts_smi.end(); i1++ )
   score1 += (*i1).fpt_fraction * getTanimoto1024( cmp_fpt, (*i1).fpt_fingerprint );
 else if ( cmp_method == "A" )
  for ( i1 = _ligand_fpts_smi.begin() ; i1 != _ligand_fpts_smi.end(); i1++ )
   score1 += (*i1).fpt_fraction * getTanimotoAve1024( cmp_fpt, (*i1).fpt_fingerprint );
 else if ( cmp_method == "C" )
  score1 = getTanimotoCnt1024( cmp_fpt, _profile_fpts_smi );
 
 return score1;
}


// ==================================================================================   getTanimotoMACCS

double Protein::getTanimotoMACCS( std::bitset<MAXMAC> cmp_fpt, std::string cmp_method )
{
 double score1 = 0.0;
 
 list<ligand_fpt_mac>::iterator i1;
 
 if ( cmp_method == "T" )
  for ( i1 = _ligand_fpts_mac.begin() ; i1 != _ligand_fpts_mac.end(); i1++ )
   score1 += (*i1).fpt_fraction * getTanimoto166( cmp_fpt, (*i1).fpt_fingerprint );
 else if ( cmp_method == "A" )
  for ( i1 = _ligand_fpts_mac.begin() ; i1 != _ligand_fpts_mac.end(); i1++ )
   score1 += (*i1).fpt_fraction * getTanimotoAve166( cmp_fpt, (*i1).fpt_fingerprint );
 else if ( cmp_method == "C" )
  score1 = getTanimotoCnt166( cmp_fpt, _profile_fpts_mac );
 
 return score1;
}


// ==================================================================================   getPocketProperty

double Protein::getPocketProperty( int prop1 )
{
 switch (prop1)
 {
  case  1 : return _ave_mw;
  case  2 : return _std_mw;
  case  3 : return _ave_logp;
  case  4 : return _std_logp;
  case  5 : return _ave_psa;
  case  6 : return _std_psa;
  case  7 : return _ave_hbd;
  case  8 : return _std_hbd;
  case  9 : return _ave_hba;
  case 10 : return _std_hba;
  case 11 : return _ave_mr;
  case 12 : return _std_mr;
  
  default : return 0;
 }
}


// ==================================================================================   matchLigands

int Protein::matchLigands( std::string m1_name, int m1_tot, std::string m1_ids[], bool m1_order )
{
 int m1_pairs = 0;
 
 map<string,bool> m1_check;
 
 for ( int it1 = 0; it1 < m1_tot; it1++ )
  m1_check.insert( pair<string,bool>(m1_ids[it1], true) );
 
 vector<string> kdata;
 
 string line1;
 
 if ( (m1_name.substr(m1_name.length()-3,3)).compare(0, 3, ".gz") == 0 )
 {
  igzstream m1_file( m1_name.c_str() );
  
  if ( ! m1_file.good() ) { cout << "Cannot open " << m1_name << endl; exit(EXIT_FAILURE); }
  
  while ( getline( m1_file, line1 ) )
   if ( line1.length() > 3 )
    if ( line1.compare(0, 3, "ALN") == 0 )
     kdata.push_back( line1 );
  
  m1_file.close();
 }
 else
 {
  ifstream m1_file( m1_name.c_str() );
  
  if ( !m1_file.is_open() )  { cout << "Cannot open " << m1_name << endl; exit(EXIT_FAILURE); }
  
  while ( getline( m1_file, line1 ) )
   if ( line1.length() > 3 )
    if ( line1.compare(0, 3, "ALN") == 0 )
     kdata.push_back( line1 );
  
  m1_file.close();
 }
 
 vector<string>::iterator line2;
 
 for ( line2 = kdata.begin(); line2 < kdata.end(); line2++ )
 {
  string tmp_id1;
  string tmp_id2;
  
  if ( m1_order )
  {
   tmp_id1 = (*line2).substr(4,7);
   tmp_id2 = (*line2).substr(12,7);
  }
  else
  {
   tmp_id1 = (*line2).substr(12,7);
   tmp_id2 = (*line2).substr(4,7);
  }
  
  double tmp_tc = atof((*line2).substr(19,6).c_str());
  int    tmp_np = atoi((*line2).substr(25,4).c_str());
  
  if ( tmp_np >= 3 && m1_check.count(tmp_id2) )
   for ( int it1 = 0; it1 < _cmps_nc; it1++ )
    if ( _cmps_id[it1].compare(tmp_id1) == 0)
    {
     compound_match tmp_pr;
     
     tmp_pr.tanimoto = tmp_tc;
     tmp_pr.total_pairs = tmp_np;
     
     istringstream iss1((*line2).substr(30,(*line2).length()-30));
     
     std::string iss2;
     
     int iss3 = 0;
     int iss4 = 0;
     
     while (iss1)
     {
      iss1 >> iss2;
      
      if ( iss3 == 0 )
       iss3 = atoi(iss2.c_str());
      else
      {
       iss4 = atoi(iss2.c_str());
       
       if ( m1_order )
        tmp_pr.atom_pairs.push_back( std::make_pair(iss3,iss4) );
       else
        tmp_pr.atom_pairs.push_back( std::make_pair(iss4,iss3) );
       
       iss3 = 0;
       iss4 = 0;
      }
     }
     
     (_cmps_matching[it1]).insert( pair<string,compound_match>(tmp_id2, tmp_pr) );
     
     m1_pairs++;
    }
 }
 
 return m1_pairs;
}


// ==================================================================================   getCmpsTotal

int Protein::getCmpsTotal( void )
{
 return _cmps_nc;
}


// ==================================================================================   getCmpsList

int Protein::getCmpsList( std::string tab1[] )
{
 for ( int i1 = 0; i1 < _cmps_nc; i1++ )
  tab1[i1] = _cmps_id[i1];
 
 return _cmps_nc;
}


// ==================================================================================   getLigandName

std::string Protein::getLigandName( int l1_number )
{
 return _cmps_id[l1_number];
}


// ==================================================================================   getLigandNum

int Protein::getLigandNum( std::string l1_name )
{
 int lig1 = _cmps_nc;
 
 for ( int it1 = 0; it1 < _cmps_nc; it1++ )
  if ( _cmps_id[it1] == l1_name )
  {
   lig1 = it1;
   
   break;
  }
 
 if ( lig1 == _cmps_nc )
  return -1;
 
 return lig1;
}


// ==================================================================================   getLigandCoords

int Protein::getLigandCoords( int l1_number, double l1_coords[][3] )
{
 int lig1 = 0;
 
 vector<CoordsLigand>::iterator it1;
 
 for ( it1 = _cmps_xyz[l1_number].begin(); it1 < _cmps_xyz[l1_number].end(); it1++ )
 {
  for ( int it2 = 0; it2 < 3; it2++ )
   l1_coords[lig1][it2] = (*it1).getCoords(it2+1);
  
  lig1++;
 }
 
 return lig1;
}


// ==================================================================================   getLigandAtomNames

int Protein::getLigandAtomNames( int l1_number, string l1_names[] )
{
 int lig1 = 0;
 
 vector<CoordsLigand>::iterator it2;
 
 for ( it2 = _cmps_xyz[l1_number].begin(); it2 < _cmps_xyz[l1_number].end(); it2++ )
  l1_names[lig1++] = (*it2).getAtomName();
 
 return lig1;
}


// ==================================================================================   getLigandBonds

int Protein::getLigandBonds( int l1_number, int l1_bonds[][3] )
{
 int lig1 = 0;
 
 vector<ligand_bond>::iterator it1;
 
 for ( it1 = _cmps_bnd[l1_number].begin(); it1 < _cmps_bnd[l1_number].end(); it1++ )
 {
  l1_bonds[lig1][0] = (*it1).atom1;
  l1_bonds[lig1][1] = (*it1).atom2;
  l1_bonds[lig1][2] = (*it1).bond_type;
  
  lig1++;
 }
 
 return lig1;
}


// ==================================================================================   getBindingResiduesTotal

int Protein::getBindingResiduesTotal( void )
{
 return _binres_tot;
}


// ==================================================================================   getBindingResCoords

void Protein::getBindingResCoords( int r1_number, double r1_coords[], int r1_type )
{
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 for ( int it2 = 0; it2 < 3; it2++ )
 {
       if ( r1_type == 0 )
   r1_coords[it2] = (*it1).second.residue_coordsCA[it2];
  else if ( r1_type == 1 )
   r1_coords[it2] = (*it1).second.residue_coordsCB[it2];
  else if ( r1_type == 2 )
   r1_coords[it2] = (*it1).second.residue_coordsSC[it2];
 }
}


// ==================================================================================   getResidueCoords

void Protein::getResidueCoords( int r1_number, vector<std::string> &r1_lines )
{
 vector<CoordsProtein>::iterator it1;
 
 for ( it1 = _protein_xyz.begin(); it1 < _protein_xyz.end(); it1++ )
  if ( (*it1).getResidueNumber() == r1_number )
  {
   double tmp_xyz[3];
   
   for ( int it2 = 0; it2 < 3; it2++ )
    tmp_xyz[it2] = (*it1).getCoords(it2+1);
   
   std::ostringstream tmp_line;
   
   tmp_line << "REMARK RESA " << setw(5) << (*it1).getAtomNumber()
                              << setw(5) << (*it1).getAtomName()
                              << setw(4) << (*it1).getResidueName()
                              << " A"
                              << setw(4) << (*it1).getResidueNumber()
                              << fixed << setw(12) << setprecision(3) << tmp_xyz[0]
                              << fixed << setw(8)  << setprecision(3) << tmp_xyz[1]
                              << fixed << setw(8)  << setprecision(3) << tmp_xyz[2];
   
   r1_lines.push_back( tmp_line.str() );
  }
}


// ==================================================================================   getBindingResNum

int Protein::getBindingResNum( int r1_number )
{
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 return (*it1).second.residue_number;
}


// ==================================================================================   getBindingResCode

std::string Protein::getBindingResCode( int r1_number )
{
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 return (*it1).second.residue_code;
}


// ==================================================================================   getBindingResGroup

int Protein::getBindingResGroup( int r1_number )
{
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 return (*it1).second.residue_group;
}


// ==================================================================================   getSeqProfile

void Protein::getSeqProfile( int r1_number, double r1_profile[] )
{
 for ( int it1 = 0; it1 < 20; it1++ )
  r1_profile[it1] = _protein_prf1[r1_number][it1];
}


// ==================================================================================   getSecProfile

void Protein::getSecProfile( int r1_number, double r1_profile[] )
{
 for ( int it1 = 0; it1 < 3; it1++ )
  r1_profile[it1] = _protein_prf2[r1_number][it1];
}


// ==================================================================================   getHPProfile

void Protein::getHPProfile( int r1_number, double r1_profile[] )
{
 for ( int it1 = 0; it1 < 20; it1++ )
  r1_profile[it1] = _protein_prf3[r1_number][it1];
}


// ==================================================================================   getBindingProb

double Protein::getBindingProb( int r1_number )
{
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 return (*it1).second.residue_score;
}


// ==================================================================================   getNeighborDist

void Protein::getNeighborDist( int r1_number, vector<double> &r1_distances )
{
 r1_distances.clear();
 
 map<int,binding_residue>::iterator it1;
 
 it1 = _binding_res.find( r1_number );
 
 for ( int it2 = 0; it2 < _binres_tot; it2++ )
  if ( it2 != r1_number )
  {
   map<int,binding_residue>::iterator it3;
   
   it3 = _binding_res.find( it2 );
   
   double rr1 = sqrt( pow( (*it1).second.residue_coordsCA[0] - (*it3).second.residue_coordsCA[0], 2.0) + 
                      pow( (*it1).second.residue_coordsCA[1] - (*it3).second.residue_coordsCA[1], 2.0) + 
                      pow( (*it1).second.residue_coordsCA[2] - (*it3).second.residue_coordsCA[2], 2.0) );
   
   r1_distances.push_back( rr1 );
  }
 
 std::sort(r1_distances.begin(), r1_distances.end(), std::less<double>());
}


// ==================================================================================   getResEntropy

double Protein::getResEntropy( int r1_number )
{
 map<string,tpl_alignment>::iterator it1;
 
 int freq1[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
 int freq2 = 0;
 
 for ( it1 = _protein_aln1.begin(); it1 != _protein_aln1.end(); it1++ )
 {
  vector< pair<int,string> >::iterator it2;
  
  for ( it2 = (it1->second).tpl_pairs.begin(); it2 < (it1->second).tpl_pairs.end(); it2++ )
   if ( (*it2).first == r1_number )
   {
    freq1[one2num((*it2).second[0])]++;
    
    freq2++;
   }
 }
 
 double ent1 = 0.0;
 
 if ( freq2 )
  for ( int it2 = 0; it2 < 20; it2++ )
  {
   double ent2 = ( (double) freq1[it2] ) / ( (double) freq2 );
   
   if ( ent2 > 0.0 )
    ent1 += ent2 * log2(ent2);
  }
 
 ent1 *= ( -1.0 );
 
 return ent1;
}


// ==================================================================================   getLigandMatch

int Protein::getLigandMatch( int l1_number, std::string l1_name, pair<int,int> l1_match[], double & l1_score )
{
 map<string,compound_match>::iterator it1;
 
 it1 = _cmps_matching[l1_number].find( l1_name );
 
 int it3 = 0;
 
 if ( it1 != _cmps_matching[l1_number].end() )
 {
  vector< pair<int,int> >::iterator it2;
  
  for ( it2 = (*it1).second.atom_pairs.begin(); it2 < (*it1).second.atom_pairs.end(); it2++ )
   l1_match[it3++] = (*it2);
  
  l1_score = (*it1).second.tanimoto;
 }
 
 return it3;
}


// ==================================================================================   alignmentSuperpose

double Protein::alignmentSuperpose( map< pair<int,int>, pair_aligned > &sup_aln, Protein * sup_ref )
{
 double mob_xyz[MAXLIG][3];
 double ref_xyz[MAXLIG][3];
 
 int pairs_tot = 0;
 
 map< pair<int,int>, pair_aligned >::iterator it1;
 
 map<int,binding_residue>::iterator it4;
 
 for ( it1 = sup_aln.begin(); it1 != sup_aln.end(); it1++ )
 {
  double tmp_xyz[3];
  
  sup_ref->getBindingResCoords( (it1->second).res1_index, tmp_xyz, 0 );
  
  for ( int it3 = 0; it3 < 3; it3++ )
   ref_xyz[pairs_tot][it3] = tmp_xyz[it3];
  
  it4 = _binding_res.find( (it1->second).res2_index );
  
  for ( int it3 = 0; it3 < 3; it3++ )
   mob_xyz[pairs_tot][it3] = (*it4).second.residue_coordsCA[it3];
  
  pairs_tot++;
 }
 
 double weights[MAXLIG];
 
 for ( int it2 = 0; it2 < pairs_tot; it2++ )
  weights[it2] = 1.0;
 
 int mode = 1;
 double rms1 = 0.0;
 int ier = 0;
 
 u3b_(&weights, &mob_xyz, &ref_xyz, &pairs_tot, &mode, &rms1, &_sup_u, &_sup_t, &ier);
 
 rms1 = sqrt( rms1 / (double) pairs_tot );
 
 for ( it4 = _binding_res.begin(); it4 != _binding_res.end(); it4++ )
 {
  double tmp_xyz[3];
  
  tmp_xyz[0] = _sup_t[0] + _sup_u[0][0] * (it4->second).residue_coordsCA[0] + _sup_u[1][0] * (it4->second).residue_coordsCA[1] + _sup_u[2][0] * (it4->second).residue_coordsCA[2];
  tmp_xyz[1] = _sup_t[1] + _sup_u[0][1] * (it4->second).residue_coordsCA[0] + _sup_u[1][1] * (it4->second).residue_coordsCA[1] + _sup_u[2][1] * (it4->second).residue_coordsCA[2];
  tmp_xyz[2] = _sup_t[2] + _sup_u[0][2] * (it4->second).residue_coordsCA[0] + _sup_u[1][2] * (it4->second).residue_coordsCA[1] + _sup_u[2][2] * (it4->second).residue_coordsCA[2];
  
  for ( int it3 = 0; it3 < 3; it3++ )
   (it4->second).residue_coordsCA[it3] = tmp_xyz[it3];
  
  tmp_xyz[0] = _sup_t[0] + _sup_u[0][0] * (it4->second).residue_coordsCB[0] + _sup_u[1][0] * (it4->second).residue_coordsCB[1] + _sup_u[2][0] * (it4->second).residue_coordsCB[2];
  tmp_xyz[1] = _sup_t[1] + _sup_u[0][1] * (it4->second).residue_coordsCB[0] + _sup_u[1][1] * (it4->second).residue_coordsCB[1] + _sup_u[2][1] * (it4->second).residue_coordsCB[2];
  tmp_xyz[2] = _sup_t[2] + _sup_u[0][2] * (it4->second).residue_coordsCB[0] + _sup_u[1][2] * (it4->second).residue_coordsCB[1] + _sup_u[2][2] * (it4->second).residue_coordsCB[2];
  
  for ( int it3 = 0; it3 < 3; it3++ )
   (it4->second).residue_coordsCB[it3] = tmp_xyz[it3];
  
  tmp_xyz[0] = _sup_t[0] + _sup_u[0][0] * (it4->second).residue_coordsSC[0] + _sup_u[1][0] * (it4->second).residue_coordsSC[1] + _sup_u[2][0] * (it4->second).residue_coordsSC[2];
  tmp_xyz[1] = _sup_t[1] + _sup_u[0][1] * (it4->second).residue_coordsSC[0] + _sup_u[1][1] * (it4->second).residue_coordsSC[1] + _sup_u[2][1] * (it4->second).residue_coordsSC[2];
  tmp_xyz[2] = _sup_t[2] + _sup_u[0][2] * (it4->second).residue_coordsSC[0] + _sup_u[1][2] * (it4->second).residue_coordsSC[1] + _sup_u[2][2] * (it4->second).residue_coordsSC[2];
  
  for ( int it3 = 0; it3 < 3; it3++ )
   (it4->second).residue_coordsSC[it3] = tmp_xyz[it3];
 }
 
 vector<CoordsProtein>::iterator it5;
 
 for ( it5 = _protein_xyz.begin(); it5 < _protein_xyz.end(); it5++ )
 {
  double tmp_xyz[3];
  
  tmp_xyz[0] = _sup_t[0] + _sup_u[0][0] * (*it5).getCoords(1) + _sup_u[1][0] * (*it5).getCoords(2) + _sup_u[2][0] * (*it5).getCoords(3);
  tmp_xyz[1] = _sup_t[1] + _sup_u[0][1] * (*it5).getCoords(1) + _sup_u[1][1] * (*it5).getCoords(2) + _sup_u[2][1] * (*it5).getCoords(3);
  tmp_xyz[2] = _sup_t[2] + _sup_u[0][2] * (*it5).getCoords(1) + _sup_u[1][2] * (*it5).getCoords(2) + _sup_u[2][2] * (*it5).getCoords(3);
  
  (*it5).setCoords( tmp_xyz[0], tmp_xyz[1], tmp_xyz[2] );
 }
 
 double tmp_xyz[3];
 
 tmp_xyz[0] = _sup_t[0] + _sup_u[0][0] * _pocket_center[0] + _sup_u[1][0] * _pocket_center[1] + _sup_u[2][0] * _pocket_center[2];
 tmp_xyz[1] = _sup_t[1] + _sup_u[0][1] * _pocket_center[0] + _sup_u[1][1] * _pocket_center[1] + _sup_u[2][1] * _pocket_center[2];
 tmp_xyz[2] = _sup_t[2] + _sup_u[0][2] * _pocket_center[0] + _sup_u[1][2] * _pocket_center[1] + _sup_u[2][2] * _pocket_center[2];
 
 for ( int it3 = 0; it3 < 3; it3++ )
  _pocket_center[it3] = tmp_xyz[it3];
 
 for ( it1 = sup_aln.begin(); it1 != sup_aln.end(); it1++ )
 {
  double tmp1_xyz[3];
  
  sup_ref->getBindingResCoords( (it1->second).res1_index, tmp1_xyz, 0 );
  
  double tmp2_xyz[3];
  
  it4 = _binding_res.find( (it1->second).res2_index );
  
  for ( int it3 = 0; it3 < 3; it3++ )
   tmp2_xyz[it3] = (*it4).second.residue_coordsCA[it3];
  
  (it1->second).sco_dst = sqrt( pow( tmp1_xyz[0]-tmp2_xyz[0], 2.0 ) + pow( tmp1_xyz[1]-tmp2_xyz[1], 2.0 ) + pow( tmp1_xyz[2]-tmp2_xyz[2], 2.0 ) );
 }
 
 return rms1;
}


// ==================================================================================   getPMbin

int Protein::getPMbin( vector<double> &pm_bin, int group1, int group2, int type1, int type2 )
{
 pm_bin.clear();
 
 for ( int it1 = 0; it1 < _binres_tot; it1++ )
 {
  map<int,binding_residue>::iterator it2;
  
  it2 = _binding_res.find( it1 );
  
  for ( int it3 = it1; it3 < _binres_tot; it3++ )
  {
   map<int,binding_residue>::iterator it4;
   
   it4 = _binding_res.find( it3 );
   
   if ( ( (*it2).second.residue_group == group1 && (*it4).second.residue_group == group2 ) || ( (*it2).second.residue_group == group2 && (*it4).second.residue_group == group1 ) )
   {
    double xyz1[3] = {0,0,0};
    double xyz2[3] = {0,0,0};
    
    if ( type1 == 0 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsCA[it5];
    else if ( type1 == 1 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsCB[it5];
    else if ( type1 == 2 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsSC[it5];
    
    if ( type2 == 0 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsCA[it5];
    else if ( type2 == 1 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsCB[it5];
    else if ( type2 == 2 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsSC[it5];
    
    double rr1 = sqrt( pow( xyz1[0] - xyz2[0], 2.0) + pow( xyz1[1] - xyz2[1], 2.0) + pow( xyz1[2] - xyz2[2], 2.0) );
    
    if ( type2 == 0 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsCA[it5];
    else if ( type2 == 1 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsCB[it5];
    else if ( type2 == 2 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz1[it5] = (*it2).second.residue_coordsSC[it5];
    
    if ( type1 == 0 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsCA[it5];
    else if ( type1 == 1 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsCB[it5];
    else if ( type1 == 2 )
     for ( int it5 = 0; it5 < 3; it5++ )
      xyz2[it5] = (*it4).second.residue_coordsSC[it5];
    
    double rr2 = sqrt( pow( xyz1[0] - xyz2[0], 2.0) + pow( xyz1[1] - xyz2[1], 2.0) + pow( xyz1[2] - xyz2[2], 2.0) );
    
    bool w1 = true;
    
    if ( (*it2).second.residue_code == "G" )
     if ( type1 == 1 )
      w1 = false;
    
    if ( (*it4).second.residue_code == "G" )
     if ( type2 == 1 )
      w1 = false;
    
    bool w2 = true;
    
    if ( (*it2).second.residue_code == "G" )
     if ( type2 == 1 )
      w2 = false;
    
    if ( (*it4).second.residue_code == "G" )
     if ( type1 == 1 )
      w2 = false;
    
    if ( it1 == it3 )
    {
     if ( type1 == type2 )
      w1 = false;
     
     w2 = false;
    }
    else
     if ( type1 == type2 )
      w2 = false;
    
    if ( w1 )
     pm_bin.push_back( rr1 );
    
    if ( w2 )
     pm_bin.push_back( rr2 );
   }
  }
 }
 
 std::sort(pm_bin.begin(), pm_bin.end(), std::less<double>());
 
 return pm_bin.size();
}


// ==================================================================================   alignmentDump

void Protein::alignmentDump( map< pair<int,int>, pair_aligned > &d1_aln, Protein * target1, std::string d1_out, double d1_score, double d1_cost, double d1_rmsd, double d1_prob, double d1_tau, double d1_tls, double d1_pms )
{
 ofstream outpdb( d1_out.c_str() );
 
 outpdb << "REMARK SCORE" << fixed << setw(7) << setprecision(3) << d1_score << endl
        << "REMARK NRES " << fixed << setw(7) << d1_aln.size() << endl
        << "REMARK COST " << fixed << setw(7) << setprecision(3) << d1_cost << endl
        << "REMARK RMSD " << fixed << setw(7) << setprecision(3) << d1_rmsd << endl
        << "REMARK PROB " << fixed << setw(7) << setprecision(3) << d1_prob << endl
        << "REMARK TAU  " << fixed << setw(7) << setprecision(3) << d1_tau << endl
        << "REMARK TLS  " << fixed << setw(7) << setprecision(3) << d1_tls << endl
        << "REMARK PMS  " << fixed << setw(7) << setprecision(3) << d1_pms << endl;
 
 map< pair<int,int>, pair_aligned >::iterator it1;
 
 for ( it1 = d1_aln.begin(); it1 != d1_aln.end(); it1++ )
 {
  map<int,binding_residue>::iterator it2;
  
  it2 = _binding_res.find( (it1->second).res2_index );
  
  outpdb << "REMARK ALN" << fixed << setw(5) << setprecision(0) << (it1->first).first
                         << fixed << setw(2) << setprecision(0) << target1->getBindingResCode( (it1->second).res1_index )
                         << fixed << setw(5) << setprecision(0) << (it1->first).second
                         << fixed << setw(2) << setprecision(0) << (*it2).second.residue_code
                         << fixed << setw(8) << setprecision(3) << (it1->second).sco_dst
                         << fixed << setw(8) << setprecision(3) << (it1->second).sco_svr
                         << fixed << setw(8) << setprecision(4) << (it1->second).sco_svc << endl;
 }
 
 outpdb << "REMARK MTX_T" << fixed << setw(9) << setprecision(3) << _sup_t[0] << setw(9) << _sup_t[1] << setw(9) << _sup_t[2] << endl
        << "REMARK MTX_U" << fixed << setw(9) << setprecision(3) << _sup_u[0][0] << setw(9) << _sup_u[0][1] << setw(9) << _sup_u[0][2] << endl
        << "REMARK MTX_U" << fixed << setw(9) << setprecision(3) << _sup_u[1][0] << setw(9) << _sup_u[1][1] << setw(9) << _sup_u[1][2] << endl
        << "REMARK MTX_U" << fixed << setw(9) << setprecision(3) << _sup_u[2][0] << setw(9) << _sup_u[2][1] << setw(9) << _sup_u[2][2] << endl;
 
 vector<CoordsProtein>::iterator it3;
 
 for ( it3 = _protein_xyz.begin(); it3 < _protein_xyz.end(); it3++ )
 {
  double tmp_xyz[3];
  
  for ( int it4 = 0; it4 < 3; it4++ )
   tmp_xyz[it4] = (*it3).getCoords(it4+1);
  
  outpdb << "ATOM  " << setw(5) << (*it3).getAtomNumber()
                     << setw(5) << (*it3).getAtomName()
                     << setw(4) << (*it3).getResidueName()
                     << setw(2) << (*it3).getChainID()
                     << setw(4) << (*it3).getResidueNumber()
                     << fixed << setw(12) << setprecision(3) << tmp_xyz[0]
                     << fixed << setw(8)  << setprecision(3) << tmp_xyz[1]
                     << fixed << setw(8)  << setprecision(3) << tmp_xyz[2] << endl;
 }
 
 outpdb << "END" << endl;
 
 outpdb.close();
}
