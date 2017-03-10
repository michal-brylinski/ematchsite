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


#ifndef __PROTEIN_H_
#define __PROTEIN_H_

#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<bitset>
#include<map>
#include<list>
#include<vector>
#include<utility>
#include<algorithm>
#include<cmath>
#include<gzstream.h>

#include "size.h"
#include "coords.h"
#include "data.h"
#include "tanimoto.h"
#include "rmsd.h"

using namespace std;

struct binding_residue
{
 int         residue_number;
 std::string residue_code;
 int         residue_group;
 bool        residue_binding;
 double      residue_score;
 double      residue_distance;
 double      residue_fraction;
 double      residue_coordsCA[3];
 double      residue_coordsCB[3];
 double      residue_coordsSC[3];
};

struct ligand_fpt_smi
{
 int                 ftp_number;
 double              fpt_fraction;
 std::string         fpt_pdbid;
 std::bitset<MAXSMI> fpt_fingerprint;
};

struct ligand_fpt_mac
{
 int                 ftp_number;
 double              fpt_fraction;
 std::string         fpt_pdbid;
 std::bitset<MAXMAC> fpt_fingerprint;
};

struct ligand_bond
{
 int atom1;
 int atom2;
 int bond_type;
};

struct compound_match
{
 double                  tanimoto;
 int                     total_pairs;
 vector< pair<int,int> > atom_pairs;
};

struct tpl_alignment
{
 int                        tpl_res;
 int                        tpl_len;
 double                     tpl_tms;
 double                     tpl_rms;
 double                     tpl_sid;
 vector< pair<int,string> > tpl_pairs;
};

struct pair_aligned
{
 int res1_index;
 int res2_index;
 
 double sco_svc;
 double sco_svr;
 double sco_dst;
};

class Protein {
        
  private:
    
    int                        _protein_na;                   // number of protein atoms
    int                        _protein_nr;                   // number of protein residues
    
    vector<CoordsProtein>      _protein_xyz;                  // protein coords
    
    string                     _protein_seq1;                 // aa sequence
    char                       _protein_seq2[MAXPRO];         // aa sequence
    int                        _protein_seq3[MAXPRO];         // aa sequence numbering
    double                     _protein_prf1[MAXPRO][20];     // sequence profiles
    double                     _protein_prf2[MAXPRO][3];      // secondary structure profiles
    double                     _protein_prf3[MAXPRO][20];     // hydrophobicity profiles
    map<string,tpl_alignment>  _protein_aln1;                 // template structure alignments
    
    int                        _pocket_number;                // pocket number
    int                        _ligclust_tot_smi;             // number of ligand clusters SMILES
    int                        _ligclust_tot_mac;             // number of ligand clusters MACCS
    
    int                        _protein_ntpl;                 // number of templates
    
    double                     _pocket_center[4];             // ligand geometric center
    
    int                        _binres_tot;                   // number of predicted binding residues
    
    map<int,binding_residue>   _binding_res;                  // ligand-binding residues
    
    list<ligand_fpt_smi>       _ligand_fpts_smi;              // fingerprints Daylight
    list<ligand_fpt_mac>       _ligand_fpts_mac;              // fingerprints MACCS
    
    double                     _profile_fpts_smi[MAXSMI];     // fingerprint profile Daylight
    double                     _profile_fpts_mac[MAXMAC];     // fingerprint profile MACCS
    
    double                     _ave_mw;                       // weighted average mw
    double                     _ave_logp;                     // weighted average logp
    double                     _ave_psa;                      // weighted average psa
    double                     _ave_mr;                       // weighted average mr
    double                     _ave_hbd;                      // weighted average hb donors
    double                     _ave_hba;                      // weighted average hb acceptors
    
    double                     _std_mw;                       // weighted stdev mw
    double                     _std_logp;                     // weighted stdev logp
    double                     _std_psa;                      // weighted stdev psa
    double                     _std_mr;                       // weighted stdev mr
    double                     _std_hbd;                      // weighted stdev hb donors
    double                     _std_hba;                      // weighted stdev hb acceptors
    
    double                     _pocket_plb;                   // protein-ligand binding index
    
    double                     _svm_pocket_confidence;        // pocket confidence
    
    int                        _cmps_nc;                      // number of template-bound compounds
    
    int                        _cmps_nla[MAXCMP];             // number of ligand heavy atoms
    int                        _cmps_nlb[MAXCMP];             // number of ligand bonds
    
    string                     _cmps_id[MAXCMP];              // ligand PDB-IDs
    
    vector<CoordsLigand>       _cmps_xyz[MAXCMP];             // ligand heavy atom coords
    
    vector<ligand_bond>        _cmps_bnd[MAXCMP];             // ligand bonds
    
    double                     _cmps_cen[MAXCMP][3];          // ligand geometric centers
    
    bitset<MAXSMI>             _cmps_fpt_smi[MAXCMP];         // SMILES fingerprints
    bitset<MAXMAC>             _cmps_fpt_mac[MAXCMP];         // MACCS fingerprints
    
    double                     _cmps_mw[MAXCMP];              // molecular weight
    double                     _cmps_logp[MAXCMP];            // logp
    double                     _cmps_psa[MAXCMP];             // polar surface area
    double                     _cmps_mr[MAXCMP];              // molar refractivity
    
    int                        _cmps_hbd[MAXCMP];             // HB donors
    int                        _cmps_hba[MAXCMP];             // HB acceptors
    
    map<string,compound_match> _cmps_matching[MAXCMP];        // compound matching
    
    double                     _sup_t[3];                     // translation vector
    double                     _sup_u[3][3];                  // rotation matrix
    
  public:
    
    Protein( void );
    
    ~Protein();
    
    bool loadStructure( std::string );
    
    bool loadPsipred( std::string );
    
    bool loadSequence( std::string );
    
    bool loadPocket( std::string, int );
    
    bool loadLigands( std::string, int );
    
    bool loadAlignments( std::string );
    
    bool loadLibrary( std::string, int, std::string [], int & );
    
    double getTanimotoSMILES( std::bitset<MAXSMI>, std::string );
    
    double getTanimotoMACCS( std::bitset<MAXMAC>, std::string );
    
    double getPocketProperty( int );
    
    int matchLigands( std::string, int, std::string [], bool );
    
    int getCmpsTotal( void );
    
    int getCmpsList( std::string [] );
    
    std::string getLigandName( int );
    
    int getLigandNum( std::string );
    
    int getLigandCoords( int, double [][3] );
    
    int getLigandAtomNames( int, string [] );
    
    int getLigandBonds( int, int [][3] );
    
    int getBindingResiduesTotal( void );
    
    void getBindingResCoords( int, double [], int );
    
    void getResidueCoords( int, vector<std::string> & );
    
    int getBindingResNum( int );
    
    std::string getBindingResCode( int );
    
    int getBindingResGroup( int );
    
    void getSeqProfile( int, double [] );
    
    void getSecProfile( int, double [] );
    
    void getHPProfile( int, double [] );
    
    double getBindingProb( int );
    
    void getNeighborDist( int, vector<double> & );
    
    double getResEntropy( int );
    
    int getLigandMatch( int, std::string, std::pair<int,int> [], double & );
    
    double alignmentSuperpose( map< pair<int,int>, pair_aligned > &, Protein * );
    
    int getPMbin( vector<double> &, int, int, int, int );
    
    void alignmentDump( map< pair<int,int>, pair_aligned > &, Protein *, std::string, double, double, double, double, double, double, double );
};

#endif
