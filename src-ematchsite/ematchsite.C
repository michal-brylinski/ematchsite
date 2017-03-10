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


#include "ematchsite.h"

using namespace std;

int main(int argc, char *argv[])
{
 time_t t_start, t_end, t_bench1, t_bench2;
 time(&t_start);
 
 cout << "------------------------------------------------------------" << endl
      << "                         ematchsite" << endl
      << "                        version 1.1" << endl
      << "   sequence order independent alignment of binding sites" << endl << endl
      << "       report bugs and issues to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 5 )
 {
  cout << " ematchsite -i <input arguments>" << endl
       << "            -k <kcombu alignments>" << endl
       << "            -o <output filename>" << endl << endl
       << " input file options:" << endl << endl
       << "            -structureA, -structureB   <target structure (A/B)>" << endl
       << "            -profilesA, -profilesB     <sequence profile (A/B)>" << endl
       << "            -secstrA, -secstrB         <secondary structure profiles (A/B)>" << endl
       << "            -pocketsA, -pocketsB       <eFindSite pockets (A/B)>" << endl
       << "            -numberA, -numberB         <eFindSite pocket number (A/B, default 1)>" << endl
       << "            -alignmentsA, -alignmentsB <eFindSite alignments (A/B)>" << endl
       << "            -ligandsA, -ligandsB       <eFindSite ligands (A/B)>" << endl << endl
       << " virtual screening options:" << endl << endl
       << "            -m <scoring function (default sum)>" << endl << endl
       << "                tst - classical Tanimoto coeff for Daylight" << endl
       << "                tsa - average Tanimoto coeff for Daylight" << endl
       << "                tsc - continuous Tanimoto coeff for Daylight" << endl
       << "                tmt - classical Tanimoto coeff for MACCS" << endl
       << "                tma - average Tanimoto coeff for MACCS" << endl
       << "                tmc - continuous Tanimoto coeff for MACCS" << endl
       << "                sum - data fusion using sum rule (default)" << endl
       << "                max - data fusion using max rule" << endl
       << "                min - data fusion using min rule" << endl
       << "                svm - machine learning using SVM" << endl << endl;
  exit(EXIT_SUCCESS);
 }
 
 string args_name;
 bool args_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 string kcombu_name;
 bool kcombu_opt = false;
 
 string structureA_name;
 bool structureA_opt = false;
 
 string structureB_name;
 bool structureB_opt = false;
 
 string profilesA_name;
 bool profilesA_opt = false;
 
 string profilesB_name;
 bool profilesB_opt = false;
 
 string secstrA_name;
 bool secstrA_opt = false;
 
 string secstrB_name;
 bool secstrB_opt = false;
 
 string pocketsA_name;
 bool pocketsA_opt = false;
 
 string pocketsB_name;
 bool pocketsB_opt = false;
 
 string alignmentsA_name;
 bool alignmentsA_opt = false;
 
 string alignmentsB_name;
 bool alignmentsB_opt = false;
 
 string ligandsA_name;
 bool ligandsA_opt = false;
 
 string ligandsB_name;
 bool ligandsB_opt = false;
 
 int numberA = 1;
 bool numberA_opt = false;
 
 int numberB = 1;
 bool numberB_opt = false;
 
 std::string smethod = "sum";
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-i")           && i < argc ) { args_name          = string(argv[i+1]); args_opt        = true; }
  if ( !strcmp(argv[i],"-k")           && i < argc ) { kcombu_name        = string(argv[i+1]); kcombu_opt      = true; }
  if ( !strcmp(argv[i],"-o")           && i < argc ) { output_name        = string(argv[i+1]); output_opt      = true; }
  if ( !strcmp(argv[i],"-structureA")  && i < argc ) { structureA_name    = string(argv[i+1]); structureA_opt  = true; }
  if ( !strcmp(argv[i],"-structureB")  && i < argc ) { structureB_name    = string(argv[i+1]); structureB_opt  = true; }
  if ( !strcmp(argv[i],"-profilesA")   && i < argc ) { profilesA_name     = string(argv[i+1]); profilesA_opt   = true; }
  if ( !strcmp(argv[i],"-profilesB")   && i < argc ) { profilesB_name     = string(argv[i+1]); profilesB_opt   = true; }
  if ( !strcmp(argv[i],"-secstrA")     && i < argc ) { secstrA_name       = string(argv[i+1]); secstrA_opt     = true; }
  if ( !strcmp(argv[i],"-secstrB")     && i < argc ) { secstrB_name       = string(argv[i+1]); secstrB_opt     = true; }
  if ( !strcmp(argv[i],"-pocketsA")    && i < argc ) { pocketsA_name      = string(argv[i+1]); pocketsA_opt    = true; }
  if ( !strcmp(argv[i],"-pocketsB")    && i < argc ) { pocketsB_name      = string(argv[i+1]); pocketsB_opt    = true; }
  if ( !strcmp(argv[i],"-alignmentsA") && i < argc ) { alignmentsA_name   = string(argv[i+1]); alignmentsA_opt = true; }
  if ( !strcmp(argv[i],"-alignmentsB") && i < argc ) { alignmentsB_name   = string(argv[i+1]); alignmentsB_opt = true; }
  if ( !strcmp(argv[i],"-ligandsA")    && i < argc ) { ligandsA_name      = string(argv[i+1]); ligandsA_opt    = true; }
  if ( !strcmp(argv[i],"-ligandsB")    && i < argc ) { ligandsB_name      = string(argv[i+1]); ligandsB_opt    = true; }
  if ( !strcmp(argv[i],"-numberA")     && i < argc ) { numberA            = atoi(argv[i+1]);   numberA_opt     = true; }
  if ( !strcmp(argv[i],"-numberB")     && i < argc ) { numberB            = atoi(argv[i+1]);   numberB_opt     = true; }
  if ( !strcmp(argv[i],"-m")           && i < argc ) { smethod            = string(argv[i+1]);                         }
 }
 
 char * path1;
 
 path1 = getenv("EM_LIB"); if ( path1==NULL ) { cout << "EM_LIB is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EM_MOD"); if ( path1==NULL ) { cout << "EM_MOD is not set" << endl; exit(EXIT_FAILURE); }
 
 string lib_path;
 lib_path = getenv("EM_LIB");
 
 string model_path;
 model_path = getenv("EM_MOD");
 
 ifstream f01( (model_path+"/matchSVC.model").c_str() );
 ifstream f02( (model_path+"/matchSVR.model").c_str() );
 ifstream f03( (model_path+"/matchSVM.scale").c_str() );
 ifstream f04( (model_path+"/scoringSVM.model").c_str() );
 ifstream f05( (model_path+"/scoringSVM.scale").c_str() );
 ifstream f06( (model_path+"/probSVC.model").c_str() );
 ifstream f07( (model_path+"/probSVC.scale").c_str() );
 
 ifstream l01( (lib_path+"/escreen.gz").c_str() );
 ifstream l02( (lib_path+"/cmps.lst").c_str() );
 
 if ( !f01 || !f02 || !f03 || !f04 || !f05 || !f06 || !f07 )
 {
  cout << "Could not find SVM models in " << model_path << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !l01 || ( !kcombu_opt && !l02 ) )
 {
  cout << "Incomplete eMatchSite library at " << lib_path << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( args_opt )
 {
       if ( structureA_opt )
  {
   cout << "Argument -structureA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( structureB_opt )
  {
   cout << "Argument -structureB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( profilesA_opt )
  {
   cout << "Argument -profilesA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( profilesB_opt )
  {
   cout << "Argument -profilesB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( secstrA_opt )
  {
   cout << "Argument -secstrA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( secstrB_opt )
  {
   cout << "Argument -secstrB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( pocketsA_opt )
  {
   cout << "Argument -pocketsA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( pocketsB_opt )
  {
   cout << "Argument -pocketsB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( alignmentsA_opt )
  {
   cout << "Argument -alignmentsA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( alignmentsB_opt )
  {
   cout << "Argument -alignmentsB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( ligandsA_opt )
  {
   cout << "Argument -ligandsA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( ligandsB_opt )
  {
   cout << "Argument -ligandsB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( numberA_opt )
  {
   cout << "Argument -numberA is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
  else if ( numberB_opt )
  {
   cout << "Argument -numberB is not allowed when using input arguments file (" << args_name << ")" << endl;
   exit(EXIT_FAILURE);
  }
 }
 
 if ( args_opt )
 {
  string line1;
  
  ifstream args_file( args_name.c_str() );
  
  if ( !args_file.is_open() )  { cout << "Cannot open " << args_name << endl; exit(EXIT_FAILURE); }
  
  while (getline(args_file,line1))
  {
   if ( line1.size() > 10 )
   {
    if ( line1.substr(0,10) == "structureA")
    {
     structureA_name = line1.substr(10,line1.length()-10);
     
     structureA_name.erase(std::remove(structureA_name.begin(), structureA_name.end(), ' '), structureA_name.end());
     structureA_name.erase(std::remove(structureA_name.begin(), structureA_name.end(), '\t'), structureA_name.end());
     
     ifstream inp1( (structureA_name).c_str() );
     
     if ( inp1 )
      structureA_opt = true;
    }
    else if ( line1.substr(0,10) == "structureB")
    {
     structureB_name = line1.substr(10,line1.length()-10);
     
     structureB_name.erase(std::remove(structureB_name.begin(), structureB_name.end(), ' '), structureB_name.end());
     structureB_name.erase(std::remove(structureB_name.begin(), structureB_name.end(), '\t'), structureB_name.end());
     
     ifstream inp1( (structureB_name).c_str() );
     
     if ( inp1 )
      structureB_opt = true;
    }
   }
   
   if ( line1.size() > 9 )
   {
    if ( line1.substr(0,9) == "profilesA")
    {
     profilesA_name = line1.substr(9,line1.length()-9);
     
     profilesA_name.erase(std::remove(profilesA_name.begin(), profilesA_name.end(), ' '), profilesA_name.end());
     profilesA_name.erase(std::remove(profilesA_name.begin(), profilesA_name.end(), '\t'), profilesA_name.end());
     
     ifstream inp1( (profilesA_name).c_str() );
     
     if ( inp1 )
      profilesA_opt = true;
    }
    else if ( line1.substr(0,9) == "profilesB")
    {
     profilesB_name = line1.substr(9,line1.length()-9);
     
     profilesB_name.erase(std::remove(profilesB_name.begin(), profilesB_name.end(), ' '), profilesB_name.end());
     profilesB_name.erase(std::remove(profilesB_name.begin(), profilesB_name.end(), '\t'), profilesB_name.end());
     
     ifstream inp1( (profilesB_name).c_str() );
     
     if ( inp1 )
      profilesB_opt = true;
    }
   }
   
   if ( line1.size() > 7 )
   {
    if ( line1.substr(0,7) == "secstrA")
    {
     secstrA_name = line1.substr(7,line1.length()-7);
     
     secstrA_name.erase(std::remove(secstrA_name.begin(), secstrA_name.end(), ' '), secstrA_name.end());
     secstrA_name.erase(std::remove(secstrA_name.begin(), secstrA_name.end(), '\t'), secstrA_name.end());
     
     ifstream inp1( (secstrA_name).c_str() );
     
     if ( inp1 )
      secstrA_opt = true;
    }
    else if ( line1.substr(0,7) == "secstrB")
    {
     secstrB_name = line1.substr(7,line1.length()-7);
     
     secstrB_name.erase(std::remove(secstrB_name.begin(), secstrB_name.end(), ' '), secstrB_name.end());
     secstrB_name.erase(std::remove(secstrB_name.begin(), secstrB_name.end(), '\t'), secstrB_name.end());
     
     ifstream inp1( (secstrB_name).c_str() );
     
     if ( inp1 )
      secstrB_opt = true;
    }
   }
   
   if ( line1.size() > 8 )
   {
    if ( line1.substr(0,8) == "pocketsA")
    {
     pocketsA_name = line1.substr(8,line1.length()-8);
     
     pocketsA_name.erase(std::remove(pocketsA_name.begin(), pocketsA_name.end(), ' '), pocketsA_name.end());
     pocketsA_name.erase(std::remove(pocketsA_name.begin(), pocketsA_name.end(), '\t'), pocketsA_name.end());
     
     ifstream inp1( (pocketsA_name).c_str() );
     
     if ( inp1 )
      pocketsA_opt = true;
    }
    else if ( line1.substr(0,8) == "pocketsB")
    {
     pocketsB_name = line1.substr(8,line1.length()-8);
     
     pocketsB_name.erase(std::remove(pocketsB_name.begin(), pocketsB_name.end(), ' '), pocketsB_name.end());
     pocketsB_name.erase(std::remove(pocketsB_name.begin(), pocketsB_name.end(), '\t'), pocketsB_name.end());
     
     ifstream inp1( (pocketsB_name).c_str() );
     
     if ( inp1 )
      pocketsB_opt = true;
    }
   }
   
   if ( line1.size() > 11 )
   {
    if ( line1.substr(0,11) == "alignmentsA")
    {
     alignmentsA_name = line1.substr(11,line1.length()-11);
     
     alignmentsA_name.erase(std::remove(alignmentsA_name.begin(), alignmentsA_name.end(), ' '), alignmentsA_name.end());
     alignmentsA_name.erase(std::remove(alignmentsA_name.begin(), alignmentsA_name.end(), '\t'), alignmentsA_name.end());
     
     ifstream inp1( (alignmentsA_name).c_str() );
     
     if ( inp1 )
      alignmentsA_opt = true;
    }
    else if ( line1.substr(0,11) == "alignmentsB")
    {
     alignmentsB_name = line1.substr(11,line1.length()-11);
     
     alignmentsB_name.erase(std::remove(alignmentsB_name.begin(), alignmentsB_name.end(), ' '), alignmentsB_name.end());
     alignmentsB_name.erase(std::remove(alignmentsB_name.begin(), alignmentsB_name.end(), '\t'), alignmentsB_name.end());
     
     ifstream inp1( (alignmentsB_name).c_str() );
     
     if ( inp1 )
      alignmentsB_opt = true;
    }
   }
   
   if ( line1.size() > 8 )
   {
    if ( line1.substr(0,8) == "ligandsA")
    {
     ligandsA_name = line1.substr(8,line1.length()-8);
     
     ligandsA_name.erase(std::remove(ligandsA_name.begin(), ligandsA_name.end(), ' '), ligandsA_name.end());
     ligandsA_name.erase(std::remove(ligandsA_name.begin(), ligandsA_name.end(), '\t'), ligandsA_name.end());
     
     ifstream inp1( (ligandsA_name).c_str() );
     
     if ( inp1 )
      ligandsA_opt = true;
    }
    else if ( line1.substr(0,8) == "ligandsB")
    {
     ligandsB_name = line1.substr(8,line1.length()-8);
     
     ligandsB_name.erase(std::remove(ligandsB_name.begin(), ligandsB_name.end(), ' '), ligandsB_name.end());
     ligandsB_name.erase(std::remove(ligandsB_name.begin(), ligandsB_name.end(), '\t'), ligandsB_name.end());
     
     ifstream inp1( (ligandsB_name).c_str() );
     
     if ( inp1 )
      ligandsB_opt = true;
    }
   }
   
   if ( line1.size() > 7 )
   {
    if ( line1.substr(0,7) == "numberA")
    {
     string number_tmp = line1.substr(7,line1.length()-7);
     
     number_tmp.erase(std::remove(number_tmp.begin(), number_tmp.end(), ' '), number_tmp.end());
     number_tmp.erase(std::remove(number_tmp.begin(), number_tmp.end(), '\t'), number_tmp.end());
     
     numberA = atoi(number_tmp.c_str());
     
     numberA_opt = true;
    }
    else if ( line1.substr(0,7) == "numberB")
    {
     string number_tmp = line1.substr(7,line1.length()-7);
     
     number_tmp.erase(std::remove(number_tmp.begin(), number_tmp.end(), ' '), number_tmp.end());
     number_tmp.erase(std::remove(number_tmp.begin(), number_tmp.end(), '\t'), number_tmp.end());
     
     numberB = atoi(number_tmp.c_str());
     
     numberB_opt = true;
    }
   }
  }
  
  args_file.close();
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !structureA_opt )
 {
  cout << "Provide target structure (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !structureB_opt )
 {
  cout << "Provide target structure (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !profilesA_opt )
 {
  cout << "Provide sequence profiles (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !profilesB_opt )
 {
  cout << "Provide sequence profiles (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !secstrA_opt )
 {
  cout << "Provide secondary structure profiles (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !secstrB_opt )
 {
  cout << "Provide secondary structure profiles (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !pocketsA_opt )
 {
  cout << "Provide pockets predicted by eFindSite (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !pocketsB_opt )
 {
  cout << "Provide pockets predicted by eFindSite (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !alignmentsA_opt )
 {
  cout << "Provide alignments constructed by eFindSite (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !alignmentsB_opt )
 {
  cout << "Provide alignments constructed by eFindSite (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !ligandsA_opt )
 {
  cout << "Provide ligands identified by eFindSite (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !ligandsB_opt )
 {
  cout << "Provide ligands identified by eFindSite (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( numberA < 1 )
 {
  cout << "Argument numberA must be a positive integer" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( numberB < 1 )
 {
  cout << "Argument numberB must be a positive integer" << endl;
  exit(EXIT_FAILURE);
 }
 
 ModelSVM * model_svm;
 
 model_svm = new ModelSVM( false, false, false, false );
 
 cout << "Loading SVM models " << flush;
 
 time(&t_bench1);
 
 model_svm->loadModel( 1, model_path+"/matchSVC.model" ); cout << "." << flush;
 model_svm->loadModel( 2, model_path+"/matchSVR.model" ); cout << "." << flush;
 model_svm->loadModel( 3, model_path+"/scoringSVM.model" ); cout << "." << flush;
 model_svm->loadModel( 4, model_path+"/probSVC.model" ); cout << "." << flush;
 
 model_svm->loadScale( 1, model_path+"/matchSVM.scale" ); cout << "." << flush;
 model_svm->loadScale( 2, model_path+"/matchSVM.scale" ); cout << "." << flush;
 model_svm->loadScale( 3, model_path+"/scoringSVM.scale" ); cout << "." << flush;
 model_svm->loadScale( 4, model_path+"/probSVC.scale" ); cout << "." << flush;
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 Protein * proteinA;
 
 proteinA = new Protein();
 
 Protein * proteinB;
 
 proteinB = new Protein();
 
 Screen * cmp_library;
 
 cmp_library = new Screen();
 
 cout << "Loading screening library ..." << flush;
 
 time(&t_bench1);
 
 if ( cmp_library->loadLibrary( (lib_path+"/escreen.gz").c_str() ) )
 {
  cout << "Cannot load eMatchSite library: " << lib_path << "/escreen.gz" << endl;
  exit(EXIT_FAILURE);
 }
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Reading input data ." << flush;
 
 time(&t_bench1);
 
 if ( proteinA->loadStructure(structureA_name) )
 {
  cout << "Cannot read target structure (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadStructure(structureB_name) )
 {
  cout << "Cannot read target structure (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinA->loadPsipred(secstrA_name) )
 {
  cout << "Cannot read secondary structure profiles (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadPsipred(secstrB_name) )
 {
  cout << "Cannot read secondary structure profiles (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinA->loadSequence(profilesA_name) )
 {
  cout << "Cannot read sequence profiles (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadSequence(profilesB_name) )
 {
  cout << "Cannot read sequence profiles (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinA->loadPocket(pocketsA_name, numberA) )
 {
  cout << "Cannot read eFindSite pockets (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadPocket(pocketsB_name, numberB) )
 {
  cout << "Cannot read eFindSite pockets (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinA->loadLigands(ligandsA_name, numberA) )
 {
  cout << "Cannot read eFindSite ligands (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadLigands(ligandsB_name, numberB) )
 {
  cout << "Cannot read eFindSite ligands (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinA->loadAlignments(alignmentsA_name) )
 {
  cout << "Cannot read eFindSite alignments (A)" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "." << flush;
 
 if ( proteinB->loadAlignments(alignmentsB_name) )
 {
  cout << "Cannot read eFindSite alignments (B)" << endl;
  exit(EXIT_FAILURE);
 }
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 string cmps_idA[MAXCMP];
 string cmps_idB[MAXCMP];
 
 int cmps_ncA = proteinA->getCmpsList( cmps_idA );
 int cmps_ncB = proteinB->getCmpsList( cmps_idB );
 
 int cmps_npA = 0;
 int cmps_npB = 0;
 
 time(&t_bench1);
 
 if ( kcombu_opt )
 {
  cout << "Reading kcombu alignments ." << flush;
  
  cmps_npA = proteinA->matchLigands( kcombu_name, cmps_ncB, cmps_idB, true );
  
  cout << "." << flush;
  
  cmps_npB = proteinB->matchLigands( kcombu_name, cmps_ncA, cmps_idA, false );
  
  cout << "." << flush;
 }
 else
 {
  cout << "Loading eMatchSite library ." << flush;
  
  if ( proteinA->loadLibrary( lib_path, cmps_ncB, cmps_idB, cmps_npA ) )
  {
   cout << "Cannot load eMatchSite library at " << lib_path << endl;
   exit(EXIT_FAILURE);
  }
  
  cout << "." << flush;
  
  if ( proteinB->loadLibrary( lib_path, cmps_ncA, cmps_idA, cmps_npB ) )
  {
   cout << "Cannot load eMatchSite library at " << lib_path << endl;
   exit(EXIT_FAILURE);
  }
 
  cout << "." << flush;
 }
 
 time(&t_bench2);
 
 if ( cmps_npA > 0 && cmps_npB > 0 )
 {
  if ( cmps_npA > 1 )
   cout << " " << fixed << setprecision(3) << cmps_npA << setprecision(0) << " pairs (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
  else
   cout << " " << fixed << setprecision(3) << cmps_npA << setprecision(0) << " pair (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 }
 else
 {
  cout << " no similar ligands found (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
  
  delete proteinA;
  delete proteinB;
  
  delete cmp_library;
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_FAILURE);
 }
 
 cout << "Constructing binding site alignment ." << flush;
 
 time(&t_bench1);
 
 map< pair<int,int>, pair_score > site_pairwise;
 
 SitePairMake( site_pairwise, proteinA, proteinB );
 
 cout << "." << flush;
 
 SitePairSVM( site_pairwise, model_svm );
 
 cout << "." << flush;
 
 map< pair<int,int>, pair_aligned > site_alignment;
 
 double score_mat = SitePairMunkres( site_pairwise, site_alignment, proteinA->getBindingResiduesTotal(), proteinB->getBindingResiduesTotal() );
 
 if ( site_alignment.size() < 3 )
 {
  cout << " best solution is too short (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
  
  delete proteinA;
  delete proteinB;
  
  delete cmp_library;
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_FAILURE);
 }
 
 time(&t_bench2);
 
 cout << " " << site_alignment.size() << " residues (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Superposing pockets ..." << flush;
 
 time(&t_bench1);
 
 double rmsd_mat = proteinB->alignmentSuperpose( site_alignment, proteinA );
 
 time(&t_bench2);
 
 cout << " " << fixed << setprecision(3) << rmsd_mat << setprecision(0) << " A (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Assigning probability score ..." << flush;
 
 time(&t_bench1);
 
 double prob_mat = SitePairScoreProb( site_alignment );
 
 time(&t_bench2);
 
 cout << " " << fixed << setprecision(3) << prob_mat << setprecision(0) << " (" << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Calculating Kendall Tau ." << flush;
 
 time(&t_bench1);
 
 cmp_library->rankLibrary(model_svm, proteinA, 0, smethod);
 
 cout << "." << flush;
 
 cmp_library->rankLibrary(model_svm, proteinB, 1, smethod);
 
 cout << "." << flush;
 
 cmp_library->calculateKendallTau();
 
 time(&t_bench2);
 
 cout << " " << fixed << setprecision(3) << cmp_library->getKendallTau() << " (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cmp_library->cleanCompounds();
 
 cout << "Calculating TLscore ..." << flush;
 
 time(&t_bench1);
 
 double tls_mat = SitePairTLscore( proteinA, proteinB );
 
 if ( tls_mat < 0.001 )
  tls_mat = 0.001;
 
 time(&t_bench2);
 
 cout << " " << fixed << setprecision(3) << tls_mat << " (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Calculating PMscore ..." << flush;
 
 time(&t_bench1);
 
 double pms_mat = SitePairPMscore( proteinA, proteinB );
 
 time(&t_bench2);
 
 cout << " " << fixed << setprecision(3) << pms_mat << " (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Final scoring ..." << flush;
 
 double ems_svm[MAXSV4] = { score_mat, rmsd_mat, prob_mat, cmp_library->getKendallTau(), log(tls_mat), pms_mat };
 
 double ems_sco = model_svm->SVMpredict( 4, ems_svm );
 
 cout << " " << fixed << setprecision(3) << ems_sco << " (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Writing output ..." << flush;
 
 time(&t_bench1);
 
 proteinB->alignmentDump( site_alignment, proteinA, output_name, ems_sco, score_mat, rmsd_mat, prob_mat, cmp_library->getKendallTau(), tls_mat, pms_mat );
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 delete proteinA;
 delete proteinB;
 
 delete cmp_library;
 
 time(&t_end);
 
 printTime( difftime(t_end, t_start) );
 
 return 0;
}
