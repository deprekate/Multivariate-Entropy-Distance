#ifndef TYPEDEFBASE_H
#define TYPEDEFBASE_H

#include <string>
#include<vector>
#include<utility>
typedef std::vector<int> Ve_I;
typedef std::vector<std::string> Ve_Str;
typedef const std::vector<int> con_Ve_I;
typedef std::vector< std::pair<int,int> > Ve_Pa_I_I;
typedef std::vector<double> Ve_D;
typedef const std::vector<double> con_Ve_D;
typedef std::vector< std::vector<int> > Ve_Ve_I;
typedef std::vector< std::vector<double> > Ve_Ve_D;
typedef std::vector< double* > Ve_D_Ptr;
typedef std::vector< std::vector<double>* > Ve_Ve_D_Ptr;
typedef const std::vector< std::vector<double> > con_Ve_Ve_D;
typedef std::vector< std::pair<std::string, int> > Ve_Pa_Str_I;
typedef std::pair<int, int> Pa_I_I;
typedef std::pair<int, double> Pa_I_D;
typedef std::pair<bool, double> Pa_B_D;
typedef std::pair<bool, int> Pa_B_I;
typedef std::pair<double, double> Pa_D_D;
typedef std::pair<double, int> Pa_D_I;
typedef const std::pair<int, int> con_Pa_I_I;
typedef std::vector< std::pair< double, double > > Ve_Pa_D_D;
typedef std::vector< std::pair< double, int > > Ve_Pa_D_I;
typedef std::pair<std::vector<double>,std::vector<int> > Pa_Ve_D_Ve_I;
#include<fstream>
typedef std::ifstream ifstream;
typedef std::ofstream ofstream;
typedef std::ostream ostream;
typedef std::istream istream;

//#include<string>
typedef std::string Str;
typedef const std::string con_Str;

#include<deque>
typedef const std::deque<int> con_De_I; 
typedef std::deque<int> De_I; 

#include<map>
typedef std::map< int, int > Map_I_I;
typedef std::map< int, std::string > Map_I_Str;

#include<set>
typedef std::set<int> Set_I;
typedef std::set<Str> Set_Str;
typedef std::set<double> Set_D;

#include<limits>
typedef std::numeric_limits<int> int_limits;

#include<list>
typedef std::list< std::pair< double, double > > List_Pa_D_D;
typedef std::list< int > List_I;

#include<iostream>
#include<iomanip>
#include<list>
#include"assert.h"
#include<algorithm>
#endif //TYPEDEFBASE_H


