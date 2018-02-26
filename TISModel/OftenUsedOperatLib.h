#ifndef OFTENUSEDOPERATLIB_H
#define OFTENUSEDOPERATLIB_H

#include"TypeDef.h"

double getDistance	
( con_Ve_D& rhs, con_Ve_D& lhs );

int getMin( Ve_D lhs, con_Ve_D& coefficent = Ve_D( ) );

void addToPreVector( Ve_D& lhs, Ve_D& rhs, Ve_I counter = Ve_I( ) );

double GetDerivate( Ve_D& v, bool l = false );

double GetAver( Ve_D& v );

void GetEDP( const char* seq, int endPos, int endPosition, Ve_D& EDP );

double GCContent( const char* seq, int begPos, int endPos );

Str toAminoSeq( const char* seq, int endPos, int endPosition);

void x2LongestORF( const char* seq,  Pa_I_I& location );

bool isSDSignal( Str& signal );
bool isTASignal( Str& signal );
bool defaultSignal( Str& signal );
double getGB( const char* seq, int beg, int end );
int getPrePhaseATG( const char* seq, int hint );

#include <string>
#include <sstream>
#include <iostream>
template <class T>

std::string to_string(T t, 
					  std::ios_base & (*f)(std::ios_base&) = std::dec)
{
   std::ostringstream oss;
   oss << f << t;
   return oss.str();
}

std::map<int,double> calculateDis2PreSTP(Ve_Location& locas);
#endif



