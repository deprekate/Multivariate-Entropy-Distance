//a public lib in our team
#ifndef OFTENUSEDOPERATLIB_H
#define OFTENUSEDOPERATLIB_H

#include"TypeDef.h"

Ve_D calculateAminoFreq( double A, double C, double G, double T );

double getDistance				
( con_Ve_D& rhs, con_Ve_D& lhs );

int getMin( Ve_D lhs, con_Ve_D& coefficent = Ve_D( ) );

void addToPreVector( Ve_D& lhs, Ve_D& rhs, Ve_I counter = Ve_I( ) );

double GetDerivate( Ve_D& v );

double GetAver( Ve_D& v );

void GetEDP( const char* seq, int endPos, int endPosition, Ve_D& EDP );

double GCContent( const char* seq, int begPos, int endPos );

Str toAminoSeq( const char* seq, int endPos, int endPosition);

std::pair<Str,Str> strali(const char* stra, const char* strb);

std::pair<Ve_Pa_Str_I, int > findHostStr( const Str& source, const Str& target );

std::pair<Ve_Pa_Str_I, int > findHostStr( const char* seq, int begPos, int endPos, const Str& target );

void x2LongestORF( const char* seq, Pa_I_I& location );
#endif



