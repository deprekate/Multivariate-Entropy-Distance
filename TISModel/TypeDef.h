#ifndef TYPEDEF_H
#define TYPEDEF_H

#include"TypeDefBase.h"

#include"Matrix.h"
typedef Matrix_T<double> M_D;
typedef std::vector<Matrix_T<double> > Ve_M_D;
typedef std::pair< Ve_D, M_D > Pa_Ve_D_M_D;
extern bool isGCRich;
class Location_T			
{
public:
	Location_T() : isPositive( false ), ATGIndexInORF( 0 )
		,ORFLength(0 ), SDPos(0), TAPos(0), tisScore(0), isLongAndNonOverlap(true)
	{
		anotherTIS = -1;
	}
	Location_T( Pa_I_I& rhs, bool positive = false ) 
		: isPositive( positive ), ATGIndexInORF( 0 ), ORFLength(0 )
		, SDPos(0), TAPos(0), tisScore(0),isLongAndNonOverlap(true)
	{
		location = rhs;
		anotherTIS = -1;
	}

	const bool operator <( const Location_T& rhs )const
	{
		return location < rhs.location;
	}

	const bool operator ==( const Location_T& rhs )const
	{
		if( isPositive != rhs.isPositive )
			return false;
		return isPositive ? location.second == rhs.location.second
			: location.first == rhs.location.first;
	}
	
	void location2FileForm( int seqLen )
	{
		if( isPositive )
		{
			location.first += 1;
			location.second += 3;
		}
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second ) - 2;
			location.second = abs( seqLen - tmp );
		}
	}

	void fileForm2Location( int seqLen )
	{
		if( isPositive )
		{
			location.first -= 1;
			location.second -= 3;
		}
	
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second );
			location.second = abs( seqLen - tmp ) - 2;
		}
	}

	Pa_I_I location;//|A|TG........|T|TA¡£
	bool isPositive, isLongAndNonOverlap;
	Str TISCode;
	int ATGIndexInORF;
	int ORFLength;
	int SDPos,TAPos, SigPos;
	Str SDSig, TASig, Sig;
	int dis2PreSTP;
	double TAScore, SDScore, ATGScore, disScore;
	double sigScore, tisScore;
	Str Class;
	Str preSeq;
	Str seqAroundTA, seqAroundSD, seqAroundTIS;
	int anotherTIS;
	double anotherTISScore;
};

class ResultLSortRule_T
{
public:
	bool operator () ( Location_T lhs, Location_T rhs )
	{
		return lhs.location.first < rhs.location.first;
	}
};

typedef std::vector< Location_T > Ve_Location; 
typedef const std::vector< Location_T > con_Ve_Location; 
typedef std::vector< Location_T > Ve_Loca;  
typedef const std::vector< Location_T > con_Ve_Loca; 
typedef std::list< Location_T > Li_Loca;  
typedef std::set< Location_T > Set_Loca;
using std::cout;
using std::endl;
using std::setw;
using std::make_pair;
#endif //TYPEDEF_H
