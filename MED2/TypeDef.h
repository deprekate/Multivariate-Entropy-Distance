#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <string>
#include"TypeDefBase.h"

#include"Matrix.h"
typedef Matrix_T<double> M_D;
typedef std::vector<Matrix_T<double> > Ve_M_D;
typedef std::pair< Ve_D, M_D > Pa_Ve_D_M_D;
extern bool isGCRich;
class Location_T								
{
public:
	Location_T() : isPositive( false ),ATGIndexInORF( 0 ), 
		RBSScore( -100000 ), ORFLength( 0 ), MinDis( 10000 ),LongestORFLength(0),
		motifWM(0), motifSpacer(0), XTGWM(0), XTGIndex(0), RBSDis(0), GC(0), tag("A"),
		canBeCenter(true), certainCDS(false)
	{
	}
	Location_T( Pa_I_I& rhs, bool positive = false ) 
		: isPositive( positive ),ATGIndexInORF(0),LongestORFLength(0),
		motifWM(0), motifSpacer(0), XTGWM(0), XTGIndex(0),RBSDis(0), GC(0), tag("A"),
		canBeCenter(true), certainCDS(false)
	{
		location = rhs;
	}

	const bool operator <( const Location_T& rhs )const
	{
		return location.second < rhs.location.second;
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
	bool isValidate()
	{
		return motifWM!=0 && motifSpacer!=0 && XTGWM!=0 && XTGIndex!= 0 && RBSScore !=-100000;
	}
	Pa_I_I location;//|A|TG........|T|TA¡£
	bool isPositive;
	bool canBeCenter;
	bool certainCDS;
	int ATGIndexInORF;
	double RBSScore;
	double RBSDis;
	double EDPScore;
	double ATGScore;
	double ORFLength;
	double LongestORFLength;
	double MinDis;
	double motifWM, motifSpacer, XTGWM, XTGIndex;
	double MEDDis2T, MEDDis2F;
	double GC;
	double TDis, FDis;
	double GBScore, preATGGBScore, posATGGBScore;
	Str tag;
	Str otherInfo;
	Str hitMotif, motif;
	int ATGNumCounted;
	int motifPosition;
	Str codingSeq;
	std::pair<Ve_Pa_Str_I, int > hosters;
	Str SDStr;
};

	class ORF_T : public Location_T
	{
	public:
		ORF_T():EDP(20,0), Dis2TCenters(100000), Dis2FCenters(100000)
		{
		}
		Ve_D EDP;
		double Dis2TCenters;
		double Dis2FCenters;
	};

	class ResultSortRule_T
{
public:
	bool operator () ( ORF_T lhs, ORF_T rhs )
	{
		return lhs.location.first < rhs.location.first;
	}
};

class ResultLSortRule_T
{
public:
	bool operator () ( Location_T lhs, Location_T rhs )
	{
		return lhs.location.first < rhs.location.first;
	}
};

class SortByLenth
{
public:
	bool operator () ( Location_T lhs, Location_T rhs )
	{
		return lhs.ORFLength < rhs.ORFLength;
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
