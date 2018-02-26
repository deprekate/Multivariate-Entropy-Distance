#include"OftenUsedOperatLib.h"

bool isSDSignal( Str& signal ){
	int i = 0;
	for(; i < signal.size(); ++i ){
		if( signal[i] == '1' )
			return false;
	}
	return true;
}
int getPrePhaseATG( const char* seq, int hint ){
	int iter = hint;
	for( iter -= 3; ; iter -= 3 )
	{
		if( iter < 0 )
		{
			return hint;
		}
		
		int subStr = ((seq[iter] - 48)<<4) + ((seq[iter+1] - 48)<<2) + seq[iter + 2] - 48;
		if( subStr == 48 || subStr == 56 
			|| subStr == 50 ){
			return hint;
		}
		if( subStr == 14 || subStr == 46 || subStr == 62 )//|| subStr == 30)
			return iter;
	}
}

double getGB( const char* seq, int beg, int end ){
	double GB13 = 0, GB123 = 0;
	int index =0 ;
	int i = beg;
	for(; i < end  - 6; ++i, ++index ){
		if( seq[i] == '1' || seq[i] == '2' ) {
			if( index % 3 == 0 || index % 3 == 2 )
				++GB13;
			++GB123;
		}
	}
	return 100*GB13/GB123;
}

bool isTASignal( Str& signal ){
	int i = 0;
	for(; i < signal.size(); ++i ){
		if( signal[i] == '1' || signal[i] == '2' )
			return false;
	}
	return true;
}
bool defaultSignal( Str& signal ){
	return true;
}
double GetAver( Ve_D& v )
{
	double aver = 0;
	int i = 0;
	for(; i < v.size(); ++i )
		aver += v[i];
	return (aver / v.size());
}

double GetDerivate( Ve_D& tmp, bool perio  )
{
	Ve_D v;
	int k = 0;
	for(; k < tmp.size()-3; ++k ){
		v.push_back( tmp[k] + tmp[k+1] + tmp[k+2]);
	}
	if( !perio )
		v = tmp;
	double aver = GetAver( v );
	double der = 0;
	int i = 0;
	for(; i < v.size(); ++i )
		der += ( v[i] - aver ) * ( v[i] - aver );
	return sqrt( der / v.size() );
}

void x2LongestORF( const char* seq, Pa_I_I& location )
{
	int STPPos = location.second;
	int ATGPos = -1;
	int tmp = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;

	for( STPPos -= 3; ; STPPos -= 3 )
	{
		if( STPPos < 0 )
		{
			location.first = ATGPos;
			return;
		}
		
		int subStr = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;
		//			TAA-300			TGA-320			TAG-302
		if( subStr == 48 || subStr == 56 || subStr == 50 ){
			location.first = ATGPos;
			return;
		}
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( subStr == 14 || subStr == 46 || subStr == 62 )//|| subStr == 30)
			ATGPos = STPPos;
	}
	location.first = ATGPos;
}

std::map<int,double> calculateDis2PreSTP(Ve_Location& locas){
 	std::map<int,double> lss;
	Ve_Location NL, PL;
 	int i = 0;
 	for(; i < locas.size(); ++i ){
 		if( locas[i].isPositive )
 			PL.push_back( locas[i] );
 		else
 			NL.push_back( locas[i] );
 	}
 	ResultLSortRule_T sortRule;
 	std::sort( PL.begin(), PL.end(), sortRule );
 	for( i = 1; i < PL.size(); ++i ){
 		PL[i].dis2PreSTP = PL[i].location.first
 			- PL[i-1].location.second;
		++lss[PL[i].dis2PreSTP];
 	}
 	std::sort( NL.begin(), NL.end(), sortRule );
 	for(  i = 0; i < NL.size()-1; ++i ){
 		NL[i].dis2PreSTP = NL[i+1].location.first
 			- NL[i].location.second;
		++lss[NL[i].dis2PreSTP];
 	}
 	locas = PL;
 	locas.insert( locas.end(), NL.begin(), NL.end() );
	return lss;
}
