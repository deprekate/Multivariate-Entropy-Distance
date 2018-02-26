//File name		:OftenUsedOperatLib.cpp
//Last modified	:27,7,2003
//Discription	:realize the functions defined in the file :OftenUsedOperatLib.cpp.

#include"OftenUsedOperatLib.h"

double getDistance
( con_Ve_D& rhs, con_Ve_D& lhs )
{
	double distance = 0;
	assert(rhs.size()==lhs.size());
	int i=0;
	for(;i<rhs.size();++i)
		distance+=pow(rhs[i]-lhs[i],2);
	return distance;
}

int getMin( Ve_D lhs, con_Ve_D& coefficent )
{
	int minPosition=-1;
	if(coefficent.size()!=0)
	{
		assert(coefficent.size()==lhs.size());
		int i=0;
		for(;i<lhs.size();++i)
			lhs[i]*=coefficent[i];
	}

	double minElement=lhs[0];
	int i=0;
	for(;i<lhs.size();++i)
	{
		if(lhs[i]<=minElement)
		{
			minElement=lhs[i];
			minPosition=i;
		}
	}
	return minPosition;
}

void addToPreVector( Ve_D& lhs, Ve_D& rhs, Ve_I& counter )
{
	assert(rhs.size()==lhs.size()&&(lhs.size()==counter.size()||counter.size()==0));
	int i=0;
	for(;i<rhs.size();++i)
	{
		lhs[i]+=rhs[i];
		if(counter.size()!=0)
			++counter[i];
	}
}

double GetAver( Ve_D& v )
{
	double aver = 0;
	int i = 0;
	for(; i < v.size(); ++i )
		aver += v[i];
	return aver / v.size();
}

double GetDerivate( Ve_D& v )
{
	double aver = GetAver( v );
	double der = 0;
	int i = 0;
	for(; i < v.size(); ++i )
		der += ( v[i] - aver ) * ( v[i] - aver );
	return sqrt( der / (v.size()-1) );
}


void GetEDP( const char* seq, int startPos, int endPosition, Ve_D& EDP )
{
	double sum = 0;
	assert( EDP.size() == 20 );
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )
	{
		++sum;
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 )
		{
		case 36 : case 37 : case 38 : case 39 :
			++EDP[0]; break;	//tmp += 'A'; break;
		case 57 : case 59 :
			++EDP[1]; break;	//tmp += 'C'; break;
		case 33 : case 35 :
			++EDP[2]; break;	//tmp += 'D'; break;
		case 32 : case 34 :
			++EDP[3]; break;	//tmp += 'E'; break;
		case 61 : case 63 :
			++EDP[4]; break;	//tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			++EDP[5]; break;	//tmp += 'G'; break;
		case 17 : case 19 :
			++EDP[6]; break;	//tmp += 'H'; break;
		case 12 : case 13 : case 15 : 
			++EDP[7]; break;	//tmp += 'I'; break;
		case 0 : case 2 : 
			++EDP[8]; break;	//tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			++EDP[9]; break;	//tmp += 'L'; break;
		case 14 : 
			++EDP[10]; break;//tmp += 'M'; break;
		case 1 : case 3 : 
			++EDP[11]; break;//tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			++EDP[12]; break;//tmp += 'P'; break;
		case 16 : case 18 : 
			++EDP[13]; break;//tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			++EDP[14]; break;//tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :  
			++EDP[15]; break;//tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 : 
			++EDP[16]; break;//tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			++EDP[17]; break;//tmp += 'V'; break;
		case 58 :
			++EDP[18]; break;//tmp += 'W'; break;
		case 49 : case 51 :
			++EDP[19]; break;//tmp += 'Y'; break;
		case 48 : case 50 : case 56 :
		{
			--sum;
			break;
		}
		default :;
		}
	}

	int j = 0;
	for(; j < 20; ++j )
	{
		if( EDP[j] != 0 )
		{
			double pi = EDP[j] / sum;
			EDP[j] = pi * log( pi );
		}
	}
	sum = 0;
	for( j = 0; j < 20; ++j )
		sum += EDP[j];

	if( sum != 0 )
	{
		for( j = 0; j < 20; ++j)		
			EDP[j] /= sum;
	}
}

double GCContent( const char* seq, int begPos, int endPos ){
	int num = 0;
	int i = begPos;
	for(; i<= endPos; ++i )
		if( seq[i] == '1' || seq[i] == '2' )
			++num;
		return num/double( endPos - begPos + 1 );
}

Str toAminoSeq( const char* seq, int startPos, int endPosition){
	double sum = 0;
	Str tmp;
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )
	{
		++sum;
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 )
		{
		case 36 : case 37 : case 38 : case 39 :
			tmp += 'A'; break;
		case 57 : case 59 :
			tmp += 'C'; break;
		case 33 : case 35 :
			tmp += 'D'; break;
		case 32 : case 34 :
			tmp += 'E'; break;
		case 61 : case 63 :
			tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			tmp += 'G'; break;
		case 17 : case 19 :
			tmp += 'H'; break;
		case 12 : case 13 : case 15 : 
			tmp += 'I'; break;
		case 0 : case 2 : 
			tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			tmp += 'L'; break;
		case 14 : 
			tmp += 'M'; break;
		case 1 : case 3 : 
			tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			tmp += 'P'; break;
		case 16 : case 18 : 
			tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :  
			tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 : 
			tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			tmp += 'V'; break;
		case 58 :
			tmp += 'W'; break;
		case 49 : case 51 :
			tmp += 'Y'; break;
		case 48 : case 50 : case 56 :
		{
			--sum;
			break;
		}
		default :;
		}
	}
	return tmp;
}

std::pair<Ve_Pa_Str_I, int > findHostStr( const Str& source, const Str& target )
{
	int ss = source.size(), ts = target.size(); 
	Ve_Pa_Str_I _vpsi_;
	std::pair<Ve_Pa_Str_I, int > hostInfo( _vpsi_, ts) ;
	int i = 0;
	for(; i <= ss - ts; ++i )
	{
		int mutNum = 0;
		int j = 0;
		for(; j < ts; ++j )
		{
			if( source[i + j] != target[j] )
				++mutNum;
		}
		if( mutNum < hostInfo.second ){
			hostInfo.first = Ve_Pa_Str_I( 1, std::pair<Str,int>( source.substr( i, ts ), i ) );
			hostInfo.second = mutNum;
		}
		else if( mutNum == hostInfo.second ){
			hostInfo.first.push_back( std::pair<Str,int>( source.substr( i, ts ), i ) );
		}
	}
	return hostInfo;
}

std::pair<Ve_Pa_Str_I, int > findHostStr( const char* source, int begPos, int endPos, const Str& target )
{
	int ts = target.size(); 
	int len = endPos - begPos + 1;
	Ve_Pa_Str_I _vpsi_;
	std::pair<Ve_Pa_Str_I, int > hostInfo( _vpsi_, ts) ;
	int i = begPos;
	for(; i <= endPos - ts && source[i] != '\0'; ++i )
	{
		int mutNum = 0;
		int j = 0;
		for(; j < ts; ++j )
		{
			if( source[i + j] != target[j] )
				++mutNum;
		}
		if( mutNum < hostInfo.second ){
			hostInfo.first = Ve_Pa_Str_I( 1, std::pair<Str,int>( Str( &source[i], &source[i] + ts ), i - begPos - len) );
			hostInfo.second = mutNum;
		}
		else if( mutNum == hostInfo.second ){
			hostInfo.first.push_back( std::pair<Str,int>( Str( &source[i], &source[i] + ts ), i - begPos -len ) );
		}
	}
	return hostInfo;
}

Ve_D calculateAminoFreq( double A, double C, double G, double T )
{
	double FN[4] = { A, C, G, T} ,F[4][4][4];
	Ve_D N(21,0);
	int i,j,k;
	i=j=k=0;

	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			for(k=0;k<4;k++)
				F[i][j][k]=FN[i]*FN[j]*FN[k];


	//A
	N[0]=F[2][1][0]+F[2][1][1]+F[2][1][2]+F[2][1][3];
	//C
	N[1]=F[3][2][3]+F[3][2][1];

	//D    
	N[2]=F[2][0][1]+F[2][0][3];

	//E
	N[3]=F[2][0][0]+F[2][0][2];
	//F
	N[4]=F[3][3][3]+F[3][3][1];
	//G
	N[5]=F[2][2][0]+F[2][2][1]+F[2][2][2]+F[2][2][3];
	//H
	N[6]=F[1][0][1]+F[1][0][3];
	//I
	N[7]=F[0][3][3]+F[0][3][1]+F[0][3][0];
	//K
	N[8]=F[0][0][0]+F[0][0][2];
	//L
	N[9]=F[3][3][0]+F[3][3][2]+F[1][3][0]+F[1][3][1]+F[1][3][2]+F[1][3][3];
	//M
	N[10]=F[0][3][2];
	//N
	N[11]=F[0][0][1]+F[0][0][3];
	//P
	N[12]=F[1][1][0]+F[1][1][1]+F[1][1][2]+F[1][1][3];
	//Q
	N[13]=F[1][0][2]+F[1][0][0];
	//R
	N[14]=F[1][2][0]+F[1][2][1]+F[1][2][2]+F[1][2][3]+F[0][2][0]+F[0][2][2];
	//S
	N[15]=F[3][1][0]+F[3][1][1]+F[3][1][2]+F[3][1][3]+F[0][2][3]+F[0][2][1];
	//T
	N[16]=F[0][1][0]+F[0][1][1]+F[0][1][2]+F[0][1][3];
	//V
	N[17]=F[2][3][0]+F[2][3][1]+F[2][3][2]+F[2][3][3];
	//W
	N[18]=F[3][2][2];
	//Y
	N[19]=F[3][0][1]+F[3][0][3];
	//STP
	N[20]=F[3][0][0]+F[3][2][0]+F[3][0][2];
	return N;
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
		if( subStr == 48 || subStr == 56 
			|| subStr == 50 ){
			location.first = ATGPos;
			return;
		}
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( subStr == 14 || subStr == 46 || subStr == 62 )//|| subStr == 30)
			ATGPos = STPPos;
	}
	location.first = ATGPos;
}



