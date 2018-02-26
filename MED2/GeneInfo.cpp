//名称			：GeneInfo.cpp。
//最后更改日期	：2003年9月26号。
//描叙			：将原核序列取出ORF,保存正负序列。

#include "GeneInfo.h"
#include"SequenceTransform.h"
#include"OftenUsedOperatLib.h"

GeneInfo_T::GeneInfo_T( const con_Str& seq1 ) 
: positiveSeq( seq1 ), negtiveSeq( seq1 ), seqLen( seq1.size() )
, ORFSet( Ve_Location() ),background(4,0)
{
	std::reverse( negtiveSeq.begin(), negtiveSeq.end() );	//将序列翻转。
	std::for_each( negtiveSeq.begin(), negtiveSeq.end(),	//转换为互补链。
		SequenceTransform_T::ToOppRule_T() );

	int i = 0;
	for(; i < seqLen; ++i ){
		++background[positiveSeq[i]-'0'];
	}
	background[0] = background[3] 
		= (background[0]+background[3]) / seqLen / 2;
	background[1] = background[2] 
		= (background[1]+background[2]) / seqLen /2 ;
	cout<<"GC%: "<<background[1]*2<<endl;
	AABKG = calculateAminoFreq( background[0], background[1]
		,background[2],background[3] );
	if( background[1] > 0.28 )
		isGCRich = true;
}

GeneInfo_T::GeneInfo_T( const GeneSeq& geneSeq ):ORFSet( Ve_Location() )
,background(4,0)
{
	positiveSeq = geneSeq.positiveSeq;
	negtiveSeq = geneSeq.negtiveSeq;
	seqLen = geneSeq.seqLen;
		int i = 0;
		for(; i < seqLen; ++i ){
		++background[positiveSeq[i]-'0'];
	}
	background[0] = background[3] 
		= (background[0]+background[3]) / seqLen / 2;
	background[1] = background[2] 
		= (background[1]+background[2]) / seqLen /2 ;
	AABKG = calculateAminoFreq( background[0], background[1]
		,background[2],background[3] );

	if( background[1] > 0.28 )
		isGCRich = true;
}

Ve_Location& GeneInfo_T::getORFLocation()
{
	seq = positiveSeq.data();
	int i = 0;
	for(; i < 2; ++i )
	{
		if ( i != 0 )
			seq = negtiveSeq.data();
		getORFSet();
		if( i == 0 )
		{
			Ve_Location::iterator iter = ORFSet.begin();
			for( ; iter != ORFSet.end(); ++iter )
				iter->isPositive = true;
		}
	}
	return ORFSet;
}

void GeneInfo_T::getORFSet()
{
	int i = 0;
	for(; i < 3; ++i )
	{
		int hint = i;
		for(; hint < seqLen; )
		{
			Pa_I_I tmp = getNextPhaseORF( seq, hint );
			if( tmp.first == -1 )
				break;
			ORFSet.push_back( Location_T( tmp ) );
			hint = tmp.second + 3;
		}
	}
}

Pa_I_I GeneInfo_T::getNextPhaseORF( const char* seq, int hint )
{	
	int TIS =  getNextPhaseTIS( seq, hint );
	for(; ; )
	{
		if( TIS == -1 )
			return std::make_pair< int, int >( -1, -1 );
		int STP = getNextPhaseSTP( seq, TIS );
		if( STP - TIS  - 1 >= boundOfORF )
			return std::make_pair< int, int >( TIS, STP - 1 );
		else
			TIS =  getNextPhaseTIS( seq, TIS + 3 );
	}
}

int GeneInfo_T::getNextTISPosition( const char* seq, int hint )
{
	assert( seqLen != 0);
	//assert( seqLen != 0&&"测试序列的长度未初始化seqLen静态变量");
	for( ; ; hint += 1 )
	{
		if( hint > seqLen - 3 )
			return - 1;
		int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( subStr == 14 || subStr == 46 || subStr == 62 )//|| subStr == 30)
			return hint;
	}
}

int GeneInfo_T::getNextSTPPosition( const char* seq, int hint )
{
	for( ; ; hint += 1 )
	{
		if( hint > seqLen - 3 )
			return - 1;
		int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
		//			TAA-300			TGA-320			TAG-302
		if( subStr == 48 || subStr == 56 || subStr == 50 )
			return hint;
	}
}

int GeneInfo_T::getNextPhaseTIS( const char* seq, int hint )
{
	int TIS = hint;
	for(; ; )
	{
		TIS = getNextTISPosition( seq, TIS ) + 1;
		if( TIS == 0 )
			return -1;
		if( ( TIS - hint - 1 ) % 3 == 0 )
			return TIS - 1;
	}
}

int GeneInfo_T::getNextPhaseSTP( const char* seq, int hint )
{
	int STP = hint;
	for(; ;  )
	{
		STP = getNextSTPPosition( seq, STP ) + 1;
		if( STP == 0 )					
			return -1;
		if( (STP - 1 - hint) % 3 == 0 )		
			return STP;
	}
}
