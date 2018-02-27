#include "GeneInfo.h"
#include"SequenceTransform.h"

GeneInfo_T::GeneInfo_T( const con_Str& seq1 ) 
: positiveSeq( seq1 ), negtiveSeq( seq1 ), seqLen( seq1.size() )
{
	std::reverse( negtiveSeq.begin(), negtiveSeq.end() );
	std::for_each( negtiveSeq.begin(), negtiveSeq.end(),
		SequenceTransform_T::ToOppRule_T() );
}

GeneInfo_T::GeneInfo_T( const GeneSeq& geneSeq )
:background(4,0)
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
	if( background[1] > 0.28 )
		isGCRich = true;
}

int GeneInfo_T::getNextTISPosition( const char* seq, int hint )
{
	assert( seqLen != 0 && "The length of the test sequence is not initialized seqLen static variables");
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

