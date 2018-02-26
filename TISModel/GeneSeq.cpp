// GeneSeq.cpp: implementation of the GeneSeq class.
//
//////////////////////////////////////////////////////////////////////

#include "GeneSeq.h"
#include"SequenceTransform.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GeneSeq::GeneSeq( const con_Str& seq1 )
:positiveSeq( seq1 ), negtiveSeq( seq1 ), seqLen( seq1.size() )
{
	std::reverse( negtiveSeq.begin(), negtiveSeq.end() );
	std::for_each( negtiveSeq.begin(), negtiveSeq.end(),
		SequenceTransform_T::ToOppRule_T() );
	double cNum = 0;
	int i = 0;
	for(; i < seqLen; ++i ){
		if( positiveSeq[i] == '1' )
			++cNum;
	}
	if( (cNum / seqLen) > 0.28 )
		isGCRich = true;
}


