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
}

GeneSeq::~GeneSeq()
{

}
