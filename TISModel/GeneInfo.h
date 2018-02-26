#ifndef GENEINFO_H
#define GENEINFO_H
#include"TypeDef.h"
#include"GeneSeq.h"

class GeneInfo_T						
{
public:
	GeneInfo_T( const con_Str& seq );
	GeneInfo_T( const GeneSeq& geneSeq );
	GeneInfo_T() : seq(NULL), seqLen(0)
	{
	}
	const Str& getPositiveSeq()const		
	{
		return positiveSeq;
	}
	const Str& getNegtiveSeq()const			
	{
		return negtiveSeq;
	}
	int getSeqLength()
	{
		return seqLen;
	}
protected:
	int getNextTISPosition( const char* seq, int hint = 0 );
	int getNextPhaseTIS( const char* seq, int hint = 0 );
	const char* seq;
	Str negtiveSeq;								
	Str positiveSeq;							
	int seqLen;									
	Ve_D background;
private:
	enum { boundOfORF = 90 };					
};

#endif// GENEINFO_H



