#ifndef SEQUENCETRANSFORM_H
#define SEQUENCETRANSFORM_H

#include"TypeDefBase.h"
extern Str genomeID;
class SequenceTransform_T
{
public:
	static void char2FileDigitalSeq( Str& in, Str& seq );
	static Str digital2CharSeq( Str& seq );
	static Str char2digitalSeq( Str& seq );
	class ToOppRule_T
	{
	public:
		void operator ()( char& c )
		{
				switch( c )
				{
					case '0' : c = '3'; break;
					case '1' : c = '2'; break;
					case '2' : c = '1'; break;
					case '3' : c = '0'; break;
					default  : ;
				}
		}
	};
private:
	SequenceTransform_T()
	{
	}
	static char char2digital( char resourse );	
};

#endif// SEQUENCETRANSFORM_H



