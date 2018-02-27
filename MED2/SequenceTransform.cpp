#include "SequenceTransform.h"
#include <ctype.h>

Str genomeID;
Str SequenceTransform_T::char2DigitalSeq( Str& seq )
{
	int seqLen = seq.size();
	Str tmp;
	tmp.reserve( seqLen );
	int i = 0;
	for(; i < seqLen; ++i)
		tmp += char2digital( seq[i] );
	return tmp;
}

void SequenceTransform_T::reverseSeq( Str& seq )
{
	std::reverse( seq.begin(), seq.end() );
}

void SequenceTransform_T::char2DigitalSeqFile( Str& in, Str& out)
{
	ifstream inFile( in.data() );
	if( !inFile.good() ){
		std::cout<<"file "<<in<<" does not exist"<<std::endl;
		//std::cout<<"文件"<<in<<"不存在"<<std::endl;
		return;
	}

	Str descri;
	std::getline( inFile, descri );
	genomeID = descri.substr( descri.find_last_of("|") );

	Str tmp;
	tmp.reserve( 10000000 );
	while( !inFile.eof() )
	{
		char tmpChar;
		inFile>>tmpChar;
		tmp += char2digital( tmpChar );
	}
	tmp.erase( tmp.end() - 1 );
	inFile.close();
	
	int seqLen = tmp.size();
	ofstream outFile( out.data() );
	outFile.write( (char*)&seqLen, sizeof(int) );
	outFile.write( tmp.data(), tmp.size() );
	outFile.close();
}
	
char SequenceTransform_T::char2digital( char res )
{
	switch( toupper(res) )
	{
	case 'A' : return '0';
	case 'C' : return '1';
	case 'G' : return '2';
	case 'T' : return '3';
	case 'N' : return '2';
	case 'X' : return '0';
	case 'H' : return '3';
	case 'M' : return '1';
	case 'K' : return '2';
	case 'D' : return '0';
	case 'R' : return '2';
	case 'Y' : return '3';
	case 'S' : return '1';
	case 'W' : return '0';
	case 'B' : return '1';
	case 'V' : return '2';
	default  : 
		return '$';
	}
}

Str SequenceTransform_T::digital2CharSeq( Str& seq )
{
	Str tmp( seq.size(), '0' );
	int i = 0;
	for(; i < seq.size(); ++i )
	{
		switch( seq[i] )
		{
		case '0' : tmp[i] = 'a'; break;
		case '1' : tmp[i] = 'c'; break;
		case '2' : tmp[i] = 'g'; break;
		case '3' : tmp[i] = 't'; break;
		default :;
		}
	}
	return tmp;
}


void SequenceTransform_T::char2FileDigitalSeq( Str& in, Str& seq )
{
	seq.erase( seq.begin(), seq.end() );
	ifstream inFile( in.data() );
	if( !inFile.good() ){
		std::cout<<"file "<<in<<" not found!"<<std::endl;
		exit(1);
	}

	Str firstline;
	std::getline( inFile, firstline );
	if( firstline.find( "gb|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "gb|" ) + 3);
	if( firstline.find( "bj|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "bj|" ) + 3);
	if( firstline.find( "emb|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "emb|" ) + 4);

	Str GB, NM;
	GB = firstline.substr( 0, firstline.find_first_of( '|' ) );
	NM = firstline.substr( firstline.find_last_of( '|' ) + 1 );
	if( GB.find( NM.substr( 0, NM.find_first_of(' ')  ) ) 
		!= Str::npos ) 
		NM= NM.substr( NM.find_first_of(' ') + 1 );
	genomeID = GB + "|" + NM;

	seq.reserve( 10000000 );
		

	while( !inFile.eof() )
	{
		char tmpChar;
		inFile>>tmpChar;
		seq += char2digital( tmpChar );
	}
	seq.erase( seq.end() - 1 );
	inFile.close();
}
