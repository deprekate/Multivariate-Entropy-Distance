
#include"SequenceTransform.h"
	
char SequenceTransform_T::char2digital( char res )
{
	switch( res)
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


Str SequenceTransform_T::digital2CharSeq( Str& seq )
{
	Str tmp( seq.size(), '0' );
	int i = 0;
	for(; i < seq.size(); ++i )
	{
		switch( seq[i] )
		{
		case '0' : tmp[i] = 'A'; break;
		case '1' : tmp[i] = 'C'; break;
		case '2' : tmp[i] = 'G'; break;
		case '3' : tmp[i] = 'T'; break;
		default :;
		}
	}
	return tmp;
}

Str SequenceTransform_T::char2digitalSeq( Str& seq ){
	Str tmp( seq.size(), '0' );
	int i = 0;
	for(; i < seq.size(); ++i ){
		tmp[i] = char2digital( seq[i]);
	}
	return tmp;
}
