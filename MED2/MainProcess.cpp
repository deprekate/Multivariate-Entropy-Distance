#include"MED.h"
#include"SequenceTransform.h"

//write predicted CDS to files.
void resultToFile( std::vector<ORF_T>& result, Str& resultFilename );
//read genome seqeunces.
void readSeq( Str& seq, Str& seqFile );
//process overlapping genes.
void processOverlapping( Str file );
//process overlapping genes for GC rich genomes.
void processOverlappingForGCRich( Str file );

bool isGCRich = false;//will be updated automatically in Class::GeneInfo_T

int main( int argc, char* args[] 
		  )
{
	if( argc != 2 )
	{
		cout<<"usage: MED2 "
			"<GenomeSequenceFileName>"<<endl;
		return 1;
	}

	cout<<"Wellcome to MED2.0"<<endl;
	cout<<"Genome file is "<<args[1]<<endl;
	Ve_Str genomes;
	parafile = Str( "out.para" );
	//read genome sequence
	Str seq;
	Str _ar_(args[1]);
	SequenceTransform_T::char2FileDigitalSeq( _ar_, seq );
	std::cout<<"sequence size: "<<seq.size()<<std::endl;
	//MED2.0 initial CDS call
	MED geneFind( seq );
	std::vector<ORF_T> result = geneFind.getGeneLocation();

	//Write initial results to a file
	Str rstFile = Str(args[1]) + "out.MED";
	int pos = rstFile.find( "." );
	rstFile = rstFile.substr( 0, pos ) + "out.MED";
	resultToFile( result, rstFile );//*/
	//Read in intial results and process overlaping genes 
	if( isGCRich )
		processOverlappingForGCRich( rstFile );
	else
		processOverlapping( rstFile );
	//Call TIS model to revise starts annotation
	system( Str("./TISModel \"" 
		+Str(args[1]) + "\" \""+ rstFile+"\"").data() );
	return 1;
}

void readSeq( Str& seq, Str& seqFile )
{
	ifstream seqIn( seqFile.data(), std::ios::binary );
	assert( seqIn.good() );
	int seqLen;
	seqIn.read( (char*)&seqLen, sizeof(int) );
	seq.resize( seqLen, '0' );
	seqIn.read( const_cast<char*>( seq.data() ), seq.size() );
	seqIn.close();
}

void resultToFile( std::vector<ORF_T>& result, Str& resultFilename )
{
	ResultSortRule_T sortRule;
	std::sort( result.begin(), result.end(), sortRule );
	ofstream outResult( resultFilename.data() );
	int i = 0;
	for(; i < result.size(); ++i )
	{
		outResult<<setw(10)<<result[i].location.first
			<<std::setw(10)<<result[i].location.second<<std::setw(4)
			<<(result[i].isPositive ? '+' : '-' )
			<<setw(15)<<result[i].ATGScore
			<<setw(15)<<result[i].EDPScore<<setw(15)<<result[i].GBScore
			<<setw(15)<<result[i].certainCDS<<endl;
	}
	outResult.close();
}

class OverRecord{
public:
	OverRecord() : Onleft(0), OnRight(0), OnMid(0),isTotalOverlapped(false){
	}
	int Onleft;
	int OnRight;
	int OnMid;
	bool isTotalOverlapped;
};

struct Gene{
	char _strand;
	int _left;
	int _right;
	int _len;
	OverRecord localStrand;
	OverRecord oppositeStrand;
	double ATG;
	double EDP;
	double GB;
	double overPer;
	bool isCertainCDS;
	const bool operator < ( const Gene& rfh ) const{
		return _left < rfh._left;
	}
};

void processOverlapping( Str file ){
	cout<<"Process overlapping genes..."<<endl;
	ifstream in( file.data() );
	std::vector<Gene> genes;
	int maxLen = 0;	
	while( !in.eof() ){
		Gene gene;
		in>>gene._left>>gene._right>>gene._strand>>gene.ATG>>gene.EDP;
		gene._len = abs( gene._right - gene._left ) + 1;
		if( gene._len > maxLen )
			maxLen = gene._len;

		Str line;
		std::getline( in, line );
		gene.overPer = 0;
		genes.push_back( gene );
	}
	genes.erase( genes.end() - 1, genes.end() );
	in.close();
	Ve_D cris;
	double cr = 0.5;
	double b = 1.5;
	int i = 0;
	for(; ;++i){
		b -= 0.01;
		if( b < cr )
			break;
		cris.push_back(b);
	}
	int preNum = genes.size();
	int m = 0;
	for(; m < cris.size(); ++m ){
		std::sort( genes.begin(), genes.end() );//see MainProcess.err
		int pos = 0, neg = 0;
		int curGene = 1;
		for(; curGene < genes.size(); ++curGene ){
			{
				int iter;
				for(iter = curGene - 1
					; genes[curGene]._left - genes[iter]._left < maxLen 
					;  ){
					 if( --iter < 0)
						 break;
				};
				for( ++iter
					; genes[iter]._left - genes[curGene]._left < maxLen
					&& iter < genes.size(); ++iter ){
					if( iter != curGene ){
						if( genes[curGene]._strand == genes[iter]._strand ){
							if( genes[iter]._left <= genes[curGene]._left
								&& genes[iter]._right >= genes[curGene]._left
								&& genes[iter]._right <= genes[curGene]._right
								){
								int overNum = genes[iter]._right - genes[curGene]._left + 1;
								if( overNum > genes[curGene].localStrand.Onleft ) 
									genes[curGene].localStrand.Onleft = overNum;
							}
							if( genes[iter]._left < genes[curGene]._left
								&& genes[iter]._right > genes[curGene]._right
								){
								genes[curGene].localStrand.isTotalOverlapped = true;
							}
							if( genes[iter]._right > genes[curGene]._right
								&& genes[iter]._left > genes[curGene]._left
								&& genes[iter]._left < genes[curGene]._right
								){
								int overNum = genes[curGene]._right - genes[iter]._left + 1;
								if( overNum > genes[curGene].localStrand.OnRight )
									genes[curGene].localStrand.OnRight = overNum; 
							}
							if( genes[iter]._right < genes[curGene]._right
								&& genes[iter]._left > genes[curGene]._left ){
								if( genes[iter]._len > genes[curGene].localStrand.OnMid )
									genes[curGene].localStrand.OnMid = genes[iter]._len;
							}
						}
						else{
							if( genes[iter]._left < genes[curGene]._left
								&& genes[iter]._right > genes[curGene]._left
								&& genes[iter]._right < genes[curGene]._right
								){
								int overNum = genes[iter]._right - genes[curGene]._left + 1;
								if( overNum > genes[curGene].localStrand.Onleft ) 
									genes[curGene].oppositeStrand.Onleft = overNum;
							}
							if( genes[iter]._left < genes[curGene]._left
								&& genes[iter]._right > genes[curGene]._right
								){
								genes[curGene].oppositeStrand.isTotalOverlapped = true;
							}
							if( genes[iter]._right > genes[curGene]._right
								&& genes[iter]._left > genes[curGene]._left
								&& genes[iter]._left < genes[curGene]._right
								){
								int overNum = genes[curGene]._right - genes[iter]._left + 1;
								if( overNum > genes[curGene].localStrand.OnRight )
									genes[curGene].oppositeStrand.OnRight = overNum; 
							}
							if( genes[iter]._right < genes[curGene]._right
								&& genes[iter]._left > genes[curGene]._left ){
								if( genes[iter]._len > genes[curGene].localStrand.OnMid )
									genes[curGene].oppositeStrand.OnMid = genes[iter]._len;
							}
						}
					}
				}
			}
		}
		int number = genes.size();
		int i = 0;
		for(; i < genes.size(); ++i ){
			int overNum = 0;
			if( genes[i].oppositeStrand.isTotalOverlapped){
				overNum = genes[i]._len;
			}
			else
				overNum = genes[i].oppositeStrand.Onleft
				+ genes[i].oppositeStrand.OnRight + genes[i].oppositeStrand.OnMid;
			genes[i].overPer = overNum /double(genes[i]._len );
			if( genes[i].overPer >cris[m] ){
				++neg;
			
				genes.erase( (std::vector<Gene>::iterator)&genes[i] );
				--i;
			}
		}
		for( i = 0; i < genes.size(); ++i ){
			int overNum = 0;
			if( genes[i].localStrand.isTotalOverlapped){
				overNum = genes[i]._len;
			}
			else
				overNum = genes[i].localStrand.Onleft
				+ genes[i].localStrand.OnRight + genes[i].localStrand.OnMid;;
			double overPer = overNum /double(genes[i]._len );
			if( overPer > cris[m] ){
				++pos;
				genes.erase( (std::vector<Gene>::iterator)&genes[i] );
				--i;
			}
			if( genes[i].overPer > overPer ){
				genes[i].overPer *= -1;
			}
			else{
				genes[i].overPer = overPer;
			}
		}
		for( i = 0; i < genes.size(); ++i ){
			genes[i].localStrand = OverRecord();
			genes[i].oppositeStrand = OverRecord();
		}

	}
	ofstream out( file.data() );
	cout<<"Totally predict "<<genes.size()<<" genes."<<endl;
	for(  i = 0; i < genes.size(); ++i ){
		bool islongAndNoOverLap = fabs(genes[i].overPer)<0.0001 
			&& genes[i]._len > 300;
			out<<setw(15)<<genes[i]._left<<setw(15)<<genes[i]._right
				<<setw(15)<<genes[i]._strand
			<<setw(15)<<(islongAndNoOverLap?"Y":"N")<<endl;
	}
	out.close();
}

void processOverlappingForGCRich( Str file ){
	cout<<"Process overlapping genes..."<<endl;
	ifstream in( file.data() );
	std::vector<Gene> genes;
	int maxLen = 0;	
	if( !in.good() ){
		cout<<"file "<<file<<" not found"<<endl;
		exit(1);
	}
	while( !in.eof() ){
		Gene gene;
		in>>gene._left>>gene._right>>gene._strand>>gene.ATG>>gene.EDP>>gene.GB>>gene.isCertainCDS;
		gene._len = abs( gene._right - gene._left ) + 1;
		if( gene._len > maxLen )
			maxLen = gene._len;

		Str line;
		std::getline( in, line );
		gene.overPer = 0;
		genes.push_back( gene );
	}
	genes.erase( genes.end() - 1, genes.end() );
	in.close();
	Ve_D cris;
	double cr = 0.5;
	double b = 1;
	int i = 0;
	for(; ;++i){
		b -= 0.01;
		if( b < cr )
			break;
		cris.push_back(b);
	}
	int preNum = genes.size();
	int m = 0;
	for(; m < cris.size(); ++m ){
		do{
			int N = genes.size();
			std::vector<int> marks(genes.size(), 0);
			std::sort( genes.begin(), genes.end() );
			int i = 1;
			for(; i < genes.size(); ++i )
				genes[i].overPer = 0;
			for( i = 1; i < genes.size(); ++i ){
				Gene& first = genes[i-1];
				Gene& second = genes[i];
				int left = first._left > second._left ? first._left : second._left;
				int right = first._right < second._right ? first._right : second._right;
				double overLap = right - left;
				double firstOverP = overLap / first._len;
				double secondOverP = overLap / second._len;
				if( secondOverP > cris[m] ){
					if( first.EDP < second.EDP )
					{
						++marks[i];
					}
				}
				if( firstOverP > cris[m] ){
					if( first.EDP > second.EDP )
					{
						++marks[i-1];
					}
				}
				if( first._strand == second._strand ){
					if( fabs(first.overPer ) < firstOverP )
						first.overPer=firstOverP;
					if( fabs(second.overPer ) < secondOverP )
						second.overPer=secondOverP;
				}
				else{
					if( fabs(first.overPer ) < firstOverP )
						first.overPer=-firstOverP;
					if( fabs(second.overPer ) < secondOverP )
						second.overPer=-secondOverP;
				}
			}
			for( i = 0; i < genes.size(); ++i ){
				if( marks[i] > 0 
					){
					genes.erase( (std::vector<Gene>::iterator)&genes[i] );	
					marks.erase( (std::vector<int>::iterator)&marks[i] );
					--i;
				}
			}
			if( N == genes.size() )
				break;
		}while( true );
	}
	ofstream out( file.data() );
	cout<<"Totally predict "<<genes.size()<<" genes."<<endl;
	for(  i = 0; i < genes.size(); ++i ){
		bool islongAndNoOverLap = fabs(genes[i].overPer)<0.0001 
			&& genes[i]._len > 300;
		out<<setw(15)<<genes[i]._left<<setw(15)<<genes[i]._right
			<<setw(15)<<genes[i]._strand
			<<setw(15)<<(islongAndNoOverLap?"Y":"N")<<endl;
	}
	out.close();
}
