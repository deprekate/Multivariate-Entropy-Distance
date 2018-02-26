#include<algorithm>
#include<functional>
#include"SequenceTransform.h"
#include"MED_start.h"
#include"Classfier.h"
//outupt revised TIS
void resultToFile( std::vector<Location_T>& result, Str resultFilename, bool isLast = false );
//outupt revised TIS with more info.
void resultToFileFinal( std::vector<Location_T>& result, Str resultFilename );
//read genomic sequence
void readSeq( Str& seq, Str& seqFile );
//read initial gene coordinates
Ve_Location getRvsLocation( Str fileName );

bool isGCRich =false;//will be updated automatically in Class::GeneInfo_T

int main( int argc, char* args[]
		  )
{
	if( argc != 3)
	{
		cout<<"usage: TISModel <GenomeSequenceFileName>"
			" <locationFileName>"<<endl;
		return 1;
	}
	//*/
	cout<<"Second run of TIS model."<<endl;
	Str seq;
	//translate nucleatides to digital sequence.
	Str _ar_(args[1]);
	SequenceTransform_T::char2FileDigitalSeq( _ar_, seq );

	GeneSeq geneSeq( seq );

	Str resultFilename = Str( args[2] );

	double ATNum = 0;
	int i = 0;
	for(; i < seq.size(); ++i ){
		if( seq[i] == '0' || seq[i] == '3' )
			++ATNum;
	}
	Ve_D bg(4,0);
	bg[0] = bg[3] = ATNum / seq.size() / 2;
	bg[1] = bg[2] = 0.5 - bg[0];

	Ve_Location rst;
	int iterNum = 0;
	Str finalParaFileName;
	do{
		rst = getRvsLocation( resultFilename );
		std::map<int,double> predis = calculateDis2PreSTP( rst );
		Ve_Str paraFiles;
		Ve_Str tags;
		cout<<"Extract TUI and TUL genes (length >300bp)"<<endl;
		std::map<Str,Ve_Location> rvsLocations 
			= classGene( rst, bg);	
		std::map<Str,Ve_Location>::iterator iter = rvsLocations.begin();
		std::vector<double> modelType;
		cout<<"Compare TUI TIS model with TUL TIS model..."<<endl;
		for( ; iter != rvsLocations.end(); ++iter ){
			cout<<"Train "<<iter->first<<" TIS model";
			Ve_Location tmp = iter->second;
			MED_start med_start=MED_start( geneSeq, tmp, CParametersDlg(4, 15) );

			modelType.push_back(med_start.getRvsLocations(2));
		}
		if( modelType[0] * modelType[1] > 0 ){
			cout<<"Lamda > 0: TUL TIS model is similar to TUI TIS model."<<endl;
			cout<<"Re-train TIS model";
			tags = Ve_Str(1,"" );
			parafile = resultFilename.substr
				(0, resultFilename.find("."))+".paratmp";
			paraFiles = Ve_Str(1,parafile);
			MED_start med_start( geneSeq, rst, CParametersDlg(4, 15) );

			med_start.getRvsLocations(20);
		}
		else
		{
			cout<<"Lamda < 0: TUL TIS model differs from TUI TIS model"<<endl;
			std::map<Str,Ve_Location>::iterator iter = rvsLocations.begin();
			for( ; iter != rvsLocations.end(); ++iter ){
				cout<<"Re-train "<<iter->first<<" TIS model";
				tags.push_back(iter->first);
				parafile = resultFilename.substr
				(0, resultFilename.find("."))+iter->first+".paratmp";
				paraFiles.push_back( parafile );
				MED_start med_start
					( geneSeq, iter->second, CParametersDlg(4, 15) );

				med_start.getRvsLocations(20);
			}
		}
		std::vector<MED_start> med_starts;
		std::vector<Ve_Location> rsts;
		finalParaFileName = resultFilename.substr
				(0, resultFilename.find(".") )+".para";
		std::ofstream out(finalParaFileName.data());
		out<<endl;
		out.close();//clear previous record if file already exist.
		for( i = 0; i < paraFiles.size(); ++i ){
			MED_start med_start( geneSeq, rst );
			med_start.readParametersFromFile( paraFiles[i].data() );
			system(Str("rm " + paraFiles[i]).data());
			Ve_Location rstTmp = rst;
			med_start.reviseGivenORFs( rstTmp );
			med_start.outReadableParameters( finalParaFileName
				, tags[i] );
			rsts.push_back( rstTmp );
		}
		cout<<"Choose final TIS..."<<endl;
		int j = 0;
		for(; j < rsts.size(); ++j ){
			Ve_D scores;
			int m = 0;
			for(; m < rsts[j].size(); ++m ){
				if( fabs(rsts[j][m].tisScore)  <100000 ){
					scores.push_back( rsts[j][m].tisScore );
				}
			}
			double aver = GetAver( scores ), sd = GetDerivate( scores );
			for( m = 0; m < rsts[j].size(); ++m ){
				rsts[j][m].tisScore = (rsts[j][m].tisScore - aver )	/sd ;
			}
		}
		rst.erase(rst.begin(), rst.end());

		for( j = 0; j < rsts[0].size(); ++j ){
			int maxIndex = -1;
			double maxScore = -1000000;
			int i = 0;
			for(; i < rsts.size(); ++i ){
				if( rsts[i][j].tisScore > maxScore ){	
					maxIndex = i;
					maxScore = rsts[i][j].tisScore;
				}
			} 
			if( rsts.size() == 2 ){
				if( fabs(rsts[maxIndex][j].tisScore 
					- rsts[::fabs(1-maxIndex)][j].tisScore) < 0.1 )
					if( rsts[maxIndex][j].isPositive )
						maxIndex = rsts[maxIndex][j].location.first 
						< rsts[::fabs(1-maxIndex)][j].location.first 
						? maxIndex : ::fabs(1-maxIndex);
					else
						maxIndex = rsts[maxIndex][j].location.second 
						> rsts[::fabs(1-maxIndex)][j].location.second
						? maxIndex : ::fabs(1-maxIndex);
			}
			rst.push_back(rsts[maxIndex][j]);
		}
		resultToFile( rst, resultFilename );
	}while( ++iterNum != 1 ); 
	resultToFile( rst, resultFilename, true );
	resultToFileFinal( rst, finalParaFileName );
	return 1;
}

void readSeq( Str& seq, Str& seqFile )
{
	std::ifstream seqIn( seqFile.data(), std::ios::binary );
	assert( seqIn.good() );
	int seqLen;
	seqIn.read( (char*)&seqLen, sizeof(int) );
	seq.resize( seqLen, '0' );
	seqIn.read( const_cast<char*>( seq.data() ), seq.size() );
	seqIn.close();
}

void resultToFile( std::vector<Location_T>& result, Str resultFilename, bool isLast )
{
	if( isLast )
		cout<<"Predicted genes have been saved in file: "<<resultFilename<<endl;
	ResultLSortRule_T sortRule;
	std::sort( result.begin(), result.end(), sortRule );
	std::ofstream outResult( resultFilename.data() );

	int  i = 0;
	for(; i < result.size(); ++i )
	{
				outResult<<std::setw(10)<<(result[i].location.first)
					<<std::setw(10)<<result[i].location.second
					<<std::setw(2)
					<<(result[i].isPositive ? '+' : '-' );
				if( !isLast )
					outResult<<std::setw(5)<<(result[i].isLongAndNonOverlap ? 'Y' : 'N' )<<endl;
				else
					outResult<<'\n';
	}
	outResult.close();
}

void resultToFileFinal( std::vector<Location_T>& result, Str resultFilename )
{
	ResultLSortRule_T sortRule;
	std::sort( result.begin(), result.end(), sortRule );
	std::ofstream outResult( resultFilename.data(), std::ios::app );
	cout<<"Detailed information for predicted TIS has been saved in "
		<<resultFilename<<endl;
	outResult<<"====================================="<<endl;
	outResult<<"Detailed information for predicted TIS"<<endl;
	outResult<<std::setw(10)<<"rhs\'"
			<<std::setw(10)<<"lhs\'"
			<<std::setw(8)
			<<"Strand"
			<<std::setw(8)<<"Signal"
			<<std::setw(8)<<"Spacer"
			<<std::setw(5)<<"XTG"
			<<std::setw(8)<<"Index"<<endl;
	int  i = 0;
	for(; i < result.size(); ++i )
	{
		outResult<<std::setw(10)<<(result[i].location.first)
			<<std::setw(10)<<result[i].location.second
			<<std::setw(8)
			<<(result[i].isPositive ? '+' : '-' )
			<<std::setw(8)<<result[i].Sig;
		if( result[i].Sig == "-" ) 
			outResult<<std::setw(8)<<"-";
		else
			outResult<<std::setw(8)<<abs(result[i].SigPos);
		outResult<<std::setw(5)<<result[i].TISCode
			<<std::setw(8)<<result[i].ATGIndexInORF<<endl;
	}
	outResult.close();
}

Ve_Location getRvsLocation( Str fileName ){
	Ve_Location results;
	std::ifstream in( fileName.data() );
	bool isTaged = false;
	while( !in.eof() )
	{
		Location_T tmp;
		Str c, d;
		in>>tmp.location.first>>tmp.location.second>>c>>d;
		string tag;
		std::getline( in, tag );

		if( d == "N" && isGCRich )
			tmp.isLongAndNonOverlap =false;

		tmp.ORFLength = ::fabs( tmp.location.first - tmp.location.second );
		if( c.empty() )
			break;
		if( c == "+" )
			tmp.isPositive = true;
		else
			tmp.isPositive = false;
		tmp.preSeq = SequenceTransform_T::char2digitalSeq(tmp.preSeq);
		results.push_back( tmp );
	}
	return results;
}

