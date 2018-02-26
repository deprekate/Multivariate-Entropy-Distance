#include"MED_start.h"
#include"SequenceTransform.h"
#include"OftenUsedOperatLib.h"
#include"time.h"
#pragma warning (disable:4786)
Str parafile;
 int postmp;
MED_start::MED_start
( GeneSeq& geneseq, Ve_Location& rvsLocations, CParametersDlg parameters )
	: GeneInfo_T( geneseq ), XTGIndex( 6, 1 ), XTG( 4, 0.01 )
{
	setRvsLocation( rvsLocations );
	setParameters( &parameters );
	ATGWM = M_D( 4, preATGBP + posATGBP + 3,0.0001 );
}

double MED_start::findMotifs(int mutBp, Ve_Str& regions
							 , int regionLen, Str file, int limitNum)
{
	Ve_I initSignal( sigLen, 0 );
	initSignal[sigLen-1] = -1;
	
	int upBoundary = pow( 4, sigLen );
	ofstream out( file.data() );
	//
	Ve_D sig;
	long double nonrandom = 0;
	int i = 0;
	for(; i < upBoundary; ++i )
	{
		++initSignal[sigLen-1];
		int j = sigLen-1;
		for(; j > 0; --j)
		{
			if( initSignal[j] == 4 )
			{
				initSignal[j] = 0;
				if( ++initSignal[j-1] != 4 )
					break;				
			}
		}
		Str candidateSignal( sigLen, '0' );
		int k = 0;
		for(; k < sigLen; ++k )
		{
			switch( initSignal[k] )
			{
			case 1 :
				candidateSignal[k] = '1';
				break;
			case 2:
				candidateSignal[k] = '2';
				break;
			case 3:
				candidateSignal[k] = '3';
			}
		}
		double count = 0;
		int i = 0;
		for(; i < regions.size(); ++i ){
			if( hasSignal( regions[i].data(), candidateSignal.data()
				, regionLen, mutBp ) )
				++count;
		}
		double w = 1;
		for( i = 0; i < candidateSignal.size(); ++i ){
			w *= background[candidateSignal[i]-'0'];
		}
		double numExp = regions.size() * w *(regionLen-sigLen +1);

		double nonR = (count-numExp)*(count-numExp)/numExp;
		nonrandom += nonR;
		if( count > limitNum ){
			if( count-numExp > 0 )
			out<<setw(20)<<candidateSignal<<setw(20)<<nonR
			<<'\t'<<count<<endl;
			sig.push_back( nonR );
		}
	}
	return nonrandom;
}

Pa_D_D MED_start::getSigScore( const char* seq, int pos, int preSigBP
							  , int regionLen, int ATGPos, Ve_D& distri, M_D& sigWM )
{
	double WMScore = 0;
	int i = 0;
	for(; i < sigWM.getColum(); ++i )
		if( pos + i - preSigBP >= 0 )
			WMScore += sigWM( seq[pos + i - preSigBP]-'0', i );
	return Pa_D_D( WMScore, distri[regionLen - ( ATGPos - pos )] );
}

void MED_start::kerneWordsDistriAndWM( int mutNum )
{
	int i = 0;
	for(; i < SDHitMotifs.size(); ++i )
	{
		Pa_Ve_D_M_D SDResult = signalDistriAndWM
			( SDHitMotifs[i].data()
			, SDRegions, SDRegionsLen
			,SDHitMotifs[i].size(), preSDBP, posSDBP
			, mutNum );
		SDSigWM.push_back( SDResult.second );
		SDSpacerLenDistri.push_back( SDResult.first );
	}
	for( i = 0; i < SDSpacerLenDistri.size(); ++i )
	{
		int iter = 0;
		for(; iter < SDSpacerLenDistri[0].size(); ++iter ){
			SDSpacerLenDistri[i][iter] = log( SDSpacerLenDistri[i][iter] );
	}
		}
	for( i = 0; i < TAHitMotifs.size(); ++i )
	{
		Pa_Ve_D_M_D TAResult = signalDistriAndWM
			( TAHitMotifs[i].data()
			, TARegions, TARegionsLen, TAHitMotifs[i].size()
			, preTABP, posTABP, mutNum );
		TASigWM.push_back( TAResult.second );
		TASpacerLenDistri.push_back( TAResult.first );
	}
	for( i = 0; i < TASpacerLenDistri.size(); ++i )
	{
		int iter = 0;
		for(; iter < TASpacerLenDistri[0].size(); ++iter ){
			TASpacerLenDistri[i][iter] = log( TASpacerLenDistri[i][iter] );
		}
	}
}

Pa_Ve_D_M_D MED_start::signalDistriAndWM( const char* signal
										 , Ve_Str& regions, int regionLen
										 ,int sigLen, int preBP
										 ,int posBP, int mutBp )
{
	int count = 0;
	Ve_D distri( regionLen - sigLen + 1, 0.1 );
	M_D WM( 4, sigLen + preBP + posBP );
 	int i = 0;
 	for(; i < regions.size(); ++i )
	{
		int iter = preBP;
		for(; iter + sigLen + posBP <= regionLen; ++iter )
		{
			if( isSignalInPos( regions[i].data(), signal, iter, mutBp ) )
			{	
				++count;
				char *seq1 = &regions[i][iter - preBP];
				int k = 0;
				for(; k < WM.getColum(); ++k )
				++WM( seq1[k]-'0', k );
				++distri[iter];
			}
		}
	}
	WM.toBeAveraged().toBeLoged();
	int iter = 0;
	for(; iter < distri.size(); ++iter ){
		distri[iter] = ( distri[iter] / count );
	}
	return std::make_pair< Ve_D, M_D >( distri, WM );
}

Ve_Location MED_start::test()
{
	for( Ve_Loca::iterator iter = rvsLocations.begin(); iter != rvsLocations.end(); ++iter )
	{
		const char* seq = iter->isPositive ? positiveSeq.data() : negtiveSeq.data();
		x2LongestORF( seq, iter->location );
		iter->ATGIndexInORF = 0;
		reviseTIS( seq, *iter, XTGIndex );
	}
	return rvsLocations;
}
void MED_start::formPreSTPScore(double aver){
	double count = 0;
	for( Ve_Loca::iterator iter = rvsLocations.begin()
		; iter != rvsLocations.end(); ++iter )
	{
		if( iter->tisScore < aver )
			continue;
		if( fabs(iter->dis2PreSTP) < 50 
			&& iter->ORFLength > learnORFLenLmt
			&& iter->isLongAndNonOverlap )
		{
			++toPreSTPScore[iter->dis2PreSTP];
			++count;
		}
	}
	if( !toPreSTPScore.empty() )
		defaultPreSTPScore = log(1/count);
	else
		defaultPreSTPScore = log(1./rvsLocations.size());
	for( std::map<int,double>::iterator it = toPreSTPScore.begin()
 		; it != toPreSTPScore.end(); ++it ){
		it->second = log(it->second/count);
 	}
}

void  MED_start::prepareForTesting(double aver)
{
	Ve_Location rvsLocations = this->rvsLocations;
	for( Ve_Loca::iterator iter = rvsLocations.begin(); iter != rvsLocations.end(); ++iter )
	{	
		if( iter ->tisScore < aver )
			continue;
		const char* seq = iter->isPositive ? positiveSeq.data() : negtiveSeq.data();
		reviseTIS( seq, *iter );
		if( fabs(iter->location.first - iter->location.second) >= learnORFLenLmt 
			&& iter->isLongAndNonOverlap )
		{
			if( seq[iter->location.first] == '0' )
				XTGIndex[iter->ATGIndexInORF] += pow(1.4, -iter->ATGIndexInORF);
			if( iter->ATGIndexInORF <= 0 )
			{
				switch( seq[iter->location.first] )
				{
				case '0' : XTG[0] += 1; break;
				case '1' : XTG[1] += .8; break;
				case '2' : XTG[2] += 1; break;
				case '3' : XTG[3] += .8; 
				}
			}
		}
	}
	double sum = 0;
	int i =0;
	for(; i < XTGIndex.size(); ++i )
		sum += XTGIndex[i];
	for( i =0; i < XTGIndex.size(); ++i )
		XTGIndex[i] /= sum;
	double all = 0;
	for( i = 0; i < 4; ++i )
		all += XTG[i];

	for( i = 0; i < 4; ++i )
	{
		ATGWM( i, preATGBP ) = log( XTG[i] / all );
		XTG[i] = 0;
	}
}

void MED_start::marskSigBK( Ve_M_D& sigWMs, int beg, int end ){
	if( !sigWMs.empty() ){
		int sigNum = sigWMs.size();
		Ve_M_D bkWMs = Ve_M_D( sigNum, M_D(4,sigWMs[0].getColum(),0.001));
		int i = 0;
		for(; i < rvsLocations.size(); ++i ){
			if( rvsLocations[i].ORFLength > learnORFLenLmt 
				&& rvsLocations[i].isLongAndNonOverlap )
			{
				const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() : negtiveSeq.data() ;
				int ATGPos = getPrePhaseATG(seq,rvsLocations[i].location.first);
				if( ATGPos < rvsLocations[i].location.first && ATGPos +beg >= 0 ){
					int j = 0;
					for(; j < sigNum; ++j ){
						int pos = ATGPos + beg;
						int endP = ATGPos + end;
						double score = -1000000;
						int posIndex = pos;
						for( ; pos <= endP; ++pos ){
							if( isSignalInPos( seq, SDHitMotifs[j].data(), pos, 1 ) )
							{
								double sigScore = 0;
								int k = 0;
								for(; k < sigWMs[j].getColum(); ++k ){
									if( pos >= 0 )
										sigScore += sigWMs[j]( seq[pos -preSDBP+ k ]-'0', k );
								}
								if( sigScore > score ){
									score = sigScore;
									posIndex = pos;
								}
							}
						}
						int iter = 0;
						for(; iter < sigWMs[j].getColum(); ++iter )
							++bkWMs[j]( seq[posIndex-preSDBP+iter]-'0', iter );
					}
				}
			}
		}
		int n = 0;
		for(; n < sigNum; ++n ){
			sigWMs[n].toBeExp();
			bkWMs[n].toBeAveraged();
			int p = 0;
			for(; p < sigWMs[n].getLine(); ++p ){
				int q = 0;
				for(; q < sigWMs[n].getColum(); ++q ){
					sigWMs[n](p,q) /= bkWMs[n](p,q);
				}
			}
			sigWMs[n].toBeAveraged();
			sigWMs[n].toBeLoged();
		}
	}
}

void MED_start::formATGWM(double aver)
{
	if( !isGCRich ){
		ATGWM = M_D( 4, ATGWM.getColum(), 0.0001 );
		int count = 0;
		int i = 0;
		for(; i < rvsLocations.size(); ++i )
		{
			if( rvsLocations[i].tisScore < aver )
				continue;
			++count;
			const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() : negtiveSeq.data() ;
			if( rvsLocations[i].ORFLength > learnORFLenLmt && 
				rvsLocations[i].location.first -preATGBP >= 0 
				&& rvsLocations[i].isLongAndNonOverlap )
			{
				int iter = 0;
				for(; iter < ATGWM.getColum(); ++iter )
					++ATGWM( seq[rvsLocations[i].location.first + iter-preATGBP] - '0', iter );
			}
		}
		ATGWM.toBeAveraged().toBeLoged();
	}
	else{
		ATGWM = M_D( 4, ATGWM.getColum(), 0.0001 );
		M_D ATGWMBK( 4, ATGWM.getColum(), 0.0001 );

		int i = 0;
		for(; i < rvsLocations.size(); ++i ){
			const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() : negtiveSeq.data() ;
			if(	rvsLocations[i].location.first -preATGBP >= 0 ){
				int iter = 0;
				for(; iter < ATGWM.getColum(); ++iter )
					++ATGWM( seq[rvsLocations[i].location.first + iter-preATGBP] - '0', iter );
			}
			int ATGPos = getNextPhaseTIS(seq,rvsLocations[i].location.first+3);
			if( ATGPos > rvsLocations[i].location.first && ATGPos-preATGBP >= 0){
				int iter = 0;
				for(; iter < ATGWM.getColum(); ++iter )
					++ATGWMBK( seq[ATGPos + iter-preATGBP] - '0', iter );
			}
		}
		ATGWM.toBeAveraged();
		ATGWMBK.toBeAveraged();
		for( i = 0; i < ATGWM.getLine(); ++i ){
			int j = 0;
			for( ; j < preATGBP; ++j )
				ATGWM( i, j ) /= background[i];
			for( ; j < ATGWM.getColum(); ++j ){
				ATGWM( i, j ) /= ATGWMBK(i,j);
			}
		}
		ATGWM( 1, preATGBP ) = 0.0001;
		ATGWM.toBeAveraged();
		for( i = 0; i < 4; ++i ){
			if( ATGWM( i, preATGBP ) < 0.001 && i != 1 ){
				cout<<"Extrem low percentage for certain TIS.\n"
					"Will omit the start codon column in the weight\n"
					"matrix of starts"<<endl;
				int j = 0;
				for(; j < 4; ++j )
				{
					ATGWM( j, preATGBP ) = 0.25;
				}
				break;
			}
		}
		ATGWM.toBeLoged();
	}
}

void MED_start::findHitMotifs( int mutNum, Ve_Str& regions
							  , Ve_Str& candidateMotifs
							  , Ve_Str& hitMotifs, int regionLen
							  , Str file
							  , bool (*filter)(Str&) )
{
	ifstream in( file.data() );
	if( !in.good() )
		return;
	std::map<int,std::string,std::greater<int> > sig;
	while( !in.eof() )
	{
		Str signal;
		double times;
		double t;
		in>>signal>>times>>t;
		if( !signal.empty() && t > 50 ){
			sig[int(times*100)] = signal;
		}
		Str _st_;
		std::getline( in, _st_ );
	};
	in.close();
	if( sig.empty() ){
	}
	system( Str("rm \"" + file +'\"' ).data() );
	std::map<int,std::string,std::greater<int> >::iterator iter = sig.begin();
	std::map<double,std::pair<Ve_D,Str>,std::greater<double> > kernel;
	int i = 0;
	for(; i < 5 && iter != sig.end(); ++iter, ++i )
	{
		Ve_D tmp = signalDistri( iter->second.data()
			, mutNum, regionLen, regions );
 		kernel[GetDerivate( tmp )] = std::make_pair<Ve_D,Str>( tmp, iter->second );
		candidateMotifs.push_back( iter->second );
	}
	std::map<double,std::pair<Ve_D,Str>,std::greater<double> >::iterator kernelIter = kernel.begin();
	int j = 0;
	for(; j < 3; ++j ) 
	{
		if( kernelIter == kernel.end() )
			break;
		std::map<double,std::pair<Ve_D,Str>,std::greater<double> >
			::iterator tmp = kernelIter;
		hitMotifs.push_back( tmp->second.second );
		++kernelIter;
		if( fabs( tmp->first - (kernelIter)->first ) > 0.01 )
			break;
	}
}

Ve_D MED_start::signalDistri(  Str signal , int mutNum,  int regionLen, Ve_Str& regions )
{
	int disCount = 0;
	Ve_D distri( regionLen - sigLen + 1, 0 );
 	int i = 0;
 	for(; i < regions.size(); ++i ){
		int iter = 0;
		for(; iter + sigLen <= regionLen; ++iter )
		{
			if( isSignalInPos( regions[i].data(), signal.data(), iter, mutNum ) ){	
				++distri[iter];
				++disCount;
			}
		}
	}
	int iter = 0;
	for(; iter < distri.size(); ++iter )
		distri[iter] = distri[iter] / disCount;
	return distri;
}

double MED_start::reviseTIS( const char* seq, Location_T& iter, Ve_D& XTGIndexTmp )
{
	
	int STP = iter.location.second;
	int ATGPosition = iter.location.first;
	int count = iter.ATGIndexInORF;
	double score = -100000;
	int predictATGPos = ATGPosition;
	iter.SDPos = 0;
	iter.TAPos = 0;
	Location_T tmpl = iter;
	tmpl.location2FileForm( seqLen );
	std::vector<std::pair<std::pair<Location_T, bool>,double> > list;
	int dis2PreSTP = iter.dis2PreSTP;
	do
	{
		int ATGtmp = ATGPosition;
		ATGPosition = getNextPhaseTIS( seq, ATGPosition );
		dis2PreSTP += ATGPosition - ATGtmp;
		if( ATGPosition  == -1 || ATGPosition >=  STP 
			|| count >= XTGIndex.size() || STP - ATGPosition <= 87 || ATGPosition<35)
			break;
		
		double SDScore = -100000;
		int SDPos = 0;
		if( SDSpacerLenDistri.size() > 0 ){
			int pos = ATGPosition - SDRegionsLen;
			for(; pos <= ATGPosition - sigLen; ++pos )
			{
				int i = 0;
				for(; i < SDSpacerLenDistri.size(); ++i )
				{
					Pa_D_D rbsScore = getSigScore( seq
						, pos, preSDBP, SDRegionsLen, ATGPosition
						, SDSpacerLenDistri[i], SDSigWM[i] );
					double SDS = (rbsScore.first 
						+ rbsScore.second) * SDWeight;
					if( SDS > SDScore ){
						SDScore = SDS; 
						SDPos = pos;
					}
				}
			}
		}
		double TAScore = -100000;
		int TAPos = 0;
		if( TASpacerLenDistri.size() > 0 ){
			for( int pos = ATGPosition - TARegionEnd - TARegionsLen
				; pos <= ATGPosition - TARegionEnd - sigLen; ++pos )
			{
				int i = 0;
				for(; i < TASpacerLenDistri.size(); ++i )
				{
					Pa_D_D taScore = getSigScore( seq, pos, preTABP
						,TARegionsLen + TARegionEnd, ATGPosition
						,TASpacerLenDistri[i], TASigWM[i] );
					double TAS = (taScore.first + taScore.second)
						* TAWeight;
					if( TAS > TAScore ){
						TAScore = TAS;
						TAPos = pos;
					}
				}
			}
		}
		double tmp = (SDScore > -100000 ? SDScore : 0)
			+ (TAScore > -100000 ? TAScore : 0);
		double sigScore = tmp;
		double ATGS = 0;
		int i = 0;
		for(; i < ATGWM.getColum(); ++i )
		{
			if( ATGPosition + i - preATGBP > 0 )
				ATGS += ATGWM( seq[ ATGPosition + i - preATGBP] - '0', i );
		}

		double ATGScore = ATGWeight*ATGS ;
		tmp += ATGScore + ATGWeight*log(XTGIndexTmp[count]);
		double DisScore = 0;
		if( !toPreSTPScore.empty() ){
			std::map<int, double>::iterator it = toPreSTPScore.begin();
			int dis = dis2PreSTP ;
			for( ; it != toPreSTPScore.end(); ++it ){
				if( it->first >= dis )
					break;
			}
			if( it == toPreSTPScore.end() || it == toPreSTPScore.begin()){
				DisScore =  ATGWeight*defaultPreSTPScore;
			}
			else{
				std::map<int, double>::iterator posIt = it;
				DisScore =  ATGWeight*(it->second + posIt->second)/2;
			}
		}
		tmp += DisScore;
		
		if( isGCRich ){	
			double codingS = 0;
			if( ATGPosition > 90 ){
				double noncodingP = 1-getCodProd( int( getGB(seq, ATGPosition - 90, ATGPosition) ) );
				double codingP = getCodProd( int( getGB(seq, ATGPosition, ATGPosition+90) ) );
				codingS = ATGWeight*(log(noncodingP)+log(codingP));	
			}
			tmp += codingS;
		}
		if( tmp > score )
		{
			score = tmp;
			predictATGPos = ATGPosition;
			iter.ATGIndexInORF = count;
			iter.TAPos = TAPos;
			iter.SDPos = SDPos;
			iter.sigScore = sigScore;
			iter.dis2PreSTP = dis2PreSTP ;
			iter.ATGScore = ATGScore;
			iter.SDScore = SDScore;
			iter.TAScore = TAScore;
			iter.disScore = DisScore;
		}
	    Location_T tmpl = iter;
		tmpl.location.first = ATGPosition;
		list.push_back
			( std::pair<std::pair<Location_T,bool>,double >
			(std::pair<Location_T,bool>(tmpl,false), tmp)
			);
		ATGPosition += 3;
		++count;
	}while( true );

	if( predictATGPos > 50 ){
		Str preSeq;
		int i = 0;
		for(; i < 50; ++i ){
			preSeq += seq[predictATGPos-50+i];
		}
		iter.preSeq = preSeq;
	}
	iter.location.first = predictATGPos;
	iter.ORFLength = STP - predictATGPos;
	iter.tisScore = score;
	return score;
}

double MED_start::reviseTIS( const char* seq, Location_T& iter )
{
	int STP = iter.location.second;
	int ATGPosition = iter.location.first;
	int count = iter.ATGIndexInORF;
	double score = -100000;
	int predictATGPos = ATGPosition;
	
	int dis2PreSTP = iter.dis2PreSTP;
	do
	{
		int ATGtmp = ATGPosition;
		ATGPosition = getNextPhaseTIS( seq, ATGPosition );
		dis2PreSTP += ATGPosition - ATGtmp;
		if( ATGPosition  == -1 || ATGPosition >=  STP 
			|| count >= XTGIndex.size() || STP - ATGPosition <= 87 )
			break;

		double SDScore = -100000;
		if( SDSpacerLenDistri.size() > 0 ){
			int pos = ATGPosition - SDRegionsLen;
			for(; pos <= ATGPosition - sigLen; ++pos )
			{
				int i = 0;
				for(; i < SDSpacerLenDistri.size(); ++i )
				{
					Pa_D_D rbsScore = getSigScore( seq
						, pos, preSDBP, SDRegionsLen, ATGPosition
						, SDSpacerLenDistri[i], SDSigWM[i] );
					double SDS = (rbsScore.first + rbsScore.second)
						*SDWeight;
					if( SDS > SDScore )
						SDScore = SDS;
				}
			}
		}
		double TAScore = -100000;
		if( TASpacerLenDistri.size() > 0 ){
			for(int pos = ATGPosition - TARegionEnd - TARegionsLen
				; pos <= ATGPosition - TARegionEnd - sigLen; ++pos )
			{
				int i = 0;
				for(; i < TASpacerLenDistri.size(); ++i )
				{
					Pa_D_D taScore = getSigScore( seq, pos, preTABP
						,TARegionsLen + TARegionEnd, ATGPosition
						,TASpacerLenDistri[i], TASigWM[i] );
					double TAS = (taScore.first + taScore.second)
						*TAWeight; 
					if( TAS > TAScore )
						TAScore = TAS;
				}
			}
		}
		double tmp = (SDScore > -100000 ? SDScore : 0)
			+ (TAScore > -100000 ? TAScore : 0);
		double ATGScore = 0;
		int i = 0;
		for(; i < ATGWM.getColum(); ++i )
		{
			if( ATGPosition + i - preATGBP > 0 )
				ATGScore += ATGWM( seq[ ATGPosition + i - preATGBP] - '0', i );
		}

		tmp += ATGWeight*ATGScore;
		
		if(isGCRich){
			double codingS = 0;
			if( ATGPosition > 90 ){
				double noncodingP = 1-getCodProd( int( getGB(seq, ATGPosition - 90, ATGPosition) ) );
				double codingP = getCodProd( int( getGB(seq, ATGPosition, ATGPosition+90) ) );
				codingS = ATGWeight*(log(noncodingP)+log(codingP));	
			}
			tmp += codingS;
		}

		if( tmp > score )
		{
			score = tmp;
			predictATGPos = ATGPosition;
			iter.ATGIndexInORF = count;
			iter.dis2PreSTP = dis2PreSTP ;
		}			
		ATGPosition += 3;
		++count;
	}while( true );
	iter.location.first = predictATGPos;
	iter.ORFLength = STP - predictATGPos;

	return score;
}

Ve_Location MED_start::doIteration()
{
	kerneWordsDistriAndWM(mutBp);
	if(isGCRich){
		marskSigBK( SDSigWM, -SDRegionsLen );
		formGB();
	}
	Ve_Location rvsLocations = this->rvsLocations;
	double averTIS = 0;
	for( Ve_Loca::iterator iter = rvsLocations.begin(); iter != rvsLocations.end(); ++iter ){
		averTIS += iter->tisScore;
	}
	averTIS = -1000000;
	formATGWM(averTIS);
	formPreSTPScore(averTIS);
	prepareForTesting(averTIS);

	return test();
}


bool MED_start::hasSignal( const char* seq, const char* signal, int SDRegionsLen, int mutBp )
{
	int iter = 0;
	for(; iter + sigLen <= SDRegionsLen - posSDBP; ++iter )
		if( isSignalInPos( seq, signal, iter, mutBp ) )
			return true;
	return false;
}

bool MED_start::isSignalInPos( const char* seq, const char* signal, int pos,int mutBp)
{
	int mutNum = 0;
	int i = pos;
	for(; i < pos + sigLen; ++i )
	{
		if( seq[i] != signal[i-pos] )
			++mutNum;
	}
	return mutNum <= mutBp ? true : false;
}

void MED_start::setHitMotifs(Ve_Str SDHitMotifs1)
{
	this->SDHitMotifs = SDHitMotifs1;
}

void MED_start::setCandidateMotifs(Ve_Str SDCandidateMotifs)
{
	this->SDCandidateMotifs = SDCandidateMotifs;
}

void MED_start::setParameters( CParametersDlg* m_parameters )
{
	learnORFLenLmt = m_parameters->m_ORFLen;
	SDRegionsLen = m_parameters->m_SDRegion;
	TARegionsLen = m_parameters->m_TARegion;
	sigLen = m_parameters->m_sigLen;
	mutBp = m_parameters->m_sigMutNum;
	preSDBP = m_parameters->m_SDPreNum;
	posSDBP = m_parameters->m_SDPosNum;
	preTABP = m_parameters->m_TAPreNum;
	posTABP = m_parameters->m_TAPosNum;
	preATGBP =  m_parameters->m_atgPreNum;
	posATGBP = m_parameters->m_atgPosNum;
	TARegionEnd = m_parameters->m_TARegionEnd;
	toPreSTPScore.clear();
}

void MED_start::resetLearnedRst(bool)
{
	ATGWM = M_D( 4, preATGBP + posATGBP + 3,0.0001 );
	XTG = Ve_D( 4, 0.01 );
	XTGIndex = Ve_D( XTGIndex.size(), 1 );
	SDRegions = Ve_Str();
	TARegions = Ve_Str();
	SDHitMotifs = Ve_Str();
	TAHitMotifs = Ve_Str();
	SDCandidateMotifs= Ve_Str();
	TACandidateMotifs = Ve_Str();
	SDSigWM.erase( SDSigWM.begin(), SDSigWM.end() );
	TASigWM.erase( TASigWM.begin(), TASigWM.end() );
	SDSpacerLenDistri.erase( SDSpacerLenDistri.begin(), SDSpacerLenDistri.end() );
	TASpacerLenDistri.erase( TASpacerLenDistri.begin(), TASpacerLenDistri.end() );
	toPreSTPScore.clear();
}

void MED_start::setRegions()
{
	int i = 0;
	for(; i < rvsLocations.size(); ++i )
	{
		if( rvsLocations[i].ORFLength >= learnORFLenLmt 
			&& rvsLocations[i].isLongAndNonOverlap)
		{
			Location_T tmp = rvsLocations[i];
			Str& seq = rvsLocations[i].isPositive ? positiveSeq : negtiveSeq;
			if( rvsLocations[i].location.first - SDRegionsLen >= 0 )
			{	
				Str temp = seq.substr( rvsLocations[i].location.first - SDRegionsLen,
					SDRegionsLen );
				if( temp.size() == SDRegionsLen )
					SDRegions.push_back( temp );
			}
			//TATA
			if( rvsLocations[i].location.first 
				- TARegionEnd - TARegionsLen>= 0 )
			{	
				Str temp = seq.substr( rvsLocations[i].location.first 
					- TARegionEnd - TARegionsLen,
					TARegionsLen );
				if( temp.size() == TARegionsLen )
					TARegions.push_back( temp );
			}
			//ATG
			if( rvsLocations[i].location.first 
				- preATGBP >= 0 )
			{	
				Str temp = seq.substr( rvsLocations[i].location.first - preATGBP,
					preATGBP + posATGBP + 3 );
				if( temp.size() == preATGBP + posATGBP + 3 )
					ATGRegions.push_back( temp );
			}
		}
	}
}

double MED_start::getRvsLocations(int iterN )
{
	int i = 0;
	for(; i < rvsLocations.size(); ++i ){
		rvsLocations[i].fileForm2Location( seqLen );
	}

	int iterNum = 0;

	double lemda = 0;
	for( double iterationSTPTag = 1000000; ; )
	{
		cout<<".";
		setRegions();
		SDWeight = findMotifs( 0, SDRegions, SDRegionsLen
			, m_motifsPath+"LMotifs.txt"
			, (SDRegionsLen-sigLen+1)*SDRegions.size()*0.00) 
			/ SDRegions.size()/(SDRegions[0].size()-sigLen+1);
		TAWeight = findMotifs( 0, TARegions, TARegionsLen
			, m_motifsPath+"UMotifs.txt"
			, (TARegionsLen-sigLen+1)*TARegions.size()*0.00)
			/ TARegions.size()/(TARegions[0].size()-sigLen+1);

		lemda = SDWeight - TAWeight;
		double sum = ( SDWeight + TAWeight);
		SDWeight /= sum;
		TAWeight /= sum;
		ATGWeight = (SDWeight>TAWeight)?SDWeight:TAWeight;
		findHitMotifs(mutBp, SDRegions
			, SDCandidateMotifs, SDHitMotifs
			,SDRegionsLen, m_motifsPath+"LMotifs.txt", isSDSignal );
		findHitMotifs(mutBp, TARegions
			, TACandidateMotifs, TAHitMotifs
			, TARegionsLen, m_motifsPath+"UMotifs.txt", isTASignal);

		rvsLocations = doIteration();

		double percent = fabs( iterationSTPTag /  getSTPTag() );
		if( ( percent > 0.99 && percent <1.01) || ++iterNum 
			== iterN ){
			parameters2File(parafile);
			break;
		}
		else
			iterationSTPTag = getSTPTag();
		resetLearnedRst(true);
	}
	cout<<endl;
	return 	lemda;
}

void MED_start::parameters2File(Str fileName)
{
	if( fileName.empty() )
		return;
	std::ofstream out( fileName.data() );
	out<<"default"<<endl;
	out<<learnORFLenLmt<<'\t'<<TARegionsLen<<'\t'<<TARegionEnd<<'\t'<<SDRegionsLen
		<<'\t'<<sigLen<<'\t'<<mutBp<<endl;
	out<<"GC_content"<<endl;
	out<<background[0]<<'\t'<<background[1]<<'\t'
		<<background[2]<<'\t'<<background[3]<<endl;
	out<<"SD_signal"<<endl;
	out<<preSDBP<<'\t'<<posSDBP<<endl;
	out<<SDHitMotifs.size()<<endl;
	int i = 0;
	for(; i < SDHitMotifs.size(); ++i ){
		out<<SequenceTransform_T::digital2CharSeq(SDHitMotifs[i])<<endl;
		(SDSigWM[i].toBeExp()).matrixOut(out);
		out<<SDSpacerLenDistri[i].size()<<endl;
		int j = 0;
		for(; j < SDSpacerLenDistri[i].size(); ++j ){
			out<<exp(SDSpacerLenDistri[i][j])<<endl;
		}
	}
	out<<"TA_signal"<<endl;
	out<<preTABP<<'\t'<<posTABP<<endl;
	out<<TAHitMotifs.size()<<endl;
	for( i = 0; i < TAHitMotifs.size(); ++i ){
		out<<SequenceTransform_T::digital2CharSeq(TAHitMotifs[i])<<endl;
		(TASigWM[i].toBeExp()).matrixOut(out);
		out<<TASpacerLenDistri[i].size()<<endl;
		int j = 0;
		for(; j < TASpacerLenDistri[i].size(); ++j ){
			out<<exp(TASpacerLenDistri[i][j])<<endl;
		}
	}
	out<<"XTG"<<endl;
	out<<preATGBP<<'\t'<<posATGBP<<endl;
	(ATGWM.toBeExp()).matrixOut(out);
	out<<XTGIndex.size()<<endl;
	int j = 0;
	for(; j < XTGIndex.size(); ++j ){
		out<<(XTGIndex[j])<<endl;
	}
	out<<"predis"<<endl;
	out<<toPreSTPScore.size()<<endl;
	std::map<int,double>::iterator it = toPreSTPScore.begin();
	for( ; it != toPreSTPScore.end(); ++it ){
		out<<setw(10)<<it->first<<setw(15)<<it->second<<endl;
	}
	out<<defaultPreSTPScore<<endl;
	out<<"weight"<<endl;
	out<<SDWeight<<'\t'<<TAWeight<<'\t'<<ATGWeight<<endl;
	out<<"TGB "<<CodGB.distri.size()<<endl;
	std::map<int,double>::iterator it1;
	for( it1 = CodGB.distri.begin()
		; it1 != CodGB.distri.end(); ++it1 ){
		out<<it1->first<<"  "<<it1->second<<endl;
	}
	out<<"FGB "<<NonCodGB.distri.size()<<endl;
	for( it1 = NonCodGB.distri.begin()
		; it1 != NonCodGB.distri.end(); ++it1 ){
		out<<it1->first<<"  "<<it1->second<<endl;
	}
}

void MED_start::readParametersFromFile( Str fileName ){
	std::ifstream in( fileName.data() );
	while( !in.eof() ){
		Str hint;
		in>>hint;
		if( hint == "default" ){
			in>>learnORFLenLmt>>TARegionsLen>>TARegionEnd>>SDRegionsLen
				>>sigLen>>mutBp;
		}
		else if( hint == "SD_signal" ){
			in>>preSDBP>>posSDBP;
			int num = 0;
			in>>num;
			int i = 0;
			for(; i < num; ++i ){
				Str sig;
				in>>sig;
				SDHitMotifs.push_back( SequenceTransform_T::
					char2digitalSeq( sig ) );
				M_D matrix;
				matrix.matrixIn(in);
				SDSigWM.push_back( matrix.toBeLoged() );
				int size;
				in>>size;
				Ve_D sp( size, 0 );
				int j = 0;
				for(; j < size; ++j ){
					in>>sp[j];
					sp[j] = log( sp[j] );
				}
				SDSpacerLenDistri.push_back( sp );
			}
		}
		else if( hint == "TA_signal" ){
			in>>preTABP>>posTABP;
			int num = 0;
			in>>num;
			int i = 0;
			for(; i < num; ++i ){
				Str sig;
				in>>sig;
				TAHitMotifs.push_back( SequenceTransform_T::
					char2digitalSeq( sig ) );
				M_D matrix;
				matrix.matrixIn(in);
				TASigWM.push_back( matrix.toBeLoged() );
				int size;
				in>>size;
				Ve_D sp( size, 0 );
				int j = 0;
				for(; j < size; ++j ){
					in>>sp[j];
					sp[j] = log( sp[j] );
				}
				TASpacerLenDistri.push_back( sp );
			}
		}
		else if( hint == "XTG" ){
			in>>preATGBP>>posATGBP;
			M_D matrix;
			matrix.matrixIn(in);
			ATGWM = matrix.toBeLoged();		
			int size;
			in>>size;
			XTGIndex = Ve_D( size, 0 );
			int j = 0;
			for(; j < size; ++j ){
				in>>XTGIndex[j];
			}
		}
		else if( hint == "GC_content" ){
			Ve_D bg(4,0);
			in>>bg[0]>>bg[1]>>bg[2]>>bg[3];
			background= bg;
		}
		else if( hint == "weight" ){
			in>>SDWeight>>TAWeight>>ATGWeight;
		}
		else if( hint == "predis" ){
			int size;
			in>>size;
			int i = 0;
			for(; i < size; ++i ){
				int dis;
				double s;
				in>>dis>>s;
				toPreSTPScore[dis] = s;
			}
			in>>defaultPreSTPScore;
		}
		else if( hint == "TGB" ){
			int num = 0;
			in>>num;
			int i = 0;
			for(; i < num; ++i ){
				int a;
				double b;
				in>>a>>b;
				CodGB.distri.insert(std::make_pair(a,b));
			}
		}
		else if( hint == "FGB" ){
			int num = 0;
			in>>num;
			int i = 0;
			for(; i < num; ++i ){
				int a;
				double b;
				in>>a>>b;
				NonCodGB.distri.insert(std::make_pair(a,b));
			}
		}
		else{
			if( !hint.empty() ){
				cout<<"ERR while reading parameters!"<<endl;
				cout<<"ERR line:\t"<<hint<<endl;
				exit( 0 );
			}
			return;
		}
	}
	in.close();
}//*/
#include <sstream>
void outMatrix(std::ostream & out,const M_D& matrix)
{
	out<<std::setw(15)<<'A'<<std::setw(15)<<'C'
		<<std::setw(15)<<'G'<<std::setw(15)<<'T'<<endl;	
	int col=0;
	for(;col<matrix.getColum();++col)
	{
		int lin=0;
		for(;lin<matrix.getLine();++lin){
			if( fabs(matrix(lin,col)) > 0.001 ){
				std::stringstream ost;
				ost <<matrix(lin,col);
				Str str = ost.str() + ".000";
				out<<setw(15)<<str.substr(0, str.find_first_of(".")+4);
			}
			else
				out<<setw(15)<<'0';
		}
		out<<std::endl;
	}
}

void MED_start::outReadableParameters( Str file ,  Str Class ){
	ofstream out( file.data(), std::ios::app );
	out<<Class<<" TIS model parameters"<<endl;
	out<<"Minimal length of ORF: "<<learnORFLenLmt<<endl;
	out<<"GC%: "<<(background[1]+background[2])<<endl;
	out<<"Positional weight Matrix for signals in lower"
		" upstream region [-20, 0)."<<endl;
	int i = 0;
	for(; i < SDHitMotifs.size(); ++i ){
		out<<"Signal: "<<SequenceTransform_T::digital2CharSeq(SDHitMotifs[i])<<endl;
		outMatrix(out,SDSigWM[i].toBeExp());
	}
	out<<"Positional weight Matrix for signals in upper"
		" upstream region [-35, -15)."<<endl;
	for( i = 0; i < TAHitMotifs.size(); ++i ){
		out<<"Signal: "<<SequenceTransform_T::digital2CharSeq(TAHitMotifs[i])<<endl;
		outMatrix(out,TASigWM[i].toBeExp());
	}
	out<<"Positional weight Matrix for TIS."<<endl;
	outMatrix(out,ATGWM.toBeExp());
	out<<"Weight for the ith candidate start codon"<<endl;
	out<<setw(10)<<"ith"<<setw(15)<<"weight"<<endl;
	int j = 0;
	for(; j < XTGIndex.size(); ++j ){
		std::stringstream ost;
		ost <<(XTGIndex[j]);
		Str str = ost.str() + ".000";
		out<<setw(10)<<j+1<<setw(15)
			<<str.substr(0, str.find_first_of(".")+4)<<endl;
	}
	out<<"Distribution of distance from a TIS to its immediately "
		"upstream stop codon"<<endl;
	std::map<int,double>::iterator it = toPreSTPScore.begin();
	out<<setw(10)<<"distance"<<setw(15)<<"%"<<endl;
	for( ; it != toPreSTPScore.end(); ++it ){
		std::stringstream ost;
		ost <<exp(it->second);
		Str str = ost.str() + ".000";
		out<<setw(10)<<it->first<<setw(15)
			<<str.substr(0, str.find_first_of(".")+4)<<endl;
	}
	std::stringstream ost;
	ost <<exp(defaultPreSTPScore);
	Str str = ost.str() + ".000";
	out<<setw(10)<<"others"<<setw(15)
		<<str.substr(0, str.find_first_of(".")+4)<<endl;

	if( isGCRich ){
		out<<"GB distribution for coding region downstream of TIS"<<endl;
		out<<setw(5)<<"GB"<<setw(15)<<"Prob."<<endl;
		std::map<int,double>::iterator it1;
		for(it1 = CodGB.distri.begin()
			; it1 != CodGB.distri.end(); ++it1 ){
			out<<setw(5)<<it1->first;
			std::stringstream ost;
			ost<<it1->second;
			Str str = ost.str() + ".000";
			out<<setw(15)<<str.substr(0, str.find_first_of(".")+4)<<endl;
		}
		out<<"GB distribution for non-coding region upstream of TIS"<<endl;
		out<<setw(5)<<"GB"<<setw(15)<<"Prob."<<endl;
		for( it1 = NonCodGB.distri.begin()
			; it1 != NonCodGB.distri.end(); ++it1 ){
			out<<setw(5)<<it1->first;
				std::stringstream ost;
			ost<<it1->second;
			Str str = ost.str() + ".000";
			out<<setw(15)<<str.substr(0, str.find_first_of(".")+4)<<endl;
		}	
	}
	out<<"Weight for motif in region [-35,-15) and [-20, 0 )"<<endl;
	out<<setw(15)<<TAWeight<<setw(15)<<SDWeight<<endl;
	out.close();
}

void MED_start::reviseGivenORFs(Ve_Location &locations){
	int i = 0;
	for(; i < locations.size(); ++i ){
		locations[i].fileForm2Location( seqLen );

	}
	for( i = 0; i < locations.size(); ++i ){
		const char* seq1 = locations[i].isPositive 
			? positiveSeq.data() : negtiveSeq.data();
		x2LongestORF( seq1, locations[i].location );
		locations[i].ATGIndexInORF = 0;
		reviseTIS( seq1, locations[i], XTGIndex );
		int code = seq1[locations[i].location.first]-'0';
		locations[i].TISCode = code == 0 ? "ATG" : (code == 2 
			? "GTG" : "TTG" );

		const Str& seq2 = locations[i].isPositive ? positiveSeq : negtiveSeq;
		Str _st_=seq2.substr(locations[i].SDPos, sigLen);
		locations[i].SDSig = SequenceTransform_T::digital2CharSeq(
			_st_);
		_st_=seq2.substr(locations[i].TAPos, sigLen);
		locations[i].TASig = SequenceTransform_T::digital2CharSeq(
			_st_);
		if( locations[i].location.first > 35 ){
			locations[i].Sig = SDWeight > TAWeight ? locations[i].SDSig 
				: locations[i].TASig;
			locations[i].SigPos = SDWeight > TAWeight ? locations[i].SDPos 
				- locations[i].location.first + sigLen: locations[i].TAPos 
				- locations[i].location.first + sigLen; 
		}else{
			locations[i].Sig ="-";
			locations[i].SigPos = 10000;
		}
		locations[i].location2FileForm( seqLen );
	}
}

void MED_start::formGB(){
	CodGB = GBDistri();
	NonCodGB = GBDistri();
	double NonCodGBnumber = 0, CodGBnumber = 0;
	int i = 0;
	for(; i < rvsLocations.size(); ++i ){
		const Location_T& rvsLocation = rvsLocations[i];
		const char* seq = rvsLocation.isPositive ? positiveSeq.data()
			: negtiveSeq.data();
		if(rvsLocation.location.first>90
			&& rvsLocations[i].isLongAndNonOverlap ){
			double gb = getGB( seq, rvsLocation.location.first - 90, rvsLocation.location.first );
			++NonCodGB.distri[ gb];
			++NonCodGBnumber;
		}
		double gb = getGB( seq, rvsLocation.location.first,  rvsLocation.location.first + 90 );
		++CodGB.distri[ gb ];
		++CodGBnumber;
	}
	std::map<int,double>::iterator it;
	for(it = NonCodGB.distri.begin()
		; it != NonCodGB.distri.end(); ++it ){
		it->second /= NonCodGBnumber;
	}
	for(  it = CodGB.distri.begin(); it != CodGB.distri.end(); ++it ){
		it->second /= CodGBnumber;
	}
}
