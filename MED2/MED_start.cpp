#include "MED_start.h"
#include "SequenceTransform.h"
#include "OftenUsedOperatLib.h"

#pragma warning (disable:4786)

MED_start::MED_start
( GeneSeq& geneseq, Ve_Location& rvsLocations, CParametersDlg parameters )
	: GeneInfo_T( geneseq ), XTGIndex( 6, 1 ), XTG( 4, 0.01 )
{
	setRvsLocation( rvsLocations );
	setParameters( &parameters );
	ATGWM = M_D( 4, preATGBP + posATGBP + 3,0.0001 );
}

void MED_start::findMotifs(int mutBp)
{
	Ve_I initSignal( sigLen, 0 );
	initSignal[sigLen-1] = -1;
	
	int upBoundary = pow( 4, sigLen );
	ofstream out( (m_motifsPath+"motifs.txt").data() );
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
		for(; i < motifRegions.size(); ++i ){
			if( hasSignal( motifRegions[i].data(), candidateSignal.data()
				, motifRegionsLen, mutBp ) )
				++count;
		}
		double w = 1;
		for( i = 0; i < candidateSignal.size(); ++i ){
			w *= background[candidateSignal[i]-'0'];
		}
		double numExp = motifRegions.size() * w *(motifRegionsLen-sigLen +1);

		double nonR = (count-numExp)*(count-numExp)/numExp;
		if( nonR > 0.001 )
			out<<setw(20)<<candidateSignal<<setw(20)<<nonR<<endl;
	}
	out.close();
}

Pa_D_D MED_start::getRBSScore( const char* seq, int pos, int ATGPos, Ve_D& distri, M_D& sigWM )
{
	double WMScore = 0;
	int i = 0;
	for(; i < sigWM.getColum(); ++i )
		if( pos + i - preSigBP >= 0 )
			WMScore += sigWM( seq[pos + i - preSigBP]-'0', i );
	return Pa_D_D( WMScore, distri[motifRegionsLen - ( ATGPos - pos )] );
}

void MED_start::kerneWordsDistriAndWM( int mutNum )
{
	int i = 0;
	for(; i < hitMotifs.size(); ++i )
	{
		Pa_Ve_D_M_D result = signalDistriAndWM( hitMotifs[i].data(), hitMotifs[i].size(), mutNum );
		SigWM.push_back( result.second );
		spacerLenDistri.push_back( result.first );
	}
	if(spacerLenDistri.size() < 1){
		printf("No root ORFs found\n");
		exit (EXIT_FAILURE);
	}
	int iter = 0;
	for(; iter < spacerLenDistri[0].size(); ++iter )
	{
		int i = 0;
		for(; i < spacerLenDistri.size(); ++i )
			spacerLenDistri[i][iter] = log( spacerLenDistri[i][iter] );
	}
}

Pa_Ve_D_M_D MED_start::signalDistriAndWM( const char* signal, int sigLen, int mutBp )
{
	int count = 0;
	Ve_D distri( motifRegionsLen - sigLen + 1, 0.1 );
	M_D WM( 4, sigLen + preSigBP + posSigBP );
 	int i = 0;
 	for(; i < motifRegions.size(); ++i )
	{
		int iter = preSigBP;
		for(; iter + sigLen + posSigBP <= motifRegionsLen - preATGBP; ++iter )
		{
			if( isSignalInPos( motifRegions[i].data(), signal, iter, mutBp ) )
			{	
				++count;
				char *seq1 = &motifRegions[i][iter - preSigBP];
				int k = 0;
				for(; k < WM.getColum(); ++k )
				++WM( seq1[k]-'0', k );
				++distri[iter];
			}
		}
	}
	

	WM.toBeAveraged().toBeLoged();
	int iter = 0;
	for(; iter < distri.size(); ++iter )
		distri[iter] = ( distri[iter] / count );
	return std::make_pair< Ve_D, M_D >( distri, WM );
}

Ve_Location MED_start::test()
{
	for( Ve_Loca::iterator iter = rvsLocations.begin(); iter != rvsLocations.end(); ++iter )
	{
		const char* seq = iter->isPositive ? positiveSeq.data() : negtiveSeq.data();
		x2LongestORF( seq, iter->location );
		iter->ORFLength = iter->location.second - iter->location.first;
		iter->ATGIndexInORF = 0;
		reviseTIS( seq, *iter, XTGIndex );
	}
	return rvsLocations;
}

void  MED_start::prepareForTesting()
{
	Ve_Location rvsLocations = this->rvsLocations;
	Ve_D XTGIndexTmp = XTGIndex;
	for( Ve_Loca::iterator iter = rvsLocations.begin(); iter != rvsLocations.end(); ++iter )
	{	
		const char* seq = iter->isPositive ? positiveSeq.data() : negtiveSeq.data();
		reviseTIS( seq, *iter, XTGIndexTmp );
		if( iter->ORFLength >= learnORFLenLmt )
		{
			if( seq[iter->location.first] == '0' )
				++XTGIndex[iter->ATGIndexInORF];
			if( iter->ATGIndexInORF <= 0 )
			{
				switch( seq[iter->location.first] )
				{
				case '0' : XTG[0] += 1; break;
				case '1' : XTG[1] += 0.8; break;
				case '2' : XTG[2] += 0.8; break;
				case '3' : XTG[3] += 0.8; 
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
	XTGIndexTmp1 = XTGIndex;
	i =XTGIndexTmp.size();
	double all = 0;
	for( i = 0; i < 4; ++i )
		all += XTG[i];

	XTGTmp = XTG;
	for( i = 0; i < 4; ++i )
	{
		XTGTmp[i] /= all;
		ATGWM( i, preATGBP ) = log( XTG[i] / all );
		XTG[i] = 0;
	}
}


void MED_start::formATGWM()
{
	ATGWM = M_D( 4, ATGWM.getColum(), 0.0001 );
	int i = 0;
	for(; i < rvsLocations.size(); ++i )
	{
		const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() : negtiveSeq.data() ;
		if( rvsLocations[i].ORFLength > learnORFLenLmt 
			&& rvsLocations[i].location.first -preATGBP >= 0 )
		{
			int iter = 0;
			for(; iter < ATGWM.getColum(); ++iter )
				++ATGWM( seq[rvsLocations[i].location.first + iter-preATGBP] - '0', iter );
		}
	}
	ATGWM.toBeAveraged().toBeLoged();
}

M_D MED_start::getATGWM()
{
	M_D ATGWM = M_D( 4, 100, 0.0001 );
	int i = 0;
	for(; i < rvsLocations.size(); ++i )
	{
		const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() : negtiveSeq.data() ;
		if( rvsLocations[i].ORFLength > learnORFLenLmt 
			&& rvsLocations[i].location.first -50 >= 0 )
		{
			int iter = 0;
			for(; iter < ATGWM.getColum(); ++iter )
				++ATGWM( seq[rvsLocations[i].location.first + iter-60] - '0', iter );
		}
	}
	return ATGWM.toBeAveraged();
}

void MED_start::findHitMotifs( int mutNum )
{
	ifstream in( (m_motifsPath+"motifs.txt").data() );
	assert( in.good() );
	std::map<double,std::string,std::greater<double> > sig;
	while( !in.eof() )
	{
		Str signal;
		double times;
		in>>signal>>times;
		if( !signal.empty() )
			sig[times] = signal;//see MED_start.err
	};
	in.close();
	system( Str("rm \"" + m_motifsPath+"motifs.txt\"" ).data() );
	std::map<double,std::string,std::greater<double> >::iterator iter = sig.begin();
	std::map<double,std::pair<Ve_D,Str>,std::greater<double> > kernel;
	int i = 0;
	for(; i < 5 && iter != sig.end(); ++iter, ++i )
	{
		Ve_D tmp = signalDistri( iter->second.data(), mutNum );
		kernel[GetDerivate( tmp )] = std::make_pair<Ve_D,Str>( tmp, iter->second );
		int j = 1;

		candidateMotifs.push_back( iter->second );
	}

	std::map<double,std::pair<Ve_D,Str>,std::greater<double> >::iterator kernelIter = kernel.begin();
	int j = 0;
	for(; j < 3; ++j ) 
	{
		if( kernelIter == kernel.end() )
			break;
		std::map<double,std::pair<Ve_D,Str>,std::greater<double> >::iterator tmp = kernelIter;
		hitMotifs.push_back( tmp->second.second );
		spacerLenDistriTmp.push_back( tmp->second.first );
		if( fabs( tmp->first - (++kernelIter)->first ) > 0.006 )
			break;
	}
}

Ve_D MED_start::signalDistri(  Str signal , int mutNum )
{
	int disCount = 0;
	Ve_D distri( motifRegionsLen - sigLen + 1, 0 );
 	int i = 0;
 	for(; i < motifRegions.size(); ++i )
	{
		int iter = 0;
		for(; iter + sigLen <= motifRegionsLen; ++iter )
		{
			if( isSignalInPos( motifRegions[i].data(), signal.data(), iter, mutNum ) )
			{	
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

int MED_start::geneStatus( Str& seq )
{
	int len = seq.size();
	if( len%3 != 0 )
	{
		cout<<"LENNOT3TIMES ERR "<<endl;
		return LENNOT3TIMES;
	}
	int i = 0;
	for(; i < len - 2; i += 3 )
	{
		int subStr = ((seq[i] - 48)<<4) + ((seq[i+1] - 48)<<2) + seq[i + 2] - 48;
		if( i == 0 && subStr != 14 && subStr != 46 && subStr != 62 ) 
		{
			cout<<"TISERR at "<<i<<"\t, as "<<seq[i]<<seq[i+1]<<seq[i+2]<<endl;
			return TISERROR;
		}
		if( i == len - 3 && subStr != 48 && subStr != 56 && subStr != 50 )
		{
			cout<<"STPERR at "<<i<<"\t, as "<<seq[i]<<seq[i+1]<<seq[i+2]<<endl;
			return STPERROR;
		}
		if( i != len - 3 && (subStr == 48 || subStr == 56 || subStr == 50 ) )
		{
			cout<<"STPINFRAME at "<<i<<"\t, as "<<seq[i]<<seq[i+1]<<seq[i+2]<<endl;
			return STPINFRAME;
		}
	}
	return ISAGENE;
}

int MED_start::theIthATGInORF( Str& seq, int ATGPos, int ORFBeg, int ORFEnd )
{
	int index = 0;
	int i = ORFBeg;
	for(; i < ORFEnd - 1; i += 3 )
	{
		int subStr = ((seq[i] - 48)<<4) + ((seq[i+1] - 48)<<2) + seq[i + 2] - 48;
		if( subStr == 14 || subStr == 46 || subStr == 62 )
		{
			++index;
			if( i == ATGPos )
				return index;
		}
	}
	return -1;
}


int MED_start::containedATGNum( const char* seq, int beg, int end )
{
	int index = 0;
	int i = beg;
	for(; i < end - 1; i += 3 )
	{
		int subStr = ((seq[i] - 48)<<4) + ((seq[i+1] - 48)<<2) + seq[i + 2] - 48;
		if( subStr == 14 || subStr == 46 || subStr == 62 )
			++index;
	}
	return index;
}

double MED_start::reviseTIS( const char* seq, Location_T& iter, Ve_D& XTGIndexTmp )
{
	int STP = iter.location.second;
	int ATGPosition = iter.location.first;
	int count = iter.ATGIndexInORF;
	double score = -100000;
	int predictATGPos = ATGPosition;
	do
	{
		ATGPosition = getNextPhaseTIS( seq, ATGPosition );
		if( ATGPosition  == -1 || ATGPosition >=  STP 
			|| count >= XTGIndex.size() || STP - ATGPosition <= 87 )
			break;
		
		double tmp = -100000;
		int pos = ATGPosition - motifRegionsLen;
		for(; pos <= ATGPosition - sigLen; ++pos )
		{
			int i = 0;
			for(; i < spacerLenDistri.size(); ++i )
			{
				Pa_D_D rbsScore = getRBSScore( seq, pos, ATGPosition, spacerLenDistri[i], SigWM[i] );
				if( rbsScore.first + rbsScore.second > tmp )
					tmp = rbsScore.first + rbsScore.second;
			}
		}
		int i = 0;
		for(; i < ATGWM.getColum(); ++i )
		{
			if( ATGPosition + i - preATGBP > 0 )
				tmp += ATGWM( seq[ ATGPosition + i - preATGBP] - '0', i );
		}
		tmp += log(XTGIndexTmp[count]);
		if( tmp > score )
		{
			score = tmp;
			predictATGPos = ATGPosition;
			iter.ATGIndexInORF = count;
		}
			
		ATGPosition += 3;
		++count;
	}while( true );
	
	iter.location.first = predictATGPos;
	iter.RBSScore = score;
	iter.ORFLength = STP - predictATGPos;
	return iter.RBSScore;
}

double MED_start::reviseTIS( const char* seq, Location_T& iter )
{
	Ve_D XTGIndexTmp = XTGIndex;
	int STP = iter.location.second;
	int ATGPosition = iter.location.first;
	int count = iter.ATGIndexInORF;
	double score = -100000;
	int predictATGPos = ATGPosition;
	do
	{
		ATGPosition = getNextPhaseTIS( seq, ATGPosition );
		if( ATGPosition  == -1 || ATGPosition >=  STP 
			|| count >= XTGIndex.size() || STP - ATGPosition <= 87 )
			break;

		double tmp = -100000;
		int sigPos = -1;
		int pos = ATGPosition - motifRegionsLen;
		for(; pos <= ATGPosition - sigLen; ++pos )
		{
			int i = 0;
			for(; i < spacerLenDistri.size(); ++i )
			{
				Pa_D_D rbsScore = getRBSScore( seq, pos, ATGPosition, spacerLenDistri[i], SigWM[i] );
				if( rbsScore.first + rbsScore.second > tmp ){
					tmp = rbsScore.first + rbsScore.second;
					iter.hitMotif = hitMotifs[i];
					sigPos = pos;
				}
			}
		}
		iter.hitMotif = SequenceTransform_T::digital2CharSeq( iter.hitMotif );
		Str tmpMotif(sigLen+6,'5');
		int i = 0 ;
		for(; i < tmpMotif.size(); ++i )
			if( sigPos+i-3 < seqLen )
				tmpMotif[i]=( seq[sigPos+i-3] );
			else
				break;
		iter.motif = SequenceTransform_T::digital2CharSeq( tmpMotif );
		double ATGScore = 0;
		for( i = 0; i < ATGWM.getColum(); ++i )
		{
			if( ATGPosition + i - preATGBP > 0 )
				ATGScore += ATGWM( seq[ ATGPosition + i - preATGBP] - '0', i );
		}

		ATGScore += log(XTGIndexTmp[count]);
		tmp += ATGScore;
		if( tmp > score )
		{
			score = tmp;
			predictATGPos = ATGPosition;
			iter.ATGIndexInORF = count;
			iter.ATGScore = ATGScore;
			iter.motifPosition = sigPos - ATGPosition;
		}
			
		ATGPosition += 3;
		++count;
	}while( true );
	iter.location.first = predictATGPos;
	iter.RBSScore = score;
	iter.ORFLength = STP - predictATGPos;
	iter.ATGNumCounted = containedATGNum( seq, iter.location.first, iter.location.second );

	return iter.RBSScore;
}
 Str parafile;
Ve_Location MED_start::doIteration()
{
	kerneWordsDistriAndWM(mutBp);
	formATGWM();
	prepareForTesting();
	return test();
}


bool MED_start::hasSignal( const char* seq, const char* signal, int motifRegionsLen, int mutBp )
{
	int iter = 0;
	for(; iter + sigLen <= motifRegionsLen - posSigBP; ++iter )
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

void MED_start::setHitmotifs(Ve_Str hitMotifs1)
{
	this->hitMotifs = hitMotifs1;
}

Ve_Str MED_start::getHitMotifs()
{
	return hitMotifs;
}

void MED_start::setCandidateMotifs(Ve_Str candidateMotifs)
{
	this->candidateMotifs = candidateMotifs;
}

Ve_Str MED_start::getCandidateMotifs()
{
	return candidateMotifs;
}

void MED_start::findHitMotif1()
{
	findMotifs(0);
	findHitMotifs(mutBp);
}

void MED_start::setParameters( CParametersDlg* m_parameters )
{
	learnORFLenLmt = m_parameters->m_ORFLen;
	motifRegionsLen = m_parameters->m_sigRegion;
	sigLen = m_parameters->m_sigLen;
	mutBp = m_parameters->m_sigMutNum;
	preSigBP = m_parameters->m_sigPreNum;
	posSigBP = m_parameters->m_sigPosNum;
	preATGBP =  m_parameters->m_atgPreNum;
	posATGBP = m_parameters->m_atgPosNum;
}

void MED_start::resetLearnedRst()
{
	ATGWM = M_D( 4, preATGBP + posATGBP + 3,0.0001 );
	XTG = Ve_D( 4, 0.01 );
	XTGIndex = Ve_D( 6, 1 );
	motifRegions.erase( motifRegions.begin(), motifRegions.end() );
	rvsLocations.erase( rvsLocations.begin(), rvsLocations.end() );
	hitMotifs.erase( hitMotifs.begin(), hitMotifs.end() );
	candidateMotifs.erase( candidateMotifs.begin(), candidateMotifs.end() );
	SigWM.erase( SigWM.begin(), SigWM.end() );
	spacerLenDistri.erase( spacerLenDistri.begin(), spacerLenDistri.end() );
	candidateSpacerLenDistri.erase( candidateSpacerLenDistri.begin(),
		candidateSpacerLenDistri.end() );
}

void MED_start::resetLearnedRst(bool)
{
	ATGWM = M_D( 4, preATGBP + posATGBP + 3,0.0001 );
	XTG = Ve_D( 4, 0.01 );
	XTGIndex = Ve_D( 6, 1 );
	motifRegions.erase( motifRegions.begin(), motifRegions.end() );
	hitMotifs.erase( hitMotifs.begin(), hitMotifs.end() );
	candidateMotifs.erase( candidateMotifs.begin(), candidateMotifs.end() );
	SigWM.erase( SigWM.begin(), SigWM.end() );
	spacerLenDistri.erase( spacerLenDistri.begin(), spacerLenDistri.end() );
	candidateSpacerLenDistri.erase( candidateSpacerLenDistri.begin(),
		candidateSpacerLenDistri.end() );
}

void MED_start::setMotifRegions()
{
	int i = 0;
	for(; i < rvsLocations.size(); ++i )
	{
		if( rvsLocations[i].ORFLength >= learnORFLenLmt )
		{
			Location_T tmp = rvsLocations[i];
			Str& seq = rvsLocations[i].isPositive ? positiveSeq : negtiveSeq;
			if( rvsLocations[i].location.first - motifRegionsLen >= 0 )
			{	
				Location_T tmp = rvsLocations[i];
				Str tmpSeq = seq.substr( tmp.location.first, tmp.location.second - tmp.location.first + 3 );
				int status;
				if( ( status = geneStatus( tmpSeq ) ) != ISAGENE )
				{
				}

				Str temp = seq.substr( rvsLocations[i].location.first - motifRegionsLen,
					motifRegionsLen );
				if( temp.size() == motifRegionsLen )
					motifRegions.push_back( temp );
			}
		}
	}
}

Ve_Location MED_start::getRvsLocations()
{
	int i = 0;
	for(; i < rvsLocations.size(); ++i ){
		rvsLocations[i].fileForm2Location( seqLen );
	}
	int iterNum = 0;
	for( double iterationSTPTag = 1000000; ; )
	{
		cout<<"do iteration "<<iterNum+1<<endl;
		setMotifRegions();
		findHitMotif1();
		rvsLocations = doIteration();

		double percent = fabs( iterationSTPTag /  getSTPTag() );
		if( ( percent > 0.99 && percent <1.01) || ++iterNum == 20 )
			break;
		else
			iterationSTPTag = getSTPTag();
		//prepair for next iteration.
		resetLearnedRst(true);
	}
	//for tmp
	unsigned long int sigLen = 6;
	Ve_I initSignal( sigLen, 0 );
	initSignal[sigLen-1] = -1;	
	int upBoundary = pow( 4, sigLen );
	ofstream out( "G:\\arhu\\tuition\\c++\\logo\\Release\\logoWM1.txt" );

	int v = 0;
	for(; v < upBoundary; ++v )
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
		Str tmp( sigLen, '0' );
		int k = 0;
		for(; k < tmp.size(); ++k )
		{
			switch( initSignal[k] )
			{
			case 1 :
				tmp[k] = '1';
				break;
			case 2:
				tmp[k] = '2';
				break;
			case 3:
				tmp[k] = '3';
			}
		}
		M_D tmpWM( 4, tmp.size() );
		for( i = 0; i < rvsLocations.size(); ++i ){
			const char* seq = rvsLocations[i].isPositive ? positiveSeq.data() 
				: negtiveSeq.data();
			reviseTIS( seq, rvsLocations[i] );
			if( rvsLocations[i].location.first - 50 >= 0 ){
				rvsLocations[i].hosters = findHostStr( seq, rvsLocations[i].location.first - 50,
					rvsLocations[i].location.first, tmp );
				int k = 0;
				for(; k < rvsLocations[i].hosters.first.size(); ++k ){
					int m = 0;
					for(; m < tmp.size(); ++m )
						++tmpWM( rvsLocations[i].hosters.first[k].first[m] - '0', m );
				}
			}
			rvsLocations[i].location2FileForm( this->seqLen );
		}
		return rvsLocations;
		int times = 10;
		int p = 0;
		for(; p < times; ++p ){
			M_D tmpWMR( 4, tmp.size() );
			for( i = 0; i < rvsLocations.size(); ++i ){
				if( rvsLocations[i].hosters.first.size() > 1 ){
					int max = 0, index = -1;
					int k = 0;
					for(; k < rvsLocations[i].hosters.first.size(); ++k ){
						double score = 0;
						int m = 0;
						for(; m < tmp.size(); ++m ){
							score += tmpWM( rvsLocations[i].hosters.first[k].first[m] - '0', m );
						}
						if( score > max ){
							max = score ;
							index = k;
						}
					}
					if( p == times - 1 )
						rvsLocations[i].hosters.first = Ve_Pa_Str_I(1,rvsLocations[i].hosters.first[index]);
					else
					{
						int m = 0;
						for(; m < tmp.size(); ++m )
							++tmpWMR( rvsLocations[i].hosters.first[index].first[m] - '0', m );
					}

				}
				if( p != times - 1 && rvsLocations[i].hosters.first.size() == 1 )
				{
					int m = 0;
					for(; m < tmp.size(); ++m ){
						++tmpWMR( rvsLocations[i].hosters.first[0].first[m] - '0', m );
				}
				}
			}
			if( p != times - 1 )
				tmpWM = tmpWMR.toBeAveraged();
		}
		{
			std::map<int,int> pos_count;
			double sum = rvsLocations.size();
			for( i = 0; i < rvsLocations.size(); ++i )
				++pos_count[-rvsLocations[i].hosters.first[0].second -sigLen];
			out<<SequenceTransform_T::digital2CharSeq(tmp)<<endl;
			(tmpWM.toBeLogo()).matrixOut(out);
			out<<pos_count.size()<<endl;
			for( std::map<int,int>::iterator iter = pos_count.begin(); iter != pos_count.end(); ++iter )
				out<<iter->second / sum<<endl;
			out<<endl;
		}
	}
	return rvsLocations;
}


void MED_start::saveParameter(std::ofstream &out)
{
	out<<"< "<<genomeID<<endl;
	out<<"XTG_Index "<<XTGIndex.size()<<endl;
	int i = 0;
	for(; i <XTGIndex.size(); ++i )
		out<<setw(10)<<XTGIndex[i];
	out<<endl;
	out<<"XTG_WM"<<endl;
	M_D ATGWMTmp = ATGWM;
	(ATGWMTmp.toExp()).matrixOut( out );
	out<<endl;
	out<<"motif_spacer(position of first nucl of signal vs posibility)_WM  "<<spacerLenDistri.size()<<endl;
	for( i = 0; i < hitMotifs.size(); ++i ){
		out<<SequenceTransform_T::digital2CharSeq(hitMotifs[i])
			<<" "<<spacerLenDistri[i].size()<<endl;
		int j = 0;
		for(; j < spacerLenDistri[i].size(); ++j )
			out<<setw(10)<<j-motifRegionsLen <<setw(10)<<spacerLenDistriTmp[i][j]<<endl;
		M_D WMTmp = SigWM[i];
		(WMTmp.toExp()).matrixOut(out);
	}
}
