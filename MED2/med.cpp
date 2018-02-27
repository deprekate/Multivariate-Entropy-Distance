#include "MED.h"
#include "OftenUsedOperatLib.h"
#include "time.h"
#include "SequenceTransform.h"

double ratio5, ratio15;
const int preATGBP= 50, posATGBP = 15;
int TTmpIter =0, FTmpIter = 0;
int centerORFLength = 0;
double TN = 0, FN =0;
int indexNUm = 10; //Only consider the first 10 candiate TISs
int initFIndex = 0;
int initTIndex = 0;

MED::MED( Str& seq )
:trueORFCenter( 10000, Ve_D(20) ), falseORFCenter( 100000, Ve_D(20) )
,geneInfo( seq ), seq( NULL ), seqLen( seq.size() ), med_start(NULL)
{
	if( !isGCRich )
		ratio5 = 0.5, ratio15 = 1.5;
	else
		ratio5 = 0.7, ratio15 = 1.6;
}

std::vector<ORF_T> MED::getGeneLocation()
{
	assert( seq == NULL && "No update sequence, please rebuild MED object" );
	if( isGCRich )
		readCenters( "EDPCentersGC.txt" );
	else{
		readCenters( "EDPCenters.txt" );
	}
	std::list<ORF_T> ORFSet;
	getORFSet( ORFSet );
	cout<<"Number of exactcted ORFs: "
		<<ORFSet.size()<<endl;
	decorateInitPoint( ORFSet );
	return identifyGeneLocation( ORFSet );
} 


void MED::getORFSet( std::list<ORF_T>& ORFSet )
{
	const Ve_Location& location = geneInfo.getORFLocation();
	for( Ve_Location::const_iterator iter = location.begin(); iter != location.end(); ++iter )
	{
		ORF_T ORF;
		ORF.location = iter -> location;
		ORF.isPositive = iter -> isPositive;
		ORF.motifWM = iter->motifWM;
		ORF.motifSpacer = iter -> motifSpacer;
		ORF.XTGWM = iter -> XTGWM;
		ORF.XTGIndex = iter -> XTGIndex;
		if( ORF.isPositive )
			seq = geneInfo.getPositiveSeq().data();
		else
			seq = geneInfo.getNegtiveSeq().data();
		ORF.EDP = aminoEDPStrengthen
			( ORF.location.second, ORF.location.first, ORF.location.second );
		ORF.ORFLength = fabs( ORF.location.second - ORF.location.first );
		ORFSet.push_back( ORF );
	}
}

void MED::readCenters(Str fileName ){
	ifstream in( fileName.data() );
	if( !in.good() ){
		cout<<fileName<<" not found"<<endl;
		exit(1);
	}
	int tNum, fNum;
	in>>tNum>>fNum;
	int i = 0;
	for(; i < tNum; ++i ){
		Str line;
		in>>line;//skipp over tag
		int j = 0;
		for(; j < 20; ++j ){
			in>>trueORFCenter[i][j];
		}
	}
	cout<<endl;
	for( i = 0; i < fNum; ++i ){
		Str line;
		in>>line;//skipp over tag
		int j = 0;
		for(; j < 20; ++j ){
			in>>falseORFCenter[i][j];
		}
	}
	TBeg = FBeg = 0;
	TEnd = tNum, FEnd = fNum;
}

void MED::decorateInitPoint( std::list<ORF_T>& ORFSet )
{
	Ve_Ve_D trueInit, falseInit;
	Ve_D tCount( TEnd,0), fCount(FEnd,0);
	int j = 0;
	for(; j < TEnd; ++j ){
		trueInit.push_back(trueORFCenter[j]);
	}
	for( j = 0; j < FEnd; ++j ){
		falseInit.push_back(falseORFCenter[j]);
	}
	std::list<ORF_T>::iterator iter;
	for( iter = ORFSet.begin(); iter != ORFSet.end(); ++iter)
	{
		if( fabs( iter->location.second - iter->location.first ) < centerORFLength )
			continue;
		double ratio = ORFDistanceRatio( iter ); 
		iter->Dis2TCenters = 100000, iter->Dis2FCenters = 100000;

		if( ratio < ratio5 )
		{
			++tCount[initTIndex];
			int i = 0;
			for(; i < 20; ++i )
				trueInit[initTIndex][i] += iter -> EDP[i];
		}
		else if( ratio > ratio15 )
		{
			++fCount[initFIndex];
			int i = 0;
			for(; i < 20; ++i )
				falseInit[initFIndex][i] += iter -> EDP[i];
		}
	}
	for( j = 0; j < trueInit.size(); ++j ){
		int i = 0;
		for(; i < 20; ++i){
			trueInit[j][i] /= tCount[j] + 1;
		}
		trueORFCenter[j] = trueInit[j];
	}

	for( j = 0; j < falseInit.size(); ++j ){
		int i = 0;
		for(; i < 20; ++i){
			falseInit[j][i] /= fCount[j] + 1;
		}
		falseORFCenter[j] = falseInit[j];
	}
}

void MED::updateATGWMAndIndex(const std::list<ORF_T>& trueORF){
	ATGWM = M_D( 4, preATGBP+3+ posATGBP, 0.0001 );
	M_D ATGWMBK( 4, ATGWM.getColum(), 0.0001 );
	std::list<ORF_T>::const_iterator iter = trueORF.begin();
	Ve_D indexTmp( indexNUm, 1 );
	for( ; iter!=trueORF.end(); ++iter){
		const ORF_T& rvsLocation = *iter;
		if( fabs( rvsLocation.location.first-rvsLocation.location.second ) 
			< centerORFLength )
			continue;
		++indexTmp[rvsLocation.ATGIndexInORF];
		const char* seq = rvsLocation.isPositive ? geneInfo.positiveSeq.data() 
			: geneInfo.negtiveSeq.data() ;
		if(	rvsLocation.location.first -preATGBP >= 0 ){
			int iter = 0;
			for(; iter < ATGWM.getColum(); ++iter )
				++ATGWM( seq[rvsLocation.location.first + iter-preATGBP] - '0', iter );
		}
		int ATGPos = geneInfo.getNextPhaseTIS(seq, rvsLocation.location.first+3);
		if( ATGPos > rvsLocation.location.first && ATGPos-preATGBP >= 0){
			int iter = 0;
			for(; iter < ATGWM.getColum(); ++iter )
				++ATGWMBK( seq[ATGPos + iter-preATGBP] - '0', iter );
		}
	}
	ATGWM.toBeAveraged();
	ATGWMBK.toBeAveraged();
	int i = 0;
	for(; i < ATGWM.getLine(); ++i ){
		int j = 0;
		for( ; j < preATGBP; ++j ){
			ATGWM( i, j ) /= geneInfo.background[i];
		}
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
	ATGWM.toBeLoged();//*/
	double sum = 0;
	for( i = 0; i < indexTmp.size(); ++i )
		sum += indexTmp[i];
	for( i = 0; i < indexTmp.size(); ++i )
		indexTmp[i]/=sum;
	if( indexTmp[0] < 0.98 )
		index = indexTmp;
}

double MED::relocateTISByATGWManIndex( const char* seq, Location_T& iter ){
	int STP = iter.location.second;
	int ATGPosition = iter.location.first;
	int count = 0;
	double score = -100000;
	int predictATGPos = ATGPosition;
	int maxATGCount = index.empty() ? indexNUm : index.size();
	do{
		ATGPosition = geneInfo.getNextPhaseTIS( seq, ATGPosition );
		if( ATGPosition  == -1 || ATGPosition >=  STP 
			|| count >=  maxATGCount || STP - ATGPosition < 60 ){
			break;
		}
		double ATGScore = 0;
		int i = 0;
		for(; i < ATGWM.getColum(); ++i ){
			if( ATGPosition + i - preATGBP > 0 )
				ATGScore += ATGWM( seq[ ATGPosition + i - preATGBP] - '0', i );
		}
		if( !index.empty() )
			ATGScore += log(index[count]);
		if(  ATGPosition >  90){
			double noncodingP = 1-getCodProd( int( getGB(seq, ATGPosition - 90, ATGPosition) ) );
			double codingP = getCodProd( int( getGB(seq, ATGPosition, ATGPosition+90) ) );
			ATGScore += log(noncodingP)+log(codingP);
		}//*/
		if( ATGScore > score ){
			score = ATGScore;
			predictATGPos = ATGPosition;
			iter.ATGIndexInORF = count;
		}
		ATGPosition += 3;
		++count;
	}while( true );
	iter.location.first = predictATGPos;
	iter.ATGScore = score;
	iter.ORFLength = STP - predictATGPos;
	return score;
}

void MED::relocateTISsByATGWManIndex( std::list<ORF_T>& ORFs ){
	for( std::list<ORF_T>::iterator it = ORFs.begin(); it != ORFs.end(); ++it ){
		const char* seq = it->isPositive 
			? geneInfo.positiveSeq.data() : geneInfo.negtiveSeq.data();
		relocateTISByATGWManIndex( seq, *it );
	}
}

void MED::refineATGWMAndIndex(std::list<ORF_T>& trueORF){
	double totalNumber = trueORF.size();
	while( true ){
		double numberOfChanged = 0;
		for( std::list<ORF_T>::iterator it = trueORF.begin(); it != trueORF.end(); ++it ){
			Location_T& iter = *it;
			const char* seq = iter.isPositive ? 
				geneInfo.positiveSeq.data():geneInfo.negtiveSeq.data();
			int origATG = iter.location.first;
			x2LongestORF( seq, iter.location );
			relocateTISByATGWManIndex( seq, iter );
			if( origATG != iter.location.first ){
				++numberOfChanged;
			}
		}
		if( numberOfChanged/totalNumber < 0.15 )
			break;
		updateATGWMAndIndex(trueORF);	
		updateGBDistri( trueORF );
	}
}

void MED::reEvaluateSeeds(std::list<ORF_T>& trueORF
						  , std::list<ORF_T>& falseORF, std::list<ORF_T>& remain 
						  , double threshold ){
	if( isGCRich ){
		int addFalseByATG = 0, addFalseByGB = 0, addFalseByBoth = 0;

		updateATGWMAndIndex( trueORF );
		updateGBDistri(trueORF);
		refineATGWMAndIndex( trueORF );

		std::vector<std::pair<double, double> > scores;
		std::list<ORF_T>::iterator it;
		for(it = trueORF.begin(); it != trueORF.end(); ++it ){
			const char* seq = it->isPositive 
				? geneInfo.positiveSeq.data() : geneInfo.negtiveSeq.data();
			relocateTISByATGWManIndex( seq, *it );
			double gb = getGB( seq, it->location.first, it->location.second );
			it->GBScore = gb;
			if( fabs(it->ATGScore)<1000 && fabs(it->GBScore)<1000 )
				scores.push_back(std::make_pair( it->ATGScore, it->GBScore ) );
		}
		std::pair<double,double> aver = std::make_pair(0,0);
		M_D S(2,2, 0), Inv_S(2,2,0);
		int i = 0;
		for(; i < scores.size(); ++i ){
			aver.first += scores[i].first;
			aver.second += scores[i].second;
		}
		aver.first /= scores.size();
		aver.second /= scores.size();
		//calculate derivation
		for( i = 0; i < scores.size(); ++i ){
			std::pair<double,double> tmp;
			tmp.first = scores[i].first - aver.first;
			tmp.second = scores[i].second - aver.second;
			S(0,0) += tmp.first*tmp.first;
			S(0,1) += tmp.first*tmp.second;
			S(1,0) += tmp.second*tmp.first;
			S(1,1) += tmp.second*tmp.second;
		}
		for( i = 0; i < S.getLine(); ++i )
		{
			int j = 0;
			for(; j < S.getColum(); ++j )
				S(i,j) /= (scores.size() - 1 );
		}
		double b2ac = S(0,1)*S(0,1) - S(0,0)*S(1,1);
		Inv_S(0,0) = -S(1,1)/b2ac;
		Inv_S(0,1) = S(0,1)/b2ac;
		Inv_S(1,0) = S(1,0)/b2ac;
		Inv_S(1,1) = -S(0,0)/b2ac;

		double x = 0;
		for(  it = remain.begin(); it != remain.end(); ){
			const char* seq = it->isPositive 
				? geneInfo.positiveSeq.data() : geneInfo.negtiveSeq.data();
			relocateTISByATGWManIndex( seq, *it );
			it->GBScore = getGB( seq, it->location.first, it->location.second );
			double x1 = it->ATGScore-aver.first, x2 = it->GBScore-aver.second;
			double dis = Inv_S(0,0)*x1*x1+2*Inv_S(0,1)*x1*x2+Inv_S(1,1)*x2*x2;
			if( x2 > 0  ){
				++it;
			}
			else if(dis>threshold && !it->certainCDS){
				if( it->EDPScore <( (1-x)*ratio5+x*ratio15) )//{ 
					it->canBeCenter = false;
				falseORF.push_back( *it );
				it = remain.erase( it );
				++addFalseByBoth;//}
			}//*/
			else{
				++it;
			}
		}//*/
	}
	//update center
	TTmpIter =0, FTmpIter = 0;
	std::list<ORF_T>::iterator it;
	for( it = trueORF.begin(); it != trueORF.end(); ++it ){
		if( it->ORFLength > centerORFLength )
			std::copy( it -> EDP.begin(), it -> EDP.end()
			, trueORFCenter[++TTmpIter].begin() );
	}
	for( it = falseORF.begin(); it != falseORF.end(); ++it ){
		if( it->ORFLength > centerORFLength 
			&& it->canBeCenter )
			std::copy( it -> EDP.begin(), it -> EDP.end()
			, falseORFCenter[++FTmpIter].begin() );
	}
}

std::vector<ORF_T> MED::identifyGeneLocation						
( std::list<ORF_T>& ORFSet )
{
	time_t start,finish;
	time(&start);
	cout<<"Select root ORFs and first running of TIS model..."<<endl;
	cout<<std::setw(10)<<"iteration"<<std::setw(25)<<"Selected as coding "<<std::setw(25)
		<<"rejected as non-coding"<<std::endl;

	Ve_Location rvsLocations;
	std::list<ORF_T> trueORF, falseORF;
	int iterNum = 0;
	do
	{
		std::list<ORF_T>::iterator iter = ORFSet.begin();
		TTmpIter = TEnd - 1, FTmpIter = FEnd - 1;
		TN = trueORF.size(), FN = falseORF.size();
		for( ; iter != ORFSet.end(); )
		{
			double ratio = ORFDistanceRatio( iter );
			ORF_T tmp = *iter;
			iter->EDPScore = ratio;
			if( ratio > ratio5 && ratio < ratio15 )
				++iter;
			else 
			{
				if( ratio <= ratio5 )
				{	
					trueORF.push_back( *iter );
					rvsLocations.push_back( *iter );//up_casting
				}
				else
				{
					falseORF.push_back( *iter );					
				}
				iter = ORFSet.erase( iter );
			}
		}
		reEvaluateSeeds( trueORF, falseORF, ORFSet, 20000);
		TBeg = TEnd, FBeg = FEnd;
		TEnd = TTmpIter + 1, FEnd = FTmpIter + 1;
		cout<<std::setw(8)<<(iterNum ++)<<std::setw(25)<<(trueORF.size()-TN)<<std::setw(20)
			<<std::setw(25)<<(falseORF.size()-FN)<<std::endl;
	}while( TEnd != TBeg || FEnd != FBeg );
	reEvaluateSeeds( trueORF, falseORF, ORFSet, 20);

	cout<<"Selected root coding ORFs number: "<<trueORF.size()<<endl;
	cout<<"Rejected non-coding ORFs number: "<<falseORF.size()<<endl;
	cout<<"Remaining ORFs number: "<<ORFSet.size()<<endl;
	time(&finish);
	if( isGCRich ){
		std::list<ORF_T> medResult;
		//prepare normalized scores for fisher discriminant
		Ve_D edpS, atgS, gbS;	
		std::list<ORF_T>::iterator iter;
		for( iter = trueORF.begin()
			; iter != trueORF.end(); ++iter ){
			if( fabs(iter->EDPScore) < 1000 )
				edpS.push_back( iter->EDPScore );
			if( fabs(iter->ATGScore) < 1000 )
				atgS.push_back( iter->ATGScore );
			if( fabs(iter->GBScore < 1000 ) )
				gbS.push_back( iter->GBScore );
			iter->certainCDS = true;
			medResult.push_back( *iter );
		}
		double edpA = GetAver(edpS), edpD = GetDerivate(edpS)
			, atgA = GetAver(atgS), atgD = GetDerivate(atgS)
			, gbA = GetAver(gbS), gbD = GetDerivate(gbS);

		//fisher descriminant! between edp, atg and GB
		cout<<"Select coding ORFs from remianing ORFs" 
			" using Fisher discriminant functions..."<<endl;//
		double atg_x1 =0.34, edp_y1 = 9.47, atg_x2 = -3.6, edp_y2 = 2.75;
		double k1 = (edp_y1-edp_y2)/(atg_x1-atg_x2);
		double b1 =  edp_y1-k1*atg_x1;
		int addFalseByFisher = 0;
		
		double edp_x1 = 9.1, gb_y1 = -0.5, edp_x2 = 0.76, gb_y2 = -3.75;
		double k2 = (gb_y2-gb_y1)/(edp_x2-edp_x1);
		double b2 = gb_y1-k2*edp_x1;

		for(  iter = ORFSet.begin()
			; iter != ORFSet.end(); ++iter ){
			double atg_x = (iter->ATGScore-atgA)/atgD
				, edp = (iter->EDPScore-edpA)/edpD
				, gb_y = (iter->GBScore-gbA)/gbD;
			if( edp - k1*atg_x - b1 < 0 
				&& gb_y -k2*edp - b2 > 0 
				)
			{
				medResult.push_back( *iter );
			}
			else{
				++addFalseByFisher;
			}
		}//*/
		//normalize ATG, EDP and GB Score;
		std::vector<ORF_T> rst;
		for( std::list<ORF_T>::iterator it = medResult.begin()
			; it != medResult.end(); ++it ){
			it->EDPScore = (it->EDPScore-edpA)/edpD;
			it->GBScore = (it->GBScore-gbA)/gbD;
			it->ATGScore = (it->ATGScore-atgA)/atgD;
			const char* seq = it->isPositive 
				? geneInfo.positiveSeq.data():geneInfo.negtiveSeq.data();
			it->location2FileForm(seqLen );
			rst.push_back( *it );
		}
		cout<<"Remaining: "<<rst.size()<<endl;
		return rst;
	}else{
		std::list<ORF_T>::iterator iter = ORFSet.begin();
		for( ; iter != ORFSet.end(); ++iter )
		{
			iter->EDPScore = iter->Dis2TCenters / iter->Dis2FCenters;
			if( iter->EDPScore < ratio15 / 2. )
				rvsLocations.push_back( *iter );
		}
		//First run of TIS model
		GeneSeq _gs_(geneInfo.positiveSeq);
		
		med_start= new MED_start( _gs_, rvsLocations );
		iterNum = 0;
		for( double iterationSTPTag = 1000000; ; )
		{
			med_start->setMotifRegions();
			med_start->findHitMotif1();
			rvsLocations = med_start->doIteration();

			double percent = fabs( iterationSTPTag /  med_start->getSTPTag() );
			if( ( percent > 0.99 && percent <1.01) || ++iterNum == 20 )
				break;
			else
				iterationSTPTag = med_start->getSTPTag();
			med_start->resetLearnedRst();
			med_start->setRvsLocation( rvsLocations );
		}
		cout<<endl;
		for( iter = trueORF.begin(); iter!=trueORF.end(); ++iter)
		{
			const char* seq = iter->isPositive ? geneInfo.positiveSeq.data() 
				:geneInfo.negtiveSeq.data();
			iter->RBSScore = med_start->reviseTIS( seq, *iter );
		}
		for( iter = ORFSet.begin(); iter!= ORFSet.end(); ++iter)
		{	
			const char* seq = iter->isPositive ? geneInfo.positiveSeq.data() 
				:geneInfo.negtiveSeq.data();
			iter->RBSScore = med_start->reviseTIS( seq, *iter );
		}
		delete med_start;

		std::vector<ORF_T> medResult;
		cout<<"Select coding ORFs from remianing ORFs" 
			" using a Fisher discriminant function..."<<endl;//
		for( iter = trueORF.begin(); iter!=trueORF.end(); ++iter)
		{
			iter->GC = GCContent( iter->isPositive ? geneInfo.positiveSeq.data() 
				:geneInfo.negtiveSeq.data(), iter->location.first, 
				iter->location.second );
			Str _st_(( iter->isPositive ? geneInfo.positiveSeq
                                : geneInfo.negtiveSeq ).substr( iter->location.first, iter->ORFLength ));
			iter->codingSeq = SequenceTransform_T::digital2CharSeq(_st_);
			medResult.push_back( *iter);
		}//*/
		for( iter = ORFSet.begin(); iter!= ORFSet.end(); ++iter)
		{
			//fisher descriminant between edp and atg
			if( isPassDecision( *iter ) )
			{
				iter->GC = GCContent( iter->isPositive ? geneInfo.positiveSeq.data() 
				:geneInfo.negtiveSeq.data(), iter->location.first, 
				iter->location.second );

				Str _st_(( iter->isPositive ? geneInfo.positiveSeq : geneInfo.negtiveSeq ).substr( iter->location.first, iter->ORFLength ));
				iter->codingSeq = SequenceTransform_T::digital2CharSeq(_st_);
				medResult.push_back( *iter);
			}
		}//*/
		int i = 0;
		for(; i < medResult.size(); ++i ){
			const char* seq = medResult[i].isPositive 
				? geneInfo.positiveSeq.data():geneInfo.negtiveSeq.data();
			x2LongestORF( seq, medResult[i].location );
			medResult[i].location2FileForm( seqLen );
		}
		cout<<"Remaining: "<<medResult.size()<<endl;
		return medResult;
	}
}


bool MED::isPassDecision( ORF_T& candidate )
{
	double RBS = candidate.RBSScore;
	double EDP = candidate.EDPScore;
	int length = candidate.ORFLength;

	if( EDP <  0.0640485 * (RBS+3.5)+ 3.61412 )	
	{ 	
		return true;
	}
	return false;
}

double MED::ORFDistanceRatio( std::list<ORF_T>::iterator& ORF )
{
	double* EDP = &(ORF->EDP[0]);
	int j = 0;
	double& TDistance = ORF->Dis2TCenters; 
	int i = TBeg;
	for(; i < TEnd; ++i )		
	{
		double distanceTmp = 0;
		double* trueCenter = &trueORFCenter[i][0];
		for( j = 0; j < 20; ++j )
			distanceTmp += fabs( EDP[j] - trueCenter[j] );
		if( distanceTmp < TDistance ){
			TDistance = distanceTmp;
			initTIndex = i;
		}
	}
	double& FDistance = ORF->Dis2FCenters;
	for( i = FBeg; i < FEnd; ++i)			
	{
		double distanceTmp = 0;
		double* falseCenter = &falseORFCenter[i][0];
		for( j = 0; j < 20; ++j )
			distanceTmp += fabs( EDP[j] - falseCenter[j] );
		if( distanceTmp < FDistance ){
			FDistance = distanceTmp;
			initFIndex = i;
		}
	}
	ORF->TDis = TDistance;
	ORF->FDis = FDistance;
	return TDistance / FDistance;
}

Ve_D MED::aminoEDPStrengthen ( int strengthenBegin, int startPos, int endPos )
{
	int endPosition = ( endPos + 1 < geneInfo.getSeqLength() ? endPos : seqLen - 1 ); 
	if( ( strengthenBegin - startPos ) % 3 != 0 || strengthenBegin >=  endPosition )
		return aminoEDPBase( startPos , endPosition );//here make efffects
	else
	{
		Ve_D init = aminoEDPBase( startPos , endPosition );
		Ve_D strengthen( 20, 0 );
		double sum = 0;
		formProbabilitySeq( strengthen, sum, strengthenBegin , endPosition );
		int i = 0;
		for(; i < 20; ++i )
		{
			strengthen[i] += init[i];
			sum += strengthen[i];
		}
		EDPKernl( strengthen, sum );
		return strengthen;
	}
}

Ve_D MED::aminoEDPBase( int startPos , int endPos )	
{
	Ve_D aminoFrequency( 20, 0 );
	double sum = 0;					
	formProbabilitySeq( aminoFrequency, sum, startPos, endPos );
	EDPKernl( aminoFrequency, sum );
	return aminoFrequency; 
}

void MED::EDPKernl( Ve_D& aminoFrequency, double sum )
{
	int j = 0;
	for(; j < 20; ++j )
	{
		if( aminoFrequency[j] != 0 )
		{
			double pi = aminoFrequency[j] / sum;
			aminoFrequency[j] = pi * log( pi );	
		}
	}

	sum = 0;
	for( j = 0; j < 20; ++j )
		sum += aminoFrequency[j];

	if( sum != 0 )
	{
		for( j = 0; j < 20; ++j)				//¾ùÒ»»¯¡£
			aminoFrequency[j] /= sum;
	}//*/
}

void MED::formProbabilitySeq
( Ve_D& aminoFrequency, double& sum, int startPos, int endPosition )
{
	assert( aminoFrequency.size() == 20 );
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )
	{
		++sum;
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 )
		{
		case 36 : case 37 : case 38 : case 39 :
			++aminoFrequency[0]; break;	//tmp += 'A'; break;
		case 57 : case 59 :
			++aminoFrequency[1]; break;	//tmp += 'C'; break;
		case 33 : case 35 :
			++aminoFrequency[2]; break;	//tmp += 'D'; break;
		case 32 : case 34 :
			++aminoFrequency[3]; break;	//tmp += 'E'; break;
		case 61 : case 63 :
			++aminoFrequency[4]; break;	//tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			++aminoFrequency[5]; break;	//tmp += 'G'; break;
		case 17 : case 19 :
			++aminoFrequency[6]; break;	//tmp += 'H'; break;
		case 12 : case 13 : case 15 : 
			++aminoFrequency[7]; break;	//tmp += 'I'; break;
		case 0 : case 2 : 
			++aminoFrequency[8]; break;	//tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			++aminoFrequency[9]; break;	//tmp += 'L'; break;
		case 14 : 
			++aminoFrequency[10]; break;//tmp += 'M'; break;
		case 1 : case 3 : 
			++aminoFrequency[11]; break;//tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			++aminoFrequency[12]; break;//tmp += 'P'; break;
		case 16 : case 18 : 
			++aminoFrequency[13]; break;//tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			++aminoFrequency[14]; break;//tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :  
			++aminoFrequency[15]; break;//tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 : 
			++aminoFrequency[16]; break;//tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			++aminoFrequency[17]; break;//tmp += 'V'; break;
		case 58 :
			++aminoFrequency[18]; break;//tmp += 'W'; break;
		case 49 : case 51 :
			++aminoFrequency[19]; break;//tmp += 'Y'; break;
		case 48 : case 50 : case 56 :
		{
			--sum;
			break;
		}
		default :;
		}
	}
}
double MED::getGB( const char* seq, int beg, int end ){
	double GB13 = 0, GB123 = 0;
	int index =0 ;
	int i = beg;
	for(; i < end  - 6; ++i, ++index ){
		if( seq[i] == '1' || seq[i] == '2' ) {
			if( index % 3 == 0 || index % 3 == 2 )
				++GB13;
			++GB123;
		}
	}
	if( GB123 == 0 ){
		int i = beg;
		for(; i < end; ++i )
			cout<<seq[i];
		cout<<endl;
		return -10000000;
	}
	else
		return 100*GB13/GB123;
}

void MED::updateGBDistri(const std::list<ORF_T>& trueORF){
	CodGB = GBDistri();
	NonCodGB = GBDistri();
	for( std::list<ORF_T>::const_iterator iter = trueORF.begin()
		; iter != trueORF.end(); ++iter ){
		const Location_T& rvsLocation = *iter;
		const char* seq = rvsLocation.isPositive ? geneInfo.positiveSeq.data()
			: geneInfo.negtiveSeq.data();
		if(rvsLocation.location.first>90){
			double gb = getGB( seq, rvsLocation.location.first - 90, rvsLocation.location.first );
			++NonCodGB.distri[ gb];
			++NonCodGB.number;
		}
		double gb = getGB( seq, rvsLocation.location.first,  rvsLocation.location.first + 90 );
		++CodGB.distri[ gb ];
		++CodGB.number;
	}
	std::map<int,double>::iterator it;
	for( it = NonCodGB.distri.begin()
		; it != NonCodGB.distri.end(); ++it ){
		it->second /= NonCodGB.number;
	}
	for(  it = CodGB.distri.begin(); it != CodGB.distri.end(); ++it ){
		it->second /= CodGB.number;
	}
}


double MED::getGBThreshold( std::list<ORF_T>& trueORF){
	const char* posSeq = geneInfo.positiveSeq.data()
		, * negSeq = geneInfo.negtiveSeq.data();
	Ve_D scores;
	std::list<ORF_T>::iterator iter;
	for( iter = trueORF.begin()
		; iter != trueORF.end(); ++iter ){
		double gb = getGB( iter->isPositive?posSeq:negSeq
				, iter->location.first, iter->location.second );
		iter->GBScore = gb;
		scores.push_back( gb );
	}
	double ave = GetAver(scores), dev = GetDerivate(scores);
	double maxScore = -10000;
	
	for(  iter = trueORF.begin(); iter != trueORF.end(); ++iter ){
		if( (iter->GBScore-ave)/dev < -3 ){
			if( iter->GBScore > maxScore )
				maxScore = iter->GBScore;
		}
	}	
	return maxScore;
}
