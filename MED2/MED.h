#ifndef MED_H
#define MED_H
#include "TypeDef.h"	// Added by ClassView
#include"MED_start.h"


class MED
{
public:
	MED( Str& seq );							
	std::vector<ORF_T> getGeneLocation();		
protected:

	const char* seq;							
	Ve_Location result;							
	GeneInfo_T geneInfo;
private:
	class GBDistri{
		public:
			GBDistri():number(0){
			}
			std::map<int,double> distri;
			double number;
			double getProba( int x ){		
				double p = 0;
				std::map<int, double>::iterator it = distri.begin();
				for( ; it != distri.end(); ++it ){
					if( it->first >= x )
						break;
				}
				if( it == distri.end() || it == distri.begin()){
					return distri.begin()->second;
				}
				else{
					std::map<int, double>::iterator posIt = it;
					--it;
					return (x - it->first ) / (posIt->first - it->first ) 
						* (posIt->second - it->second )+ it->second;
				}
			}
	};
	double getCodProd( int x ){
		double codingP = CodGB.getProba(x), noncodingP = NonCodGB.getProba(x);
		return codingP/(noncodingP+codingP);
	}
	GBDistri CodGB, NonCodGB;
	void decorateInitPoint( std::list<ORF_T>& ORFSet);		
	void getORFSet( std::list<ORF_T>& ORFSet );			
	std::vector<ORF_T> identifyGeneLocation					
		( 
		std::list<ORF_T>& ORFSet 
		);
	double ORFDistanceRatio( std::list<ORF_T>::iterator&F );	
	Ve_D aminoEDPStrengthen									
		(
		int strengthBegin, 
		int startPos = 0, 
		int endPos = int_limits::max() 
		);
	Ve_D aminoEDPBase( int startPos, int endPos );		
	void formProbabilitySeq								
		(
		Ve_D& probability,
		double& aminoSum,
		int startPos, 
		int endPosition
		); 
	void EDPKernl( Ve_D& aminoFrequency, double sum );	
	int TBeg, TEnd, FBeg, FEnd;								
	int seqLen;
	Ve_Ve_D trueORFCenter;								
	Ve_Ve_D falseORFCenter;									
	MED_start* med_start;									
	bool isPassDecision( ORF_T& candidate );
	M_D ATGWM;
	Ve_D index;
	void updateATGWMAndIndex(const std::list<ORF_T>& trueORF);
	double relocateTISByATGWManIndex( const char* seq, Location_T& iter );
	void relocateTISsByATGWManIndex( std::list<ORF_T>& ORFs );
	void refineATGWMAndIndex(std::list<ORF_T>& trueORF);
	void reEvaluateSeeds(std::list<ORF_T>& trueORF
		, std::list<ORF_T>& falseORF, std::list<ORF_T>& remain, double threshold );
	void updateGBDistri(const std::list<ORF_T>& trueORF);
	double getGB( const char* seq,  int beg, int end );
	double getGBThreshold( std::list<ORF_T>& trueORF);
	void readCenters(Str fileName );
};

#endif //GENEFINDBASE_H



