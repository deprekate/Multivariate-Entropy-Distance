#ifndef MED_START_H
#define MED_START_H

#include"TypeDef.h"
#include"GeneInfo.h"
#include "GeneSeq.h"
#include"OftenUsedOperatLib.h"
enum geneState{ ISAGENE = 1, STPINFRAME, LENNOT3TIMES, TISERROR, STPERROR };
extern Str parafile;
extern int postmp;
class CParametersDlg
{
public: CParametersDlg(int atgPreNum = 4, int atgPosNum = 15)
	{
		m_atgPosNum = atgPosNum;
		m_atgPreNum = atgPreNum;
		m_ORFLen = 300;
		m_sigLen = 5;
		m_sigMutNum = 1;
		m_SDPosNum = 2;
		m_SDPreNum = 3;
		m_SDRegion = 20;
		m_TAPosNum = 2;
		m_TAPreNum = 3;
		m_TARegionEnd = 15;
		m_TARegion = 20;
	}
		int m_atgPosNum;
		int m_atgPreNum;
		int m_ORFLen;
		int m_sigLen;
		int m_sigMutNum;
		int m_SDPosNum;
		int m_SDPreNum;
		int m_SDRegion;
		int m_TAPosNum;
		int m_TAPreNum;
		int m_TARegion;
		int m_TARegionEnd;
};

void outMatrix(std::ostream & out,const M_D& matrix);

extern Str parafile;
class MED_start : public GeneInfo_T
{
public:
		class GBDistri{
		public:
			GBDistri(){
			}
		std::map<int,double> distri;
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
		void formGB();
	double getRvsLocations(int iterNum = 1);
	void reviseGivenORFs(Ve_Location &locations);
	void setCandidateMotifs( Ve_Str candidateMotifs );
	void setHitMotifs( Ve_Str hitMotifs );
	double getSTPTag()
	{
		return ATGWM(0, 0);
	}
	MED_start( GeneSeq& geneseq, Ve_Location& locations, 
		CParametersDlg parameters = CParametersDlg());
	Ve_Location doIteration();
	void setRvsLocation( Ve_Location& rvsLocations )
	{
		this->rvsLocations = rvsLocations;
	}
	void setParameters( CParametersDlg* m_parameters );
	void outReadableParameters( Str file, Str Class = " " );
	void resetLearnedRst(bool);
	void setRegions();
	Str m_motifsPath;
	void parameters2File( Str fileName );
	double reviseTIS( const char* seq, Location_T& ATGSTP );
	void readParametersFromFile( Str fileName );
protected:
	void formPreSTPScore(double aver);
	void formATGWM(double);
	void prepareForTesting(double aver);
	double reviseTIS( const char* seq, Location_T& ATGSTP, Ve_D& XTGIndex );
	double findMotifs( int mutBp, Ve_Str& region
		, int regionLen, Str file
		, int limitNum); 
	Pa_Ve_D_M_D signalDistriAndWM( const char* signal
		, Ve_Str& regions, int regionLen
		, int sigLen, int preBP,int posBP, int mutBp );
	Ve_Location test();
	bool hasSignal( const char* seq, const char* signal, int frameLen, int mutBp );
	bool isSignalInPos( const char* seq, const char* signal, int pos, int mutBp );
	Pa_D_D getSigScore( const char* seq, int pos, int preSigBP, int regionLen
		, int ATGPos, Ve_D& distri, M_D& sigWM );
	void findHitMotifs( int mutNum, Ve_Str& region
		, Ve_Str& candidateMotifs, Ve_Str& hitMotifs
		, int regionLen, Str file
		, bool (*filter)(Str&)= defaultSignal );
	void kerneWordsDistriAndWM( int mutNum );
	Ve_D signalDistri(  Str signal, int mutNum,int regionLen, Ve_Str& regions );
	void marskSigBK( Ve_M_D& sigWMs, int beg, int end =0);
	//void MED_start::marskSigBK( Ve_M_D& sigWMs, int beg, int end =0);
private:
	int sigLen;
	int SDRegionsLen, TARegionsLen, TARegionEnd;
	int learnORFLenLmt;
	int preSDBP, posSDBP, preTABP,posTABP,preATGBP, posATGBP;
	int mutBp;
	M_D ATGWM;
	Ve_D XTG;
	Ve_D XTGIndex;
	Ve_Str SDRegions, TARegions, ATGRegions;
	Ve_Loca rvsLocations;
	Ve_Str SDHitMotifs, TAHitMotifs;
	Ve_Str SDCandidateMotifs, TACandidateMotifs;
	Ve_M_D SDSigWM, TASigWM;
	Ve_Ve_D SDSpacerLenDistri,TASpacerLenDistri;
	double SDWeight, TAWeight, ATGWeight;
	double defaultPreSTPScore;
	std::map<int,double> toPreSTPScore;
};

#endif


