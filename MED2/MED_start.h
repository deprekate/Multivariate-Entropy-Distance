#ifndef MED_START_H
#define MED_START_H

#include"TypeDef.h"
#include"GeneInfo.h"
#include "GeneSeq.h"
enum geneState{ ISAGENE = 1, STPINFRAME, LENNOT3TIMES, TISERROR, STPERROR };
extern Str parafile;
class CParametersDlg
{
public: CParametersDlg()
	{
		m_atgPosNum = 15;
		m_atgPreNum = 4;
		m_ORFLen = 300;
		m_sigLen = 5;
		m_sigMutNum = 1;
		m_sigPosNum = 2;
		m_sigPreNum = 3;
		m_sigRegion = 20;
	}
		int m_atgPosNum;
		int m_atgPreNum;
		int m_ORFLen;
		int m_sigLen;
		int m_sigMutNum;
		int m_sigPosNum;
		int m_sigPreNum;
		int m_sigRegion;
};


void x2LongestORF( const char* seq,  Pa_I_I& location );
class MED_start : public GeneInfo_T
{
public:
	void saveParameter( std::ofstream& out );
	Ve_Location getRvsLocations();
	void findHitMotif1();
	Ve_Str getCandidateMotifs();
	void setCandidateMotifs( Ve_Str candidateMotifs );
	Ve_Str getHitMotifs();
	Ve_Ve_D getCandidateMotifSpacerDistri()
	{
		return candidateSpacerLenDistri;
	}
	void setHitmotifs( Ve_Str hitMotifs );
	double getSTPTag()
	{
		return ATGWM(0, 0);
	}
	MED_start( GeneSeq& geneseq, Ve_Location& locations, 
		CParametersDlg parameters = CParametersDlg());
	Ve_Location doIteration();
	M_D getATGWM();
	void setRvsLocation( Ve_Location& rvsLocations )
	{
		this->rvsLocations = rvsLocations;
	}
	void setParameters( CParametersDlg* m_parameters );
	void resetLearnedRst();
	void resetLearnedRst(bool);
	void setMotifRegions();
	Ve_Location getRst()
	{
		return rvsLocations;
	}
	Str m_motifsPath;
	Ve_D getXTG()
	{
		return XTGTmp;
	}
	Ve_D getXTGIndex()
	{
		return XTGIndexTmp1;
	}
	double reviseTIS( const char* seq, Location_T& ATGSTP );
protected:
	void formATGWM();
	void prepareForTesting();
	double reviseTIS( const char* seq, Location_T& ATGSTP, Ve_D& XTGIndex );
	void findMotifs( int mutBp ); 
	Pa_Ve_D_M_D signalDistriAndWM( const char* signal, int sigLen, int mutBp );
	Ve_Location test();
	bool hasSignal( const char* seq, const char* signal, int frameLen, int mutBp );
	bool isSignalInPos( const char* seq, const char* signal, int pos, int mutBp );
	Pa_D_D getRBSScore( const char* seq, int pos, int ATGPos, Ve_D& distri, M_D& sigWM );
	void findHitMotifs( int mutNum );
	void kerneWordsDistriAndWM( int mutNum );
	int geneStatus( Str& seq );
	Ve_D signalDistri(  Str signal, int mutNum );
	int theIthATGInORF( Str& seq, int ATGPos, int ORFBeg, int ORFEnd );
private:
	int sigLen;
	int motifRegionsLen;
	int learnORFLenLmt;
	int preSigBP, posSigBP, preATGBP, posATGBP;
	int mutBp;
	M_D ATGWM;
	Ve_D XTG, XTGTmp;
	Ve_D XTGIndex, XTGIndexTmp1;
	Ve_Str motifRegions;
	Ve_Loca rvsLocations;
	Ve_Str hitMotifs;
	Ve_Str candidateMotifs;
	Ve_M_D SigWM;
	Ve_Ve_D spacerLenDistri, spacerLenDistriTmp, candidateSpacerLenDistri;
	int containedATGNum( const char* seq, int beg, int end );
};

#endif


