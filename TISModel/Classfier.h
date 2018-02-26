#include"TypeDef.h"
using namespace std;

struct Gene
{
	int beg;
	int end;
	char strand;
	int length;
	int PID;
	string name;
	string info;
};
int lim = 50, lim2 = 50;
#define ORFLIM 300
std::map<Str,Ve_Location>  classGene( Ve_Location ls, Ve_D bg ){
	map<pair<int,int>, Location_T> lss;
	int k = 0;
	for(; k < ls.size(); ++k ){
		lss[ls[k].location] = ls[k];
	}
	vector<pair<int,int> > PGene;
	vector<pair<int,int> > NGene;
	int i = 0;
	for(; i < ls.size(); ++i ){
		pair<int,int> gene = ls[i].location;
		{
			if( !ls[i].isPositive ){
				NGene.push_back( gene );
			}
			else{
				PGene.push_back( gene );
			}
		}
	}
	sort( PGene.begin(), PGene.end() );
	sort( NGene.begin(), NGene.end() );
	
	std::map<Str,Ve_Location> rst;

	double SDSeedNum = 0, OtSeedNum = 0; 
	for( i = 1; i < PGene.size() - 1; ++i ){
		if(PGene[i].first - PGene[i-1].second < lim 
			&& PGene[i].first - PGene[i-1].second > -lim){
			map<pair<int,int>, Location_T>::iterator 
				it = lss.find(PGene[i]);
			if( it != lss.end() && it->second.ORFLength > ORFLIM )
				rst["TUI"].push_back( it->second );
		}
		else
		{
			map<pair<int,int>, Location_T>::iterator 
				it = lss.find(PGene[i]);
			if( it != lss.end()&& it->second.ORFLength > ORFLIM ){
				rst["TUL"].push_back( it->second );
			}
		}
	}
	for( i = 0; i < NGene.size() - 1; ++i ){
		if( NGene[i+1].first - NGene[i].second < lim 
			&& NGene[i+1].first - NGene[i].second > -lim ){
			map<pair<int,int>, Location_T>::iterator 
				it = lss.find(NGene[i]);
			if( it != lss.end()&& it->second.ORFLength > ORFLIM )
				rst["TUI"].push_back( it->second );
		}
		else
		{
			map<pair<int,int>, Location_T>::iterator 
				it = lss.find(NGene[i]);
			if( it != lss.end()&& it->second.ORFLength > ORFLIM )
				rst["TUL"].push_back( it->second );
		}
	}//*/
	std::map<Str,Ve_Location>::iterator it = rst.begin();
	for( ; it != rst.end(); ++it ){
		cout<<it->first<<" genes number: "
			<<it->second.size()<<endl;
	}
	return rst;
}


