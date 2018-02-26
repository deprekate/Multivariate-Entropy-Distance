#ifndef MAXTRIX_H
#define MAXTRIX_H

#include<iomanip>
#include"TypeDefBase.h"
#include<vector>
#include <math.h>
				
template<class Type_T>
class Matrix_T 
{
public:
	Matrix_T(unsigned rows, unsigned cols,Type_T init=Type_T(0)):
	  line(rows),colum(cols),data(rows*cols, init ){};
	  Matrix_T(){}
	inline Type_T& operator()(unsigned i,unsigned j)							//operator () 
	{
		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	inline const Type_T operator()(unsigned i,unsigned j)const					//operator() 
	{			
		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	Matrix_T<Type_T>& toBeAveraged();
	Matrix_T<Type_T>& toBeLoged();
	Matrix_T<Type_T>& toBeExp();
	double getInfo(int cl1, int cl2, Ve_D& bg = Ve_D(4,1)){
		double info = 0;
		int col=cl1;
		for(;col<cl2;col++){
			int lin=0;
			for(;lin<line;lin++){
				if(data[ lin * colum + col ] != 0 ){
					info += data[ lin * colum + col ]
						*(log(data[ lin * colum + col ]/bg[lin])/log(2));
				}
			}
		}
		return info;
	}
	void matrixOut(std::ostream & out);
	void matrixIn(std::istream & in);
	inline const int getLine()const {return line;}	
	inline const int getColum()const {return colum;}
private:
	Matrix_T(std::vector<Type_T>& d,unsigned rows, unsigned cols)
		:line(rows),colum(cols),data(d){};
	std::vector<Type_T> data;
    unsigned line, colum;
};

template<class Type_T>
void Matrix_T<Type_T>::matrixOut(std::ostream & out)
{
	const Matrix_T<Type_T>& matrix = *this;
	out<<std::setw(15)<<matrix.getLine()<<std::setw(15)<<matrix.getColum()<<std::endl;
	int col=0;
	for(;col<matrix.getColum();++col)
	{
		int lin=0;
		for(;lin<matrix.getLine();++lin)
			out<<std::setw(15)<<std::setprecision(5)<<matrix(lin,col);
		out<<std::endl;
	}
}

template<class Type_T>
void Matrix_T<Type_T>::matrixIn(std::istream & in)
{
	Matrix_T<Type_T>& matrix = *this;
	in>>matrix.line>>matrix.colum;
	matrix.data.resize(matrix.line*matrix.colum);
	int col=0;
	for(;col<matrix.colum;++col){
		int lin=0;
		for(; lin<matrix.line;++lin)
			in>>matrix(lin,col);
	}
}

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeLoged()
{
	int col=0;
	for(;col<colum;col++)
	{
		int lin=0;
		for(;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = log(data[ lin * colum + col ]);
	}
	return *this;
}

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeExp()
{
	int col=0;
	for(;col<colum;col++)
	{
		int lin=0;
		for(;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = exp(data[ lin * colum + col ]);
	}
	return *this;
}

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeAveraged()									
{
	int col=0;
	for(;col<colum;++col)
	{
		Type_T sum=0;
		int lin=0;
		for(;lin<line;lin++)
			sum+=data[ lin * colum + col ];
		for(lin=0;lin<line;++lin)
			data[ lin * colum + col ]/=sum;
	}
	return *this;
}

#endif//difine MAXTRIX_H.



