#ifndef MAXTRIX_H
#define MAXTRIX_H

#include <iomanip>
#include "TypeDefBase.h"
#include <vector>
#include <math.h>
#include <execinfo.h>

template<class Type_T> class Matrix_T;
template<class Type_T> Matrix_T<Type_T> operator-(Matrix_T<Type_T>& lMatrix,Matrix_T<Type_T>& rMatrix);	
				
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
	void rand();
	Matrix_T<Type_T>& toBeAveraged();
	Matrix_T<Type_T>& toBeLogo();
	Matrix_T<Type_T>& toBeLoged();
		Matrix_T<Type_T>& toExp();
	Matrix_T<Type_T>& append(const Matrix_T<Type_T>& m);	//append a matrix to a exist one.
	Matrix_T<Type_T>& operator+=(Matrix_T<Type_T>& matrix);
	void matrixOut(std::ostream & out);
	void matrixIn(std::istream & in);
	friend Matrix_T<Type_T> operator- <>(Matrix_T<Type_T>& lMatrix,Matrix_T<Type_T>& rMatrix);	
	//friend Matrix_T<Type_T> operator-(Matrix_T<Type_T>& lMatrix,Matrix_T<Type_T>& rMatrix);	
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
		for(;lin<matrix.line;++lin)
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
Matrix_T<Type_T>& Matrix_T<Type_T>::toExp(){
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
//	toBeLogo();
	return *this;
}
template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeLogo(){
	int col=0;
	for(;col<colum;++col){
		Type_T sum=log(line)/log(2);
		int lin=0;
		for(;lin<line;lin++)
			if( (*this)(lin, col ) != 0 )
				sum+=(*this)(lin, col )*log((*this)(lin, col ))/log(2);//pi*log(pi)
		for(lin=0;lin<line;++lin)
			(*this)(lin, col ) *= sum;
	}
	return *this;
}
template<class Type_T>
Matrix_T<Type_T> operator-(Matrix_T<Type_T>& lMatrix,Matrix_T<Type_T>& rMatrix)
{
	assert(lMatrix.getColum()==rMatrix.getColum()&&lMatrix.getLine()==rMatrix.getLine());
	Matrix_T<Type_T> tmp( rMatrix.getLine(), rMatrix.getColum() );
	int i = 0;
	for(; i < rMatrix.getLine(); ++i )
	{
		int j = 0;
		for(; j < rMatrix.getColum(); ++j )
			tmp( i, j ) = lMatrix( i, j ) - rMatrix( i, j );
	}
	return tmp;
}

template<class Type_T>													
Matrix_T<Type_T>& Matrix_T<Type_T>::append(const Matrix_T<Type_T>& m)
{
	assert(this!=&m&&m.line==line);
	Matrix_T<Type_T> temp(line,colum+m.colum);
	int lin = 0;
	for(; lin < line; ++lin )
	{
		int col = 0;
		for(; col < colum; ++col )
			temp( lin, col ) = data[ lin * colum + col ];
		for( ; col < colum+m.colum; ++col )
			temp( lin, col ) = m( lin, col - colum );
	}
	return *(new Matrix_T<Type_T>(temp));
}

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::operator+=(Matrix_T<Type_T>& matrix)
{
	assert(colum==matrix.getColum()&&line==matrix.getLine());
	data+=matrix.data;
	return *this;
}

template<class Type_T>
void Matrix_T<Type_T>::rand()
{
	int i = 0;
	for(; i < data.size(); ++i )
		data[i] = random();
		//data[i] = Random();
}

#endif//difine MAXTRIX_H.



