#ifndef ADAPTIVEDBC_H
#define ADAPTIVEDBC_H
#endif

#include <iostream>
#include "ecn.h"
#include "big.h"
#include <ctime>

long long cpucycles(void);

const int Cmax = 7;
const int MatrixB = 260;
const int MatrixT = 200;
 
//const int MatrixQ = 112;

typedef struct {
	int sign;//1 represents positive, -1 represents negative
	int  b;
	int  t;
}DBCChainItem;

typedef struct {
	int sign;//1 represents positive, -1 represents negative
	int  b;
	int  t;
	int q;
}TBCChainItem;

typedef struct {
public:
	int sign;//1 represents positive, -1 represents negative
	int  b;
	int  t;
	int HammingWeight;
	double additionCost;
	int coefficient;
	int length;
	int flag;
	//Big n;
}DBC;



Big pow(int b, int t, int q=0);



Big pow1(int b, int t);



int Cmap(int c);
int CmapInverse(int c);
long long power(int b, int t, int q=0);
double value(int bmax, int tmax, int length, double C1 = 1, double C2 = 0, double C3 = 0);

int positive(int i);
Big rand(int n);
int weishu(Big n);


DBC DynamicProgrammingDBC(Big n, DBC Item[MatrixB][MatrixT][2* Cmax], 
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles);


DBC DynamicProgrammingDBCBernstein(Big n, DBC Item[MatrixB][MatrixT][2 * Cmax], 
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles);

