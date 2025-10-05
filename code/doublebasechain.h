#include "big.h"
#pragma once
const int Cmax = 7;
const int MatrixB = 300;
const int MatrixT = 200;
//const int MatrixQ = 112;


struct doubleBaseType
{
	int sign;
	int binaryExponent;
	int ternaryExponent;
};
typedef struct {
	int sign;//1 represents positive, -1 represents negative
	int  b;
	int  t;
}DBCChainItem;


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

int mici(int &n, int b, int t);

int cishu23(Big n, int &g, int &h);
//g=v_2(n),h=v_3(n)

int cishu23(Big n, int &g, int &h, int &x);
//g=v_2(n),h=v_3(n),x=2^g3^h




int gainbt(Big n, int &s, int &g, int &h);
//input gcd(n,6)=1

int gainmbNAF(Big n, int &s, int &g, int &h);
//input gcd(n,6)=1


int gaintree(Big n, int &s, int &g, int &h);
//input gcd(n,6)=1

int NAF(Big n,int * item, int &number);
//NAF

int wNAF(Big n, int * item, int &number, int Cmax=1);
//wNAF

int binary(Big n,int * item, int &number);
//NAF

int bt(Big n, doubleBaseType * item, int &number);
//binary ternary approach

int mbNAF(Big n, doubleBaseType * item, int &number);
//mbNAF

int wmbNAF(Big n, doubleBaseType * item, int &number, int Cmax=1);
int greedy(Big n, doubleBaseType *item, int &itemNumber, int Cmax = 1);



int tree(Big n, doubleBaseType * item, int &itemNumber, int Cmax = 1);

//int tree(Big n, doubleBaseType * item, int &number);
//tree approach

int ASRandom(Big n, doubleBaseType * item, int &number);
//Add/Sub-Random approach

//int triRandom(Big n, doubleBaseType * item, int &number);
//tri-Random approach
long long cpucycles(void);

Big pow(int b, int t, int q = 0);



Big pow1(int b, int t);



int Cmap(int c);
int CmapInverse(int c);
long long power(int b, int t, int q = 0);
double value(int bmax, int tmax, int length, double C1 = 1, double C2 = 0, double C3 = 0);

int positive(int i);
Big rand(int n);
int weishu(Big n);

DBC DynamicProgrammingDBC(Big n, DBC Item[MatrixB][MatrixT][2 * Cmax],
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles);

//DBC DynamicProgrammingDBC0(Big n, DBC Item[MatrixB][MatrixT][2],
//	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles);

//DBC DynamicProgrammingDBC0(Big n, DBCChainItem finalChain[500], double &costFinal, double C1, double C2, double C3);

DBC DynamicProgrammingDBC0(Big n, DBC Item[MatrixB][MatrixT][2],
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles);