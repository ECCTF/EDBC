/*
 *   Example program demonstrates 1024 bit Diffie-Hellman, El Gamal and RSA
 *   and 168 bit Elliptic Curve Diffie-Hellman
 *
 *   Requires: big.cpp ecn.cpp
 */

#include <iostream> 
#include "big.h"
#include <ctime>
#include "AdaptiveDBC.h"
 //#include "stdio.h"
#include "intrin.h"
#include "string.h"
#include <cmath>
using namespace std;

/* large 1024 bit prime p for which (p-1)/2 is also prime */
char *primetext = (char *)
"155315526351482395991155996351231807220169644828378937433223838972232518351958838087073321845624756550146945246003790108045940383194773439496051917019892370102341378990113959561895891019716873290512815434724157588460613638202017020672756091067223336194394910765309830876066246480156617492164140095427773547319";

/* NIST p192 bit elliptic curve prime 2#192-2#64-1 */

char *ecp = (char *)"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF";

/* elliptic curve parameter B */

char *ecb = (char *)"64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1";

/* elliptic curve - point of prime order (x,y) */

char *ecx = (char *)"188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
char *ecy = (char *)"07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

char *text = (char *)"MIRACL - Best multi-precision library in the World!\n";

#ifndef MR_NOFULLWIDTH
Miracl precision(50, 0);
#else 
Miracl precision(50, MAXBASE);
#endif

// If MR_STATIC is defined in mirdef.h, it is assumed to be 100

//Miracl precision(120,(1<<26));



int initialPreComputation(double *M, double *S)
{
	M[1] = 0;
	M[3] = 9;
	M[5] = 19;
	M[7] = 27;
	M[9] = 35;
	M[11] = 43;
	M[13] = 51;
	M[15] = 59, M[17] = 67, M[19] = 67, M[21] = 75;
	M[23] = 83, M[25] = 90, M[27] = 96, M[29] = 104;
	M[31] = 109;
	S[1] = 0;
	S[3] = 2;
	S[5] = 3;
	S[7] = 3;
	S[9] = 3;
	S[11] = 3;
	S[13] = 3;
	S[15] = 3, S[17] = 3, S[19] = 10, S[21] = 10;
	S[23] = 10, S[25] = 11, S[27] = 14, S[29] = 14;
	S[31] = 14;
	return 0;
}
int main()
{
	double M[32];
	double S[32];
	initialPreComputation(M, S);
	const double Infinite = 100000000;
	Big n = 3;
	//n = 31415926535897932;
	n = rand(256);
	DBC Item[MatrixB][MatrixT][2 * Cmax];
	DBC DBCfinal;
	/*Big pri[9] = { 29,61,67,71,97,127,131,227,229 };*/
	//n = 1;
	//for (int i = 0; i < 7; i++)
	//{
	//	for (int j = i + 1; j < 8; j++)
	//	{
	//		for (int k = j + 1; k < 9; k++)
	//		{
	//			n = pri[i] * pri[j] * pri[k];
	//		}
	//	}
	//}

	/*n = 29 * 61 * 67 * 127;
	n = n * 97 * 131;*/
	/*n = n * 71;
	n = n * 227 * 229;*/
	double costFinal = Infinite;
	double mA = 6;
	double A = 7;
	double SMRatio = 0.67;
	double D = 3 + 4 * SMRatio;
	double T = 9 + 3 * SMRatio;
	double Q = 0;
	DBCChainItem finalChain[MatrixB];
	long long CPUcycles = 0;
	long long CPUcycles1 = 0;


	//const int CYCLES = 200;
	//for (int i = 0; i < CYCLES; i++)
	//{
	//	cout << "Test for 256-bit integer to generate optimal EDBC." << endl;

	//	srand(time(NULL));
	//	n = rand(256);
	//	//cout << n << endl;
	//	DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, mA, A, D, T, CPUcycles);
	//	cout << "Time cost of our Dynamic algorithm:" << CPUcycles << endl;
	//	DBCfinal = DynamicProgrammingDBCBernstein(n, Item, finalChain, mA, A, D, T, CPUcycles1);
	//	cout << "Time cost of Bernstein's algorithm:" << CPUcycles1 << endl;

	//	double ratio = (double)CPUcycles1 / (double)CPUcycles;
	//	cout << "The ratio of time cost of our Dynamic Programming to Bernstein's algorithm is 1/" << CPUcycles1 / CPUcycles << endl;

	//}


	n = rand(256);
	cout << "Please input a positive integer(less than 256 bit)." << endl;
	while (cin >> n)
	{
		if (bits(n) > 256)
		{
			cout << "Please input a positive integer(less than 256 bit)." << endl;
			continue;
		}
		DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, mA, A, D, T, CPUcycles);
		cout << "An optimal EDBC with coefficient set $C^{\ast}=\{\pm1,\pm 3,\pm 5,\pm 7\} on Edwards curve E_{Ex} for a  256-bit integer" << endl;
		cout << n << " is " << endl;
		int finalChainLength = DBCfinal.HammingWeight;
		for (int i = finalChainLength - 1; i > -1; i--)
		{
			if (i == finalChainLength - 1)
			{
				if (finalChain[i].sign == 1)
				{
					//cout << "";
				}
				else if (finalChain[i].sign == 3)
					cout << "3*";
				else if (finalChain[i].sign == 5)
					cout << "5*";
				else if (finalChain[i].sign == 7)
					cout << "7*";
				else if (finalChain[i].sign == -1)
					cout << "-";
				else if (finalChain[i].sign == -3)
					cout << "-3*";
				else if (finalChain[i].sign == -5)
					cout << "-5*";
				else if (finalChain[i].sign == -7)
					cout << "-7*";
			}
			else
			{
				if (finalChain[i].sign == 1)
					cout << "+";
				else if (finalChain[i].sign == 3)
					cout << "+3*";
				else if (finalChain[i].sign == 5)
					cout << "+5*";
				else if (finalChain[i].sign == 7)
					cout << "+7*";
				else if (finalChain[i].sign == -1)
					cout << "-";
				else if (finalChain[i].sign == -3)
					cout << "-3*";
				else if (finalChain[i].sign == -5)
					cout << "-5*";
				else if (finalChain[i].sign == -7)
					cout << "-7*";
			}
			cout << "2^{" << finalChain[i].b << "}3^{" << finalChain[i].t << "}";
		}
		cout << endl;
		cout << "Please input a positive integer(less than 256 bit)." << endl;
	}



	return 0;
}



