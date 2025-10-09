

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>

#include "preComputation.h"
#include "scalarMultiplication.h"
#include "Supersingularity.h"
using namespace std;


miracl *mip;

const int infinite = 100000;
const int bLength = 770;
const int tLength = 486;
//
//long long cpucycles(void)
//{
//	return __rdtsc();
//}

int fOut(int item[], int length)
{
	int PINAFLength = 0;
	for (int i = length - 1; i > -1; i--)
	{
		if (item[i] == 1)
		{
			cout << "+2^{" << i << "}";
			PINAFLength++;
		}
		else if (item[i] == -1)
		{
			cout << "-2^{" << i << "}";
			PINAFLength++;
		}
		else if (item[i] > 0)
		{
			cout << "+" << item[i] << "*2^{" << i << "}";
			PINAFLength++;
		}
		else if (item[i] < 0)
		{
			cout << item[i] << "*2^{" << i << "}";
			PINAFLength++;
		}
	}
	cout << "Hamming weight:" << PINAFLength << endl;
	return 0;
}
int fOut(doubleBaseType item[], int length)
{

	for (int i = length - 1; i > -1; i--)
	{
		if (item[i].sign == 1)
		{
			cout << "+2^{" << item[i].binaryExponent << "}" << "3^{" << item[i].ternaryExponent << "}";

		}
		else if (item[i].sign == -1)
		{
			cout << "-2^{" << item[i].binaryExponent << "}" << "3^{" << item[i].ternaryExponent << "}";

		}
		else if (item[i].sign > 0)
		{
			cout << "+" << item[i].sign << "*2^{" << item[i].binaryExponent << "}" << "3^{" << item[i].ternaryExponent << "}";

		}
		else if (item[i].sign < 0)
		{
			cout << item[i].sign << "*2^{" << item[i].binaryExponent << "}" << "3^{" << item[i].ternaryExponent << "}";
		}
	}

	return 0;
}


int assignment(DBCChainItem *finalChain, int DBCLength, doubleBaseType *DBCTerm)
{
	for (int i = 0; i < DBCLength; i++)
	{
		DBCTerm[i].sign = finalChain[i].sign;
		DBCTerm[i].binaryExponent = finalChain[i].b;
		DBCTerm[i].ternaryExponent = finalChain[i].t;

	}
	return 0;
}



int main()
{
	int l[74] = {179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,      
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,		 
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,		 
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3};
	//srand(time(NULL));
	mip = mirsys(1024, 2);
	mip->IOBASE = 16;
	const  int Cycle = 10;

	int start, end;
	Big n = "3141592653589793238462643383279502884197169399377123340155612345";
	Big m = n;


	//Big n = 1;

	n = 1;
	for (int j = 33; j < 74; j++)
	{
		n = n*l[j];
	}
	n = n * 4;
	//n = n - 1;
	//n = n + 1;
	//n = n / 4;
	//mip->IOBASE = 10;
	cout << n<<":" << endl;

	//system("Pause");
	//cout<<"Bits: " <<bits(n)<<endl;
	Big ECC256("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF");


	Big ECC25519("7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED");
	//mip->IOBASE =10;
	Big out = 0;
	big p = mirvar(0);
	copy(ECC25519.getbig(), p);

	big x = mirvar(1);
	big y = mirvar(1);
	big z = mirvar(1);
	big d = mirvar(1);
	//copy(n.getbig(),x);
	//copy(n.getbig(),y);
	projectPoint P[65];//P[i]=iP for i\in I_w
	big X[65];
	big Y[65];
	big Z[65];
	//big T[65];
	Big initialX, initialY;
	//srand(time(NULL));
	for (int i = 0; i < 65; i++)
	{
		X[i] = mirvar(0);
		P[i].X = X[i];
		Y[i] = mirvar(0);
		P[i].Y = Y[i];
		Z[i] = mirvar(1);
		P[i].Z = Z[i];
		//T[i] = mirvar(1);
		//P[i].X3 = T[i];
	}



	int item[810];
	int length;
	//n=12312451;
	binary(n, item, length);

	prepare_monty(p);
	affinePoint P2;
	projectPoint P1, P3;
	P2.x = mirvar(2);
	P2.y = mirvar(5);
	affinePoint P4;
	P4.x = mirvar(2);
	P4.y = mirvar(5);

	nres(P2.x, P2.x);
	nres(P2.y, P2.y);

	big Z1 = mirvar(1);

	big mid = mirvar(0);
	big one = mirvar(1);
	mip->IOBASE = 10;
	int PINAFLength = 0;
 
	long long CPUcycles;
	int IsSupersingularity;
	big A = mirvar(2);
	//productTree(A, IsSupersingularity);



	doubleBaseType DBCTerm[300];
	int DBCLength;

	DBC Item[MatrixB][MatrixT][2 * Cmax];
	DBC Item0[MatrixB][MatrixT][2];
	DBCChainItem finalChain[MatrixB];
	DBC DBCfinal;
	double mAcost = 6;
	double Acost = 7;
	double SMRatio = 0.67;
	double Dcost = 3 + 4 * SMRatio;
	double Tcost = 9 + 3 * SMRatio;
	//int IsSupersingularity = 0;

	//newProductTree(n.getbig(), IsSupersingularity);

	n = 1;

	for (int j = 32; j < 74; j++)
	{
		n = n*l[j];
	}
	n = n * 4;
	DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, mAcost, Acost, Dcost, Tcost, CPUcycles);

	int DBClength=DBCfinal.HammingWeight;

	assignment(finalChain, DBCLength, DBCTerm);
	newProductTree(A, DBCTerm, DBClength,IsSupersingularity);
	//pairingMontgomeryLadder(A, IsSupersingularity);
	if (IsSupersingularity < 1)
	{
		cout << "NonSupersingular!" << endl;
	}
	else
	{
		cout << "Supersingular!" << endl;
	}

	
	return 0;
}