#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>

#include "preComputation.h"
#include "scalarMultiplication.h"

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

	//srand(time(NULL));
	mip = mirsys(1024, 2);
	mip->IOBASE = 16;
	const  int Cycle = 10;

	int start, end;
	Big n = "3141592653589793238462643383279502884197169399377123340155612345";
	Big m = n;

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

	//n = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165";


	//n = "31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867";
	//cout << "位数:"<< bits(n) << endl;

	//n = 4;
	//NAF(n, item, length);
	//cout << "位数:" << length << endl;
	/*for (int i = length-1; i >-1; i--)
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

	}
	cout << "NAF Length;" << PINAFLength << endl;*/

	P1.X = mirvar(2);
	P1.Y = mirvar(5);
	P1.Z = mirvar(1);



	////Big D="-47";
	////D=D+ECC256;
	////copy(D.getbig(),d);
	////nres(d,d);
	////cout <<D<<endl;
	////out=P1.X;
	////cout <<out<<":";
	////out=P1.Y;
	////cout << out <<":";
	////out=P1.Z;
	////cout <<out<<endl;

	nres(P1.X, P1.X);
	nres(P1.Y, P1.Y);
	nres(P1.Z, P1.Z);

	P3.X = mirvar(1);
	P3.Y = mirvar(1);
	P3.Z = mirvar(1);



	////out=P1.X;
	////cout <<out<<":";
	////out=P1.Y;
	////cout << out <<":";
	////out=P1.Z;
	////cout <<out<<endl;
	//JacobiMixedAddition(P1,P2,P3);

	//InvEdwardsMixedAddition(P1,P2,P3,d);
	//JacobiTriple(P1, P3);
	//InvEdwardsDouble(P1,P3,d);
	//JacobiTriple(P1,P3);
	//	InvEdwardsTriple(P1,P3,d);


	extProjectPoint PD[65];//P[i]=iP for i\in I_w
	projectPoint Pre[65];
	big T[65];
	//Big initialX, initialY;
	//srand(time(NULL));
	for (int i = 0; i < 65; i++)
	{
		X[i] = mirvar(0);
		PD[i].X = X[i];
		Y[i] = mirvar(0);
		PD[i].Y = Y[i];
		Z[i] = mirvar(1);
		PD[i].Z = Z[i];
		T[i] = mirvar(1);
		PD[i].T = T[i];
	}

	big X1[65];
	big Y1[65];
	big Z2[65];
	for (int i = 0; i < 65; i++)
	{
		X1[i] = mirvar(0);
		Pre[i].X = X1[i];
		Y1[i] = mirvar(0);
		Pre[i].Y = Y1[i];
		Z2[i] = mirvar(1);
		Pre[i].Z = Z2[i];
	}
	long long CPUcycles;
	const int Cmax1 = 31;
	big a = mirvar(1);
	big one1 = mirvar(1);
	//preComputationW(P2, P, CPUcycles, Cmax);
	//preComputationE(P2, P, CPUcycles, Cmax, one1, a);
	//preComputationD(P2,PD, CPUcycles, Cmax);
	//preComputationEx(P2, PD, CPUcycles, Cmax1, one1);

	Big Px;
	Big Py;
	Big Pz;
	Big Pt;
	Px = PD[31].X;

	Py = PD[31].Y;
	Pz = PD[31].Z;
	Pt = PD[31].T;
	//cout << "X:" << Px << " Y:" << Py << " Z:" << Pz << " T:" << Pt << endl;






	doubleBaseType DBCTerm[300];
	int DBCLength;

	DBC Item[MatrixB][MatrixT][2 * Cmax];
	DBC Item0[MatrixB][MatrixT][2];
	DBCChainItem finalChain[MatrixB];
	DBC DBCfinal;


	preComputationEx(P2, PD, CPUcycles, Cmax1, d);


	double TimeCost[12];
	for (int i = 0; i < 12; i++)
	{
		TimeCost[i] = 0;
	}
	//srand(time(NULL));
	n = rand(256);

	double mAcost = 6;
	double Acost = 7;
	double SMRatio = 0.67;
	double Dcost = 3 + 4 * SMRatio;
	double Tcost = 9 + 3 * SMRatio;

	



	for (int i = 0; i < Cycle; i++)
	{
		//n = rand(256);
		wNAF(n, item, length, 1);
		// preComputationEx(P2, PD, CPUcycles, 1, d);
		start = cpucycles();
		windowScalarMultiplicationEx(P2, P4, item, length, PD, 1);
		// scalarMultiplicationAd(P2, P4, item, length);
		end = cpucycles();
		TimeCost[0] += (end - start);//NAF time


	   // wNAF(n, item, length, 15);
	   // start = cpucycles();
	   // preComputationEx(P2, PD, CPUcycles, 15, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, 15);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[1] += (end - start);//wNAF time


	   // //wmbNAF
	   //wmbNAF(n, DBCTerm, DBCLength, 1);
	   // start = cpucycles();
	   // //preComputationEx(P2, PD, CPUcycles, 1, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, 1);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[2] += (end - start) ;//mbNAF time
	//

	   // int wmbNAFCmax = 7;
	   // wmbNAF(n, DBCTerm, DBCLength, wmbNAFCmax);
	   // start = cpucycles();
	   // preComputationEx(P2, PD, CPUcycles, wmbNAFCmax, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, wmbNAFCmax);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[3] += (end - start);//wmbNAF time


	   // //greedy
	   // greedy(n, DBCTerm, DBCLength, 1);
	   // start = cpucycles();
	   // //preComputationEx(P2, PD, CPUcycles, Cmax1, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, 1);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[4]+= (end - start) ;//greedy time

	   // int greedyCmax = 5;
	   // greedy(n, DBCTerm, DBCLength, greedyCmax);
	   // start = cpucycles();
	   // preComputationEx(P2, PD, CPUcycles, greedyCmax, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, greedyCmax);
	   // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[5] += (end - start) ;//window greedy time


	   // //tree
	   // tree(n, DBCTerm, DBCLength, 1);
	   // start = cpucycles();
	   //// preComputationEx(P2, PD, CPUcycles, 1, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, 1);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[6] += (end - start) ;//tree time
	//
	   // int treeCmax = 11;
	   // tree(n, DBCTerm, DBCLength, treeCmax);
	   // start = cpucycles();
	   // preComputationEx(P2, PD, CPUcycles, treeCmax, d);
	   // windowScalarMultiplicationEx(P2, P4, item, length, PD, treeCmax);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[7] += (end - start);//window tree time




	//  //canonic EDBC
	   // double costFinal = 0;
	   ////DBCfinal = DynamicProgrammingDBC0(n, finalChain, costFinal, 1, 0, 0);
	   // DBCfinal = DynamicProgrammingDBC0(n, Item0, finalChain, 1, 1, 0, 0, CPUcycles);
	   // DBCLength = DBCfinal.HammingWeight;

	   // assignment(finalChain, DBCLength, DBCTerm);
	   // start = cpucycles();		 
	   // windowScalarMultiplicationEx(P2, P4, DBCTerm, DBCLength, PD, 1);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[8] += (end - start);//window tree time
	   //							  /*fOut(DBCTerm, DBCLength);
	   //							  cout << DBCfinal.HammingWeight << endl;*/


	   //							  //int canonicCmax = Cmax;

	   // int canonicCmax = Cmax;
	   //DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, 1, 1, 0, 0, CPUcycles);
	   // DBCLength = DBCfinal.HammingWeight;

	   // assignment(finalChain, DBCLength, DBCTerm);
	   // start = cpucycles();
	   // preComputationEx(P2, PD, CPUcycles, canonicCmax, d);
	   // windowScalarMultiplicationEx(P2, P4, DBCTerm, DBCLength, PD, canonicCmax);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[9] += (end - start);//window tree time
	   //							  /*fOut(DBCTerm, DBCLength);
	   //							  cout << DBCfinal.HammingWeight << endl;*/



	   // //optimal EDBC
	   ////DBCfinal = DynamicProgrammingDBC0(n, finalChain, costFinal, 7, 7, 12);
	   //DBCfinal = DynamicProgrammingDBC0(n, Item0, finalChain, 7, 7, 7, 12, CPUcycles);
	   // DBCLength = DBCfinal.HammingWeight;
	   //
	   // assignment(finalChain, DBCLength, DBCTerm);
	   // start = cpucycles();		
	   // windowScalarMultiplicationEx(P2, P4, DBCTerm, DBCLength, PD, 1);
	   // // scalarMultiplicationAd(P2, P4, item, length);
	   // end = cpucycles();
	   // TimeCost[10] += (end - start);//window tree time
	   // /*fOut(DBCTerm, DBCLength);
	   // cout << DBCfinal.HammingWeight << endl;*/


		int optimalCmax = Cmax;

		//DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, mAcost, Acost, Dcost, Tcost, CPUcycles);
		DBCfinal = DynamicProgrammingDBC(n, Item, finalChain, 7, 8, 7, 12, CPUcycles);
		
		DBCLength = DBCfinal.HammingWeight;

		assignment(finalChain, DBCLength, DBCTerm);
		preComputationEx(P2, PD, CPUcycles, optimalCmax, d);
		start = cpucycles();
		
		windowScalarMultiplicationEx(P2, P4, DBCTerm, DBCLength, PD, optimalCmax);
		// scalarMultiplicationAd(P2, P4, item, length);
		end = cpucycles();
		TimeCost[11] += (end - start);//window tree time
									 /*fOut(DBCTerm, DBCLength);
									 cout << DBCfinal.HammingWeight << endl;*/




	}
	 
	cout << "This test uses a combination of projective coordinates and extended coordinates on Edwards curve with $a=-1$." << endl;
	cout << "The ratio of time cost of scalar multiplication using our optimal EDBC and the time cost of scalar multiplication  using NAF is  "
		<< TimeCost[11] / TimeCost[0] << endl;

	cout << "Compared to scalar multiplication using NAF, our algorithm improves " << (1 - TimeCost[11] / TimeCost[0]) * 100
		<< "\% on Edwards curve with $a=-1$ using a combination of projective coordinates and extended coordinates." << endl;

	cout << "The time cost of scalar multiplication using NAF is " << TimeCost[0] / 10000 << endl;
	cout << "The time cost of scalar multiplication using our optimal EDBC is " << TimeCost[11] / 10000 << endl;
	cout << "If one wants to know the time cost of  scalar multiplication using wNAF, greedy, tree, "
		<< " or canonic EDBC, please print the array of TimeCost in this code." << endl;




	return 0;
}