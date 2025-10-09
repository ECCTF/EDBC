#include "Supersingularity.h"
#include "preComputation.h"
#include <ctime>
 
int randomPointSupersingularityTesting(big A, int isSupersingular)
{
	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	Big n = 1;
	Big number = 1;
	for (int j = 0; j < 74; j++)
	{
		n = n*l[j];
	}
	cout << n << endl;
 
	big constantA = mirvar(1);
	big A1 = mirvar(2);
	big A2 = mirvar(1);
	big four = mirvar(4);
	nres_modsub(A, A1, A1);//A1=A-2
	nres_moddiv(A1, four, A2);
	copy(A2, constantA);//constantA=(A-2)/4
	

	projectPoint P2;
	P2.X = mirvar(123412); P2.Y = mirvar(0); P2.Z = mirvar(1);

	projectPoint P3;
	P3.X = mirvar(0); P3.Y = mirvar(0); P3.Z = mirvar(1);
	MontgomeryLadderDoubleStep(P2, P3, constantA);
	projectPoint P4;
	P4.X = mirvar(0); P4.Y = mirvar(0); P4.Z = mirvar(1);

	



	for (int i = 0; i < 33; i++)
	{
		number = n / l[i];
		MontgomeryLadder(P2, number.getbig(), P3, constantA);
		Big zero = P3.Z;
		if (zero == 0)
		{
			isSupersingular = 0;
			return -1;
		}
		big lValue = mirvar(l[i]);
		MontgomeryLadder(P3, lValue, P2, constantA);
		zero = P3.Z;
		if (zero != 0)
		{
			isSupersingular = 0;
			return -1;
		}
	}
	isSupersingular = 1;
	return 0;
}


//int ConstructMontgomeryPointCube(big A, affinePoint &P)
int ConstructMontgomeryPoint(big A, affinePoint &P)
{// 
	int nowTime = time(NULL);
	srand(nowTime);
	big u = mirvar(rand());
	big v = mirvar(0);
	big s = mirvar(0);
	big t = mirvar(0);
	big x = mirvar(0);
	big y = mirvar(5);//random choose y
	big A1 = mirvar(0);
	big A2 = mirvar(0);
	big B1 = mirvar(0);
	big B2 = mirvar(0);
	
	big a = mirvar(-1);
	big d = mirvar(1);
	big one = mirvar(1);
	big three = mirvar(3);
	big six = mirvar(6);
	big twentySeven = mirvar(27);
	big alpha = mirvar(0);
	big beta = mirvar(0);
	//y^2=x^3+ax^2+x  

	nres_modmult(u, u, s);//s=u^2 
	nres_modsub(A, s, alpha);//alpha=a-u^2
	nres_modmult(alpha, alpha, t);//t=(a-u^2)^2
	nres_modsub(three, t, B1);//beta=3-(a-u^2)^2
	nres_modmult(six, u, t);//s=6u
	nres_moddiv(B1, t, v);//v=3-(a-u^2)^2/6u	
	nres_modmult(v, v, B1);//B1=v^2
	nres_modmult(B1, twentySeven, B2);//B2=27v^2
	nres_modadd(alpha, t, s);//s=(a-u^2)^3
	nres_modadd(B2, s, A2);// A2=27v^2+(a-u^2)^3	 

	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	Big n = 1;

	for (int j = 0; j < 74; j++)
	{
		n = n*l[j];
	}
	n = n * 4;
	n = n - 1;
	n = 2 * n - 1;
	n = n / 3;//x^n=x^{1/3}

	nres_powmod(A2, n.getbig(), beta);//{27v^2+(a-u^2)^3} ^{1/3}	 
	nres_modsub(beta, alpha, t);//t={27v^2+(a-u^2)^3} ^{1/3}-a+u^2
	nres_moddiv(t, three, P.x);//beta={27v^2+(a-u^2)^3} ^{1/3}-a+u^2/3	 
	nres_modmult(u, P.x, A1);//A1=ux
	nres_modadd(A1, v, P.y);//y=ux+v

 


	return 0;
}

//int ConstructEdwardsPointCube(big A, affinePoint &P)
int ConstructEdwardsPoint(big A,  affinePoint &P)
{//ACISP 2016
	big u = mirvar(rand());
	big v = mirvar(0);
	big s = mirvar(0);
	big t = mirvar(0);
	big x = mirvar(0);
	big y = mirvar(5);//random choose y
	big A1 = mirvar(0);
	big A2 = mirvar(0);
	big B1 = mirvar(0);
	big B2 = mirvar(0);

	big a = mirvar(-1);
	big d = mirvar(1);
	big one = mirvar(1);
	big two = mirvar(2);
	big alpha = mirvar(0);
	big beta = mirvar(0);

	//y^2=x^3+ax^2+x
	 

	nres_modsub(two, A, A1);//A1=2-A
	nres_modadd(A, two, A2);//A2=A+2
	nres_moddiv(A2, A1, d);//d=A+2/2-A

	nres_modadd(a, d, B2);//B2=a+d
	
	nres_modmult(u, u, s);//s=u^2
	nres_modadd(s, s, v);//v=2u^2
	nres_modsub(B2, v, B1);//B1=a+d-2u^2

	
	
	nres_modadd(B1, B1, A1);//A1=2(a+d-2u^2) 
	nres_modmult(A1, A1, B1);//B1=4(a+d-2u^2)^2
	nres_modmult(B1, A1, B2);//B2=8(a+d-2u^2)^3
	nres_modadd(A1, A1, alpha);//alpha=4(a+d-2u^2) 

	nres_modsub(a, d, A1);//A1=a-d
	nres_modmult(A1, A1, A2);//A2=(a-d)^2

	nres_modadd(A2, A2, A1);// 
	nres_modadd(A2, A1, A2);// A2=3(a-d)^2
	nres_modsub(A2, B1, s);//s=3(a-d)^2-4(a+d-2u^2)^2

	nres_modadd(B2, B2, B1);//B1=16(a+d-2u^2)^3
	nres_modadd(s, s, A1);//A1=2s
	nres_modadd(A1, A1, A2);//A2=4s

	nres_modmult(A2, A2, A1);//A1=16s^2
	nres_modadd(A1, A1, A2);//A2=32s^2 
	nres_modadd(A2, A1, A1);//A1=48s^2
	nres_modadd(A1, A1, A2);//A2=96s^2 

	nres_modadd(A2, B1, A2);//A2=16(a+d-2u^2)^3+96s^2 

	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	Big n = 1;

	for (int j = 0; j < 74; j++)
	{
		n = n*l[j];
	}
	n = n * 4;
	n = n - 1;
	n = 2*n - 1;
	n = n / 3;//x^n=x^{1/3}
	
	nres_powmod(A2,n.getbig(),beta);//{16(a+d-2u^2)^3+96s^2} ^{1/3}

	nres_modadd(alpha, alpha, A1);//A1=8(a+d-2u^2) 
	nres_modadd(A1, A1, alpha);//alpha=16(a+d-2u^2)  
	nres_modsub(beta, alpha, t);//t=-16(a+d-2u^2) +{16(a+d-2u^2)^3+96s^2} ^{1/3}




	nres_modmult(u, t, A1);//A1=ut
	nres_modadd(A1, s, B2);//B2=ut+s

	nres_modsub(a, d, A1);//A1=a-d
	nres_modadd(A1, A1, A2);//A2=2(a-d)
	nres_modadd(A2, A2, A1);//A1=4(a-d)
	nres_modadd(A1, A1, A2);//A2=8(a-d)
	nres_modadd(A1, A2, A2);//A2=12(a-d)
	nres_modadd(A2, A2, A1);//A1=24(a-d)
	nres_modmult(A1, u, A1);//A1=24(a-d)u
	nres_modadd(t, A1, B1);//B1=t+24(a-d)u
	nres_modsub(t, A1, A2);//A2=t-24(a-d)u
	nres_modmult(B1, B2, A1);//A1=(ut+s)(t+24(a-d)u)


	nres_moddiv(one, A1, beta);//beta=1/(ut+s)(t+24(a-d))
	nres_modmult(t,B1,A1);//A1=t(t+24(a-d)u)
	nres_modmult(A1, beta, P.x);//A=t(t+24(a-d)u)

	nres_modmult(A2, B2, A1);//A1=t(ut+s)(t-24(a-d)u)
	nres_modmult(A1, beta, P.y);//A=t(t+24(a-d)u)

	//big result = mirvar(1);
	//int item[810];
	//int length;
	////n=12312451;
	//binary(n, item, length);
	//for (int i = length; i > -1; i--)
	//{
	//	nres_modmult(result, result, result);
	//	if (item[i] == 1)
	//	{
	//		nres_modmult(result, v, result);
	//	}		
	//	
	//}
	//nres_modmult(result, u, result);//result = u(v)^{p-3/4}


	return 0;
}


int productTreeSupersingularityTesting(big A, int &isSupersingular)
{
	isSupersingular = 1;
	int l[74] = {179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	Big n = 1;
	Big number = 1;
	for (int j = 0; j < 74; j++)
	{
		n = n*l[j];
	}
	cout << n << endl;

	big constantA = mirvar(1);
	big A1 = mirvar(2);
	big A2 = mirvar(1);
	big four = mirvar(4);
	nres_modsub(A, A1, A1);//A1=A-2
	nres_moddiv(A1, four, A2);
	copy(A2, constantA);//constantA=(A-2)/4

	projectPoint P2;
	P2.X = mirvar(123412); P2.Y = mirvar(0); P2.Z = mirvar(1);

	projectPoint P3;
	P3.X = mirvar(1); P3.Y = mirvar(0); P3.Z = mirvar(0);
	MontgomeryLadderDoubleStep(P2, P3, constantA);
	projectPoint P4;
	P4.X = mirvar(1); P4.Y = mirvar(0); P4.Z = mirvar(0);


	projectPoint P[65];
	for (int i = 0; i < 65; i++)
	{
		P[i].X = mirvar(1); P[i].Y = mirvar(0); P[i].Z = mirvar(0);
	}

	//Level 0
	number = 1;
	for (int i = 34; i < 74; i++)
	{
		number = number*l[i];
	}

	MontgomeryLadder(P2, number.getbig(), P[1], constantA);
	//P[1]=l[34]l[35]\ldotsl[73] P=n/(l[1]l[2]\ldotsl[33])
 

	//Level 1
	number = 1;
	for (int i = 16; i < 34; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[1], number.getbig(), P[2], constantA);
	//P[2]=n/(l[0]l[1]\ldotsl[115])P

	number = 1;
	for (int i = 0; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[1], number.getbig(), P[3], constantA);
	//P[3]=n/(l[16]l[17]\ldotsl[33])P



	//Level 2
	number = 1;
	for (int i = 8; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[2], number.getbig(), P[4], constantA);
	//P[4]=n/(l[0]l[1]\ldotsl[7])P

	number = 1;
	for (int i = 0; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[2], number.getbig(), P[5], constantA);
	//P[5]=n/(l[8]l[9]\ldotsl[15])P


 
	number = 1;
	for (int i = 24; i < 33; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[3], number.getbig(), P[6], constantA);
	//P[6]=n/(l[16]l[17]\ldotsl[23])P

	number = 1;
	for (int i = 16; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[3], number.getbig(), P[7], constantA);
	//P[7]=n/(l[24]l[25]\ldotsl[32])P

	//Level 3
	number = 1;
	for (int i = 4; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[4], number.getbig(), P[8], constantA);
	//P[8]=n/(l[0]l[1]\ldotsl[3])P

	number = 1;
	for (int i = 0; i < 4; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[4], number.getbig(), P[9], constantA);
	//P[9]=n/(l[4]l[5]\ldotsl[7])P

	number = 1;
	for (int i = 12; i < 15; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[5], number.getbig(), P[10], constantA);
	//P[10]=n/(l[8]l[9]\ldotsl[11])P
	number = 1;
	for (int i = 8; i < 12; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[5], number.getbig(), P[11], constantA);
	//P[11]=n/(l[12]l[13]\ldotsl[15])P


	number = 1;
	for (int i = 20; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[6], number.getbig(), P[12], constantA);
	//P[12]=n/(l[16]l[17]\ldotsl[19])P
	number = 1;
	for (int i = 16; i < 20; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[6], number.getbig(), P[13], constantA);
	//P[13]=n/(l[20]l[21]\ldotsl[23])P

	number = 1;
	for (int i = 28; i < 32; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[7], number.getbig(), P[14], constantA);
	//P[14]=n/(l[24]l[25]\ldotsl[27])P

	number = 1;
	for (int i = 24; i < 28; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[7], number.getbig(), P[15], constantA);
	//P[15]=n/(l[28]l[29]\ldotsl[32])P

	//////////////////////
	//Level 4
	number = 1;
	for (int i = 2; i < 4; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[8], number.getbig(), P[16], constantA);
	//P[16]=n/(l[0]l[1])P

	number = 1;
	for (int i = 0; i < 2; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[8], number.getbig(), P[17], constantA);
	//P[17]=n/(l[2]l[3])P



	number = 1;
	for (int i = 6; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[9], number.getbig(), P[18], constantA);
	//P[18]=n/(l[4]l[5])P

	number = 1;
	for (int i = 4; i < 6; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[9], number.getbig(), P[19], constantA);
	//P[19]=n/(l[6]l[7])P

	number = 1;
	for (int i = 10; i < 12; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[10], number.getbig(), P[20], constantA);
	//P[20]=n/(l[8]l[9])P
	number = 1;
	for (int i = 8; i < 10; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[10], number.getbig(), P[21], constantA);
	//P[21]=n/(l[10]l[11])P

	number = 1;
	for (int i = 14; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[11], number.getbig(), P[22], constantA);
	//P[22]=n/(l[12]l[13]\ldotsl[15])P
	number = 1;
	for (int i = 12; i < 14; i++)
	{ 
		number = number*l[i];
	}
	MontgomeryLadder(P[11], number.getbig(), P[23], constantA);
	//P[23]=n/(l[14]l[15]\ldotsl[15])P


	number = 1;
	for (int i = 18; i < 20; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[12], number.getbig(), P[24], constantA);
	//P[24]=n/(l[16]l[17])P
	number = 1;
	for (int i = 16; i < 18; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[12], number.getbig(), P[25], constantA);
	//P[25]=n/(l[18]l[19])P

	number = 1;
	for (int i = 22; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[13], number.getbig(), P[26], constantA);
	//P[26]=n/(l[20]l[21])P
	number = 1;
	for (int i = 20; i < 22; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[13], number.getbig(), P[27], constantA);
	//P[27]=n/(l[22]l[23])P

	number = 1;
	for (int i = 26; i < 28; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[14], number.getbig(), P[28], constantA);
	//P[28]=n/(l[24]l[25])P

	number = 1;
	for (int i = 24; i < 26; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[14], number.getbig(), P[29], constantA);
	//P[29]=n/(l[26]l[27])P

	number = 1;
	for (int i = 30; i < 33; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[15], number.getbig(), P[30], constantA);
	//P[30]=n/(l[28]l[29]\ldotsl[32])P

	number = 1;
	for (int i = 28; i < 30; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[15], number.getbig(), P[31], constantA);
	//P[31]=n/(l[30]l[31][32])P

	//////////////////////
	//Level 5
	number = 1;
	number = l[1]; 
	MontgomeryLadder(P[16], number.getbig(), P[32], constantA);
	//P[32]=n/(l[0])P

	number = 1;	 
	number = number*l[0];
	MontgomeryLadder(P[16], number.getbig(), P[33], constantA);
	//P[33]=n/(l[1])P

	number = 1;
	number = l[3];
	MontgomeryLadder(P[17], number.getbig(), P[34], constantA);
	//P[34]=n/(l[2])P

	number = 1;
	number = number*l[2];
	MontgomeryLadder(P[17], number.getbig(), P[35], constantA);
	//P[35]=n/(l[3])P



	number = 1;
	number = l[5];
	MontgomeryLadder(P[18], number.getbig(), P[36], constantA);
	//P[36]=n/(l[4])P

	number = 1;
	number = number*l[4];
	MontgomeryLadder(P[18], number.getbig(), P[37], constantA);
	//P[37]=n/(l[5])P

	number = 1;
	number = l[7];
	MontgomeryLadder(P[19], number.getbig(), P[38], constantA);
	//P[38]=n/(l[6])P

	number = 1;
	number = number*l[6];
	MontgomeryLadder(P[19], number.getbig(), P[39], constantA);
	//P[39]=n/(l[7])P

	number = 1;
	number = l[9];
	MontgomeryLadder(P[20], number.getbig(), P[40], constantA);
	//P[40]=n/(l[8])P

	number = 1;
	number = number*l[8];
	MontgomeryLadder(P[20], number.getbig(), P[41], constantA);
	//P[41]=n/(l[9])P

	number = 1;
	number = l[11];
	MontgomeryLadder(P[21], number.getbig(), P[42], constantA);
	//P[42]=n/(l[10])P

	number = 1;
	number = number*l[10];
	MontgomeryLadder(P[21], number.getbig(), P[43], constantA);
	//P[43]=n/(l[11])P

	number = 1;
	number = l[13];
	MontgomeryLadder(P[22], number.getbig(), P[44], constantA);
	//P[42]=n/(l[12])P

	number = 1;
	number = number*l[12];
	MontgomeryLadder(P[22], number.getbig(), P[45], constantA);
	//P[43]=n/(l[13])P

	number = 1;
	number = l[15];
	MontgomeryLadder(P[23], number.getbig(), P[46], constantA);
	//P[44]=n/(l[14])P

	number = 1;
	number = number*l[14];
	MontgomeryLadder(P[23], number.getbig(), P[47], constantA);
	//P[45]=n/(l[15])P

	number = 1;
	number = l[17];
	MontgomeryLadder(P[24], number.getbig(), P[48], constantA);
	//P[48]=n/(l[16])P

	number = 1;
	number = number*l[16];
	MontgomeryLadder(P[24], number.getbig(), P[49], constantA);
	//P[49]=n/(l[17])P

	number = 1;
	number = l[19];
	MontgomeryLadder(P[25], number.getbig(), P[50], constantA);
	//P[50]=n/(l[18])P

	number = 1;
	number = number*l[18];
	MontgomeryLadder(P[25], number.getbig(), P[51], constantA);
	//P[51]=n/(l[19])P


	number = 1;
	number = l[21];
	MontgomeryLadder(P[26], number.getbig(), P[52], constantA);
	//P[52]=n/(l[20])P

	number = 1;
	number = number*l[20];
	MontgomeryLadder(P[26], number.getbig(), P[53], constantA);
	//P[53]=n/(l[21])P
	
	number = 1;
	number = l[23];
	MontgomeryLadder(P[27], number.getbig(), P[54], constantA);
	//P[54]=n/(l[22])P

	number = 1;
	number = number*l[22];
	MontgomeryLadder(P[27], number.getbig(), P[55], constantA);
	//P[55]=n/(l[23])P

	number = 1;
	number = l[25];
	MontgomeryLadder(P[28], number.getbig(), P[56], constantA);
	//P[56]=n/(l[24])P

	number = 1;
	number = number*l[24];
	MontgomeryLadder(P[28], number.getbig(), P[57], constantA);
	//P[57]=n/(l[25])P

	number = 1;
	number = l[27];
	MontgomeryLadder(P[29], number.getbig(), P[58], constantA);
	//P[58]=n/(l[26])P

	number = 1;
	number = number*l[26];
	MontgomeryLadder(P[29], number.getbig(), P[59], constantA);
	//P[59]=n/(l[27])P


	number = 1;
	number = l[29];
	MontgomeryLadder(P[30], number.getbig(), P[60], constantA);
	//P[60]=n/(l[28])P

	number = 1;
	number = number*l[28];
	MontgomeryLadder(P[30], number.getbig(), P[61], constantA);
	//P[61]=n/(l[29])P

	number = 1;
	number = l[31]* l[32];
	MontgomeryLadder(P[31], number.getbig(), P[62], constantA);
	//P[62]=n/(l[30])P

	number = 1;
	number = l[30] * l[32];
	MontgomeryLadder(P[31], number.getbig(), P[63], constantA);
	//P[63]=n/(l[31])P

	number = 1;
	number = l[30] * l[31];
	MontgomeryLadder(P[31], number.getbig(), P[64], constantA);
	//P[64]=n/(l[32])P

	
	for (int i = 32; i < 64; i++)
	{
		if(P[i].Z ==0)
		{
			isSupersingular = 0;
			return -1;
		}
	}
	number=l[32];
	MontgomeryLadder(P[64], number.getbig(), P[0], constantA);
	
	if (P[0].Z != 0)
	{
		isSupersingular = 0;
		return -1;
	}
	return 0;
}

int productTree(big A, int &IsSupersingularity)
{
	productTreeSupersingularityTesting(A, IsSupersingularity);
	return IsSupersingularity;
}

int Scalarmultiplication1to32(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsDouble(P, PProjective1, d); //PProjective1=2P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^2P
	EdwardsDouble(PProjective2, PProjective1, d); //PProjective1=2^3P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^4P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^5P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^5P
	EdwardsAddition(PExtP2,P,P1,d);//P1=-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^6P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^7P
	EdwardsDouble(PProjective2, PExtP2, d); //PExtP2=2^8P
	EdwardsAddition(PExtP2, P1, PMiddle1, d);//PMiddle1=2^8-2^5+1

	EdwardsDouble(PExtP2, PExtP1, d); //PExtP1=2^9P
	EdwardsAddition(PExtP1, P1, P2, d);//P2=2^9-2^5+1

	EdwardsDouble(PExtP1, PExtP2, d); //PExtP2=2^{10}P
	EdwardsAddition(PExtP2, PMiddle1, P1, d);//P1=2^{10}+2^8-2^5+1

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^{11}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{12}P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{13}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{13}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{14}P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{15}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{15}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{16}P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{17}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{17}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^{18}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{18}+2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{19}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{20}P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{21}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{21}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{22}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{23}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{23}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{23}+2^{18}+2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{24}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{25}P
	EdwardsDouble(PProjective2, PProjective2, d); //PProjective2=2^{26}P
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{27}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{27}-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{27}+2^{23}+2^{18}+2^{15}+2^{13}+2^9-2^5+1		

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{27}3^{1}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{28}3^{1}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{28}3^{1}+2^{27}+2^{23}+2^{18}+2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{29}3^{1}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{29}3^{1}+2^{27}-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{30}3^{1}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{31}3^{1}P	
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{32}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{32}3^{1}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=-2^{32}3^{1}+2^{29}3^{1}+2^{27}-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=-2^{32}3^{1}+2^{28}3^{1}+2^{27}+2^{23}+2^{18}+2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{33}3^{1}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{34}3^{1}P	
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{35}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{35}3^{1}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{35}3^{1}-2^{32}3^{1}+2^{29}3^{1}+2^{27}-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{35}3^{1}-2^{32}3^{1}+2^{28}3^{1}+2^{27}+2^{23}+2^{18}+2^{15}+2^{13}+2^9-2^5+1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^{36}3^{1}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{36}3^{1}+2^{35}3^{1}-2^{32}3^{1}+2^{29}3^{1}+2^{27}-2^{23}-2^{21}-2^{17}-2^{15}-2^{13}+2^{10}+2^8-2^5+1
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{37}3^{1}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{38}3^{1}P	
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{39}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{39}3^{1}P	 
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{39}3^{1}...
	
	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^{40}3^{1}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{40}3^{1}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{40}3^{2}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{40}3^{2}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{41}3^{2}P
	EdwardsDouble(PProjective1, PProjective2, d); //PProjective2=2^{42}3^{2}P	
	EdwardsDouble(PProjective2, PExtP1, d); //PExtP1=2^{43}3^{2}P	 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{43}3^{2}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{44}3^{2}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{45}3^{2}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{45}3^{2}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{46}3^{2}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{47}3^{2}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{47}3^{2}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{47}3^{3}P
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^{47}3^{3}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{48}3^{3}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{48}3^{3}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{48}3^{3}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{49}3^{3}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{50}3^{3}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{50}3^{3}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{51}3^{3}P
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{51}3^{3}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{51}3^{4}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{52}3^{4}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{53}3^{4}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{53}3^{4}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{53}3^{5}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{54}3^{5}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{54}3^{5}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{55}3^{5}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{56}3^{5}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{57}3^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{58}3^{5}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{58}3^{5}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{58}3^{6}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{58}3^{6}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{58}3^{6}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{58}3^{7}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{58}3^{7}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{58}3^{7}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{59}3^{7}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{60}3^{7}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{60}3^{7}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{60}3^{7}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{61}3^{7}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{61}3^{7}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=-2^{61}3^{7}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{61}3^{7}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{61}3^{8}P
	EdwardsDouble(PProjective1, PProjective1, d); //PExtP1=2^{62}3^{8}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{63}3^{8}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{63}3^{8}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{63}3^{8}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{64}3^{8}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{64}3^{8}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=-2^{64}3^{8}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{65}3^{8}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{65}3^{8}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{66}3^{8}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{67}3^{8}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{67}3^{8}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{67}3^{8}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{68}3^{8}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{69}3^{8}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{70}3^{8}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{71}3^{8}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{71}3^{8}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{71}3^{8}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{72}3^{8}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{72}3^{8}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{72}3^{8}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{73}3^{8}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{74}3^{8}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{74}3^{8}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{74}3^{8}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{75}3^{8}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{75}3^{8}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{75}3^{8}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{75}3^{9}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{75}3^{9}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=-2^{75}3^{9}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{75}3^{10}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{75}3^{10}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=-2^{75}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{76}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{77}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{78}3^{10}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{78}3^{10}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{78}3^{10}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{78}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{79}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{80}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{81}3^{10}P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{81}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{82}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{83}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{84}3^{10}P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{84}3^{10}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{85}3^{10}P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{85}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{86}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{87}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{88}3^{10}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{88}3^{10}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{88}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{89}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{90}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{91}3^{10}P		
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{91}3^{10}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{91}3^{10}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{91}3^{10}...


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{92}3^{10}P		
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{93}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{94}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{95}3^{10}P		
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{95}3^{10}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{95}3^{10}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{95}3^{11}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{95}3^{11}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^{95}3^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{96}3^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{97}3^{11}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{97}3^{11}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{97}3^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{98}3^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{99}3^{11}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{99}3^{11}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{99}3^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{100}3^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{101}3^{11}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{101}3^{11}...

	EdwardsDouble(PExtP1,PExtP1, d); //PExtP1=2^{102}3^{11}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{102}3^{11}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{103}3^{11}P
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{103}3^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{104}3^{11}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{105}3^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{106}3^{11}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{106}3^{11}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{106}3^{11}...
	EdwardsAddition(P2, PExtP2, P2, d);//P1=-2^{106}3^{11}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{107}3^{11}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{107}3^{11}P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{107}3^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{108}3^{11}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{109}3^{11}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{110}3^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{111}3^{11}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{111}3^{11}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{111}3^{11}...
	return 0;
}

int Scalarmultiplication1to16(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsDouble(P, PProjective1, d); //PProjective1=2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^3P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^3P
	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P
	EdwardsAddition(PExtP2, PMiddle1, P1, d);//P1=-2^3-1
		
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{4}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^4P
	EdwardsAddition(P1, PExtP2, P2, d);//P2=-2^{4}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{6}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{6}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{6}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{6}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{8}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{9}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{11}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{11}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{12}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{13}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{14}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{14}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{15}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{15}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{15}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{15}3^{1}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{16}3^{1}P		
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{16}3^{1}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{16}3^{1}...
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{17}3^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{18}3^{1}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{19}3^{1}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{19}3^{1}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{20}3^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{21}3^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{22}3^{1}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{23}3^{1}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{23}3^{1}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{24}3^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{25}3^{1}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{26}3^{1}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{26}3^{1}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{26}3^{1}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{26}3^{1}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{26}3^{2}P
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{26}3^{2}...

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{26}3^{3}P
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{26}3^{3}...
		
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{27}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{27}3^{3}P		
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{27}3^{3}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{28}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{28}3^{3}P		
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{28}3^{3}...
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{29}3^{3}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{30}3^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{30}3^{3}P		
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{30}3^{3}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{31}3^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{31}3^{3}P		
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{31}3^{3}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{32}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{33}3^{3}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{34}3^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{34}3^{3}P		
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{34}3^{3}...
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{35}3^{3}P	 
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{36}3^{3}P	 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{36}3^{3}...
	
	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{36}3^{4}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{37}3^{4}P	 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{37}3^{4}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{38}3^{4}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{39}3^{4}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{40}3^{4}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{40}3^{4}P		
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{40}3^{4}...
	
	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{40}3^{5}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{40}3^{5}P		
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{40}3^{5}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{41}3^{5}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{42}3^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{43}3^{5}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{43}3^{5}P		
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{43}3^{5}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{44}3^{5}P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{45}3^{5}P		
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^{45}3^{5}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{45}3^{6}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{46}3^{6}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{46}3^{6}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{46}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{47}3^{7}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{48}3^{7}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{48}3^{7}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{49}3^{7}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{49}3^{7}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{50}3^{7}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{50}3^{7}...


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{51}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{52}3^{7}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{53}3^{7}P			
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{53}3^{7}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{53}3^{7}...

	return 0;
}

int Scalarmultiplication17to32(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsDouble(P, PProjective1, d); //PProjective1=2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^2P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^3P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^4P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^5P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^5P	
	EdwardsAddition(PExtP2, P, P2, d);//P2=-2^5-1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{6}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{7}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^7P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^{7}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{7}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{8}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^8P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{8}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{9}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{10}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{11}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{12}P
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^{12}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{13}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{14}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{14}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{15}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{15}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{15}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{15}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{16}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{17}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{17}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{17}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{18}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{18}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{18}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{19}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{20}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{20}P
	EdwardsAddition(P2, PExtP2, P2, d);//P1=-2^{20}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{21}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{21}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{21}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{22}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{23}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{23}P
	EdwardsAddition(P2, PExtP2, P2, d);//P1=-2^{23}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{24}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{24}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{24}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{25}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{25}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{25}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{26}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{26}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{26}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{27}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{28}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{28}P
	EdwardsAddition(P2, PExtP2, P2, d);//P1=-2^{28}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{29}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{29}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{29}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{29}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{29}3^{1}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{29}3^{2}P
	EdwardsTriple(PProjective1, PExtP1, d); //PProjective1=2^{29}3^{3}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{29}3^{3}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{29}3^{3}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{30}3^{3}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{31}3^{3}P		
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{31}3^{3}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{31}3^{3}...

	EdwardsTriple(PExtP1, PExtP1, d); //PProjective1=2^{31}3^{4}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{31}3^{4}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{32}3^{4}P		 
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{33}3^{4}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{33}3^{4}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{34}3^{4}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{34}3^{4}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{34}3^{4}...


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{35}3^{4}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{36}3^{4}P			
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{36}3^{4}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{36}3^{4}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{36}3^{4}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{37}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{38}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{39}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{40}3^{4}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{41}3^{4}P		 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{41}3^{4}...
	 

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{42}3^{4}P		 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{42}3^{4}...	 

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^{42}3^{5}P		 
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^{42}3^{5}... 

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{43}3^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{44}3^{5}P		 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{44}3^{5}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{45}3^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{46}3^{5}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{46}3^{5}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{46}3^{5}...

	EdwardsDouble(PExtP1,PExtP1, d); //PExtP1=2^{47}3^{5}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{47}3^{5}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{47}3^{5}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{48}3^{5}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{48}3^{5}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{48}3^{6}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{49}3^{6}P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{50}3^{6}P		 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{50}3^{6}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{51}3^{6}P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{52}3^{6}P		 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{52}3^{6}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{53}3^{6}P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{54}3^{6}P		 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{54}3^{6}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{55}3^{6}P		 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{55}3^{6}...
	return 0;
}

int Scalarmultiplication1to8(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsDouble(P, P2, d); //P2=2^1P	
	EdwardsDouble(P2, PExtP1, d); //PExtP1=2^2P
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^2+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^2+2^1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{3}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{5}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{6}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{8}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{9}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{9}... 

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{9}3^{1}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{9}3^{2}P
	EdwardsTriple(PProjective1, PExtP1, d); //PProjective1=2^{9}3^{3}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{9}3^{3}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{9}3^{4}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{9}3^{5}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{9}3^{6}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{9}3^{7}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{10}3^{7}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{10}3^{7}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{10}3^{7}...
 
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{11}3^{7}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{12}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{13}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{14}3^{7}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{15}3^{7}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{15}3^{7}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{15}3^{7}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{16}3^{7}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{17}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{18}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{19}3^{7}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{20}3^{7}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{21}3^{7}P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{21}3^{7}...

	EdwardsAddition(P1, P2, P2, d);//P1=2^{21}3^{7}...
	return 0;
}

int Scalarmultiplication9to16(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);


	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P
	EdwardsDouble(P, P2, d); //P2=2P

	EdwardsTriple(P2, PExtP1, d); //PExtP1=2^{1}3^{1}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{1}3^{1}P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P2=-2^{1}3^{1}-1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^{2}3^{1}P
	EdwardsAddition(P2, PExtP1, P1, d);//P2=2^{2}3^{1}+2^1
	
	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{1}3^{2}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{2}3^{2}P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{2}3^{2}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{2}3^{2}...

	
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{3}3^{2}P			
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^{3}3^{2}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{3}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{4}3^{3}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{5}3^{3}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{5}3^{3}...
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{6}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{7}3^{3}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{8}3^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{8}3^{3}P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{8}3^{3}...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{8}3^{4}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{8}3^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{9}3^{5}P			
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{9}3^{5}...


	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{8}3^{6}P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{8}3^{7}P	
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{8}3^{8}P	
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{8}3^{9}P	
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{8}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{9}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{10}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{11}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{12}3^{10}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{13}3^{10}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PProjective1=2^{14}3^{10}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{14}3^{10}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{15}3^{10}P
	EdwardsDouble(PProjective1, PExtP1, d); //PProjective1=2^{16}3^{10}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{16}3^{10}...

	EdwardsAddition(P1, P2, P2, d);//P2= 
	return 0;
}
 
int Scalarmultiplication17to24(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);


	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{2}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{3}P
	EdwardsAddition(PMiddle1, PExtP1, P1, d);//P1=2^{3}-1
	EdwardsAddition(PMiddle1, PExtP2, P2, d);//P2=-2^{3}-1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{4}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{5}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{6}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{6}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{7}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{8}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{9}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{9}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{9}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{10}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{11}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{12}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{12}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{13}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{14}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{14}... 

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{15}P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{15}... 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{15}... 

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{16}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{17}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{18}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{18}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{18}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{19}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{20}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{20}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{20}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{21}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{21}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{21}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{22}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{22}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{22}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{23}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{24}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{24}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{24}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{25}P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^{25}...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{26}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{26}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{26}...


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{27}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{28}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{29}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{29}P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{29}...
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{29}...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{30}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{31}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{32}P		
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{32}...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{32}...
	return 0;
}

int Scalarmultiplication25to32(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^{1}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{1}3^{1}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{2}3^{1}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{3}3^{1}P	 
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^{3}3^{1}+1
 
	EdwardsTriple(P, PProjective1, d); //PProjective1=2^{3}3^{2}P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{4}3^{2}P	
	EdwardsNegative(PExtP1,P1, d);//P1=-2^{4}3^{2}P
	 

	EdwardsTriple(P, PProjective1, d); //PProjective1=2^{4}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{5}3^{3}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{6}3^{3}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{7}3^{3}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{8}3^{3}P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{8}3^{3}P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^{8}3^{3}...
	 
	EdwardsDouble(P, PProjective1, d); //PProjective1=2^{9}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{10}3^{3}P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{11}3^{3}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{11}3^{3}...

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^{12}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{13}3^{3}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{14}3^{3}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{15}3^{3}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{16}3^{3}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{16}3^{3}...

	EdwardsTriple(P, PProjective1, d); //PProjective1=2^{16}3^{4}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{17}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{18}3^{4}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{19}3^{4}P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{19}3^{4}...

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^{20}3^{4}P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{21}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{22}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{23}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{24}3^{4}P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{25}3^{4}P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{26}3^{4}P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{26}3^{4}...

	EdwardsAddition(P2, P1, P1, d);//P1=2^{26}3^{4}...
	return 0;
}



int Scalarmultiplication1to4(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^13^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^23^1+1
	 
	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^23^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^2P
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^33^2+1
	 
	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^33^3P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^43^3P
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^43^3...
	 
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^53^3P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^3...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^63^3P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^3...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^63^4P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^73^4P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^83^4P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^83^4...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^93^4P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^93^4...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^103^4P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{10}3^4...
	return 0;
}

int Scalarmultiplication5to8(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^03^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^2P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^2P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^13^2+1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^13^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^3P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^3P
	EdwardsAddition(P, PExtP1, P2, d);//P2=-2^23^3+1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^33^3P 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^33^3...


	EdwardsTriple(PExtP1, PExtP1, d); //PProjective1=2^33^4P 	
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^33^4...


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^43^4P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^53^4P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^4...
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^4...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^63^4P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^63^4...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^4...

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^73^4P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^83^4P	
	EdwardsDouble(PProjective1, PExtP1, d); //PProjective1=2^93^4P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^93^4...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{10}3^4P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{10}3^4...

	return 0;
}

int Scalarmultiplication9to12(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P


	EdwardsDouble(P, PProjective1, d); //PProjective1=2^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^2P
	EdwardsNegative(PExtP1, P1, d);//P1=-2^2P
	 

	EdwardsTriple(P1, PProjective1, d); //PProjective1=2^23^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^33^1P
	EdwardsAddition(PMiddle1, PExtP2, P2, d);//P2=-2^33^1-1

	EdwardsDouble(PExtP2, PProjective1, d); //PProjective1=2^43^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^53^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^63^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^73^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^83^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^93^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{10}3^1P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{10}3^1-2^2

	EdwardsTriple(P1, PProjective1, d); //PProjective1=2^{10}3^2P	
	EdwardsTriple(PProjective1, PExtP1, d); //PExtP1=2^{10}3^3P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{10}3^3P	
	EdwardsAddition(P1, PExtP2, P1, d);//P2=-2^{10}3^3+2^{10}3^1-2^2


	EdwardsTriple(P1, PExtP1, d); //PExtP1=2^{10}3^4P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{10}3^4+...

	EdwardsAddition(P2, P1, P1, d);//P1=2^{10}3^4+...
	return 0;
}


int Scalarmultiplication13to16(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^1P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P1=-2^13^1-1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^13^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^2P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^2P
	EdwardsAddition(PMiddle1, PExtP2, P2, d);//P2=-2^23^2-1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^23^3P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^3P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^33^3P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^33^3...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^43^3P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^53^3P	 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^3...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^63^3P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^73^3P	 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^73^3...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^83^3P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{9}3^3P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{10}3^3P	 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{10}3^3...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{11}3^3P	 
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{11}3^3...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{11}3^3...
	 
	return 0;
}

int Scalarmultiplication17to20(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsDouble(P, P1, d); //P1=2^1P
	
	EdwardsDouble(P1, PProjective1, d); //PProjective1=2^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^3P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^3P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^3+2

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^4P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^5P	
	EdwardsAddition(PMiddle1, PExtP2, P2, d);//P2=2^5-1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^53^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{6}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{7}3^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^83^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^83^1P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^83^1...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^93^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{10}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{7}3^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{11}3^1P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{11}3^1...
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^{11}3^1...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^{11}3^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{12}3^2P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{12}3^2...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{12}3^2...
 
	EdwardsAddition(P1, P2, P1, d);//P1=2^{12}3^2...

	return 0;
}

int Scalarmultiplication21to24(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);
	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P
	
	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, P1, d); //P1=2^13^1P
 
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^23^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^1P	
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^33^1+1
	
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^43^1P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^43^1+...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^53^1P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^53^1+...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^63^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^73^1P	
	EdwardsAddition(P1, PExtP1, P1, d);//P2=2^73^1+...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^83^1P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^83^1+...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^93^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{10}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{11}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{12}3^1P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{12}3^1P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{12}3^1...
	
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{13}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{14}3^1P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^{14}3^1P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=-2^{14}3^1...
	
	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{15}3^1P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{15}3^1...

	EdwardsAddition(P1, P2, P1, d);//P1=...

	return 0;
}


int Scalarmultiplication25to28(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);
	
	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{1}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^23^1+1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^23^2P	
	EdwardsTriple(PProjective1, PExtP1, d); //PExtP1=2^23^3P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^3P
	EdwardsAddition(P, PExtP1, P2, d);//P2=-2^23^3...

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^23^4P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=2^{1}3^5P	
	EdwardsTriple(PProjective1, PExtP1, d); //PExtP1=2^23^6P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^23^6...
	
	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^23^7P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^23^7...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^33^7P	
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^33^7...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^43^7P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^43^7...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^43^7...
		 
	return 0;
}


int Scalarmultiplication29to32(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{1}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{2}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^33^1P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^33^1+1

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^43^1P	
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^43^1+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^53^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^63^1P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^63^1...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^1...

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^73^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{8}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^93^1P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^93^1...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{10}3^1P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{10}3^1...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{10}3^1... 

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^{11}3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{12}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{13}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^{14}3^1P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^{14}3^1...

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^{15}3^1P		
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^{15}3^1...
	
	return 0;
}


int Scalarmultiplication1to2(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{1}3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^1P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P1=-2^23^1-1 

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^33^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^43^1P	
	EdwardsAddition(PMiddle1, PExtP1, P2, d);//P2=2^43^1-1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^53^1P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^63^1P
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^63^1...
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^1...
	
	return 0;
}

int Scalarmultiplication3to4(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{1}3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{2}3^1P	
	EdwardsDouble(PProjective1, P2, d); //P2=2^33^1P	
	EdwardsAddition(P, P2, P1, d);//P1=2^33^1+1 

	EdwardsTriple(P2, PProjective1, d); //PProjective1=2^33^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^{4}3^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //P2=2^53^2P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^2+2^33^1+1 
	
	EdwardsAddition(P2, P1, P2, d);//P2=2^33^1+P1

	return 0;
}

int Scalarmultiplication5to6(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=3^2P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^2P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^2P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P1=-2^13^2-1

	EdwardsTriple(PExtP1, P2, d); //P2=2^13^3P	

	EdwardsDouble(P2, PProjective1, d); //PProjective1=2^23^3P	
	EdwardsDouble(PProjective1, PExtP1, d); //P2=2^33^3P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^33^3-2^13^2-1

	EdwardsAddition(P2, P1, P2, d);//P2=2^13^3+P1

	return 0;
}

int Scalarmultiplication7to8(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^13^1P		
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^1P
	EdwardsAddition(PMiddle1, PExtP2, P2, d);//P2=-2^23^1-1

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^33^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^33^1P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P1=-2^33^1-1

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^33^2P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^33^2-2^33^1-1

	EdwardsDouble(P2, PProjective1, d); //PProjective1=2^43^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //P2=2^53^2P	
	EdwardsAddition(P1, PExtP2, P1, d);//P1=2^53^2-2^33^1-1
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^53^2+2^33^2-2^33^1-1

	return 0;
}


int Scalarmultiplication9to10(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P
	EdwardsDouble(P, PExtP1, d); //PExtP1=2^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^1P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^1+1


	 

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^3P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^4P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^5P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^5-2^1+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^6P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^7P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^8P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^8+2^5-2^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P1=2^8+1

	return 0;
}


int Scalarmultiplication11to12(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);
		
	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P


	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=3^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^1P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^1P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^13^1+1
	EdwardsNegative(P1, P2, d);//P1=2^13^1-1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^13^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^23^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^33^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^43^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^53^2P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^2-2^13^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^53^2+2^13^1-1

	return 0;
}

int Scalarmultiplication13to14(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P


	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^3P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^3P
	EdwardsAddition(P, PExtP2, P2, d);//P2=-2^3+1

	EdwardsTriple(PExtP1, PExtP1, d); //PProjective1=2^33^1P
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^33^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^33^1-2^3+1 

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^33^2P
	EdwardsTriple(PProjective1, PExtP1, d); //PExtP1=2^33^3P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^33^3+2^33^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^33^3+2^33^1-2^3+1 

	return 0;
}

int Scalarmultiplication15to16(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P


	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=3^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^13^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^1P	
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^23^1+1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^33^1P	
	EdwardsAddition(PMiddle1, PExtP1, P2, d);//P2=2^33^1-1 

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^33^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^43^2P
	EdwardsTriple(PProjective1, PExtP1, d); //PExtP1=2^53^2P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^53^2-2^23^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^53^2+2^33^1-1  

	return 0;
}


int Scalarmultiplication17to18(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P
	EdwardsDouble(P, PProjective1, d); //PProjective1=2^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^2P	 
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^2P
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^2+1
	EdwardsAddition(P, PExtP2, P2, d);//P2=-2^2+1

	EdwardsDouble(PExtP1,PExtP1, d); //PExtP1=2^3P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^3+2^2+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^4P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^5P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^6P	
	EdwardsAddition(P2, PExtP2, P2, d);//P2=2^6-2^2+1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^6P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^7P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^8P	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^8+2^5-2^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^8+1

	return 0;
}

int Scalarmultiplication19to20(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsDouble(P, PProjective1, d); //PProjective1=2^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^2P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^3P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^4P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^5P	 	
	EdwardsAddition(PMiddle1, PExtP1, P2, d);//P2=2^5-1


	 

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^53^1P 
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^63^1P		
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^63^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^1+2^5-1

	return 0;
}

int Scalarmultiplication21to22(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^1P	
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^13^1+1


	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^1P	
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^23^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^23^1+2^13^1+1




	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^23^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^33^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^43^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^53^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^63^1P	 	
	EdwardsAddition(P1, PExtP1, P1, d);//P1=2^63^1-2^23^1+1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^1+2^23^1+2^13^1+1

	return 0;
}

int Scalarmultiplication23to24(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^1P	
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^1P
	EdwardsAddition(PMiddle1, PExtP2, P1, d);//P1=-2^13^1-1

	EdwardsTriple(PExtP1, PProjective1, d); //PProjective1=2^13^2P 
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^23^2P
	EdwardsDouble(PProjective1, P2, d); //P2=2^33^2		

	EdwardsTriple(P2, PProjective1, d); //PProjective1=2^43^2P 	
	EdwardsDouble(PProjective1, P2, d); //PExtP1=2^53^2P	
	EdwardsAddition(P, PExtP1, P1, d);//P1=2^53^2-2^13^1-1	
	EdwardsAddition(P1, P2, P2, d);//P2=P1+2^33^2		

	return 0;
}

int Scalarmultiplication25to26(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^13^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^23^1+ 1

	EdwardsDouble(PExtP1, PExtP1, d); //PProjective1=2^33^1P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^33^1+2^23^1+ 1
	
	EdwardsDouble(P, PProjective1, d); //PProjective1=2^43^1P
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^53^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^63^1P
	EdwardsAddition(PMiddle1, PExtP1, P1, d);//P1=2^63^1-1

	 
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^63^1+2^33^1+2^23^1+ 1

	return 0;
}

int Scalarmultiplication27to28(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^13^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^23^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^23^1P
	EdwardsAddition(P, PExtP2, P2, d);//P2=-2^23^1+ 1

 

	EdwardsTriple(P, PProjective1, d); //PProjective1=2^23^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^2P
	EdwardsAddition(PMiddle1, PExtP1, P1, d);//P1=2^33^2-1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^33^2-2^23^1+ 1


	EdwardsDouble(P, PProjective1, d); //PProjective1=2^43^2P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^53^2P
	EdwardsAddition(PMiddle1, PExtP1, P1, d);//P1=2^53^2+2^33^2-1
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^53^2+2^33^2-2^23^1+ 1

	return 0;
}


int Scalarmultiplication29to30(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 
	EdwardsTriple(PProjective1, PProjective1, d); //PProjective1=3^2P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^2P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^2P
	EdwardsAddition(P, PExtP2, P1, d);//P1=-2^13^2+ 1

	EdwardsTriple(PExtP1, PExtP1, d); //PExtP1=2^13^3P	
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^13^3+ 1
 
	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^23^3P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^3P
	EdwardsAddition(P1, PExtP2, P1, d);//P1=2^33^3-2^13^2+ 1
	EdwardsAddition(P, PExtP1, P2, d);//P2=2^33^3+2^13^3+ 1

	return 0;
}

int Scalarmultiplication31to32(extProjectPoint P, extProjectPoint P1, extProjectPoint P2, big d)
{//d is the curve parameter
	extProjectPoint PMiddle1;
	PMiddle1.X = mirvar(0); PMiddle1.Y = mirvar(0); PMiddle1.Z = mirvar(1), PMiddle1.T = mirvar(1);
	extProjectPoint PMiddle2;
	PMiddle2.X = mirvar(0); PMiddle2.Y = mirvar(0); PMiddle2.Z = mirvar(1), PMiddle2.T = mirvar(1);

	extProjectPoint PExtP1;
	PExtP1.X = mirvar(0); PExtP1.Y = mirvar(0); PExtP1.Z = mirvar(1), PExtP1.T = mirvar(1);
	extProjectPoint PExtP2;
	PExtP2.X = mirvar(0); PExtP2.Y = mirvar(0); PExtP2.Z = mirvar(1), PExtP2.T = mirvar(1);

	projectPoint PProjective1;
	PProjective1.X = mirvar(0); PProjective1.Y = mirvar(0); PProjective1.Z = mirvar(1);
	projectPoint PProjective2;
	PProjective2.X = mirvar(0); PProjective2.Y = mirvar(0); PProjective2.Z = mirvar(1);

	//EdwardsNegative(P, PMiddle1, d);//PMiddle1=-P

	EdwardsTriple(P, PProjective1, d); //PProjective1=3^1P 	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^13^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^13^1P
	EdwardsAddition(P, PExtP2, P2, d);//P2=-2^13^1+ 1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^23^1P	
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^33^1P
	EdwardsNegative(PExtP1, P1, d);//P1=-2^33^1

	EdwardsDouble(PExtP1, PExtP1, d); //PExtP1=2^43^1P
	EdwardsNegative(PExtP1, PExtP2, d);//PExtP2=-2^43^1P
	EdwardsAddition(P2, PExtP2, P2, d);//P2=-2^43^1-2^13^1+ 1

	EdwardsDouble(PExtP1, PProjective1, d); //PProjective1=2^53^1P	
	EdwardsDouble(PProjective1, PProjective1, d); //PProjective1=2^63^1P
	EdwardsDouble(PProjective1, PExtP1, d); //PExtP1=2^73^1P
	EdwardsAddition(P2, PExtP1, P2, d);//P2=2^73^1-2^43^1-2^13^1+ 1

	EdwardsAddition(P1, P2, P1, d);//P1=P2-2^33^1
	return 0;
}

int newProductTree(big A, doubleBaseType DBCTerm[], int DBCLength, int &IsSupersingularity)
{
	affinePoint P;
	P.x = mirvar(1); P.y = mirvar(1);
	ConstructEdwardsPoint(A, P);
	//construct a point on Edwards curve

	IsSupersingularity = 1;
	big A1 = mirvar(0); 
	big A2 = mirvar(0);
	big d = mirvar(1); 
	big two = mirvar(2);
	big one = mirvar(1);
 

	//y^2=x^3+ax^2+x

	nres_modsub(two, A, A1);//A1=2-A
	nres_modadd(A, two, A2);//A2=A+2
	nres_moddiv(A2, A1, d);//d=A+2/2-A


	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	Big n = 1;

	for (int j = 32; j < 74; j++)
	{
		n = n*l[j];
	}
	n = n * 4;
 
	//DBC Item[MatrixB][MatrixT][2 * Cmax];
	DBCChainItem finalChain[MatrixB];
	long long CPUcycles;


	affinePoint P2;
	
	extProjectPoint PD[65];//P[i]=iP for i\in I_w
	projectPoint Pre[65];
	projectPoint P1[2];
	
	 
	big X1[2];
	big Y1[2];
	big Z1[2];
	for (int i = 0; i < 2; i++)
	{
		X1[i] = mirvar(0);
		P1[i].X = X1[i];
		Y1[i] = mirvar(0);
		P1[i].Y = Y1[i];
		Z1[i] = mirvar(1);
		P1[i].Z = Z1[i];
		 
	}

	P2.x = mirvar(2);
	P2.y = mirvar(5);
	copy(P.x, P2.x);
	copy(P.y, P2.y);
	nres(P2.x, P2.x);
	nres(P2.y, P2.y);
	affinePoint P4;
	P4.x = mirvar(2);
	P4.y = mirvar(5);
	
	big X[65];
	big Y[65];
	big Z[65];
	big T[65];
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

	preComputationEx(P2, PD, CPUcycles, 7, d);
	windowScalarMultiplicationEx(P2, P4, DBCTerm, DBCLength, PD, 7);//(p+1)/l_1...l_32
	//Scalar512Multiplication1to32(P2, P4);	 
	copy(P4.x,PD[0].X);
	copy(P4.y, PD[0].Y);
	nres_modmult(PD[0].X, PD[0].X, PD[0].T);
	copy(one, PD[0].Z);
	Scalarmultiplication1to32(PD[0],PD[2],PD[3],d);
	Scalarmultiplication1to16(PD[2], PD[4], PD[5], d);
	Scalarmultiplication17to32(PD[3], PD[6], PD[7], d);

	Scalarmultiplication1to8(PD[4], PD[8], PD[9], d);
	Scalarmultiplication9to16(PD[5], PD[10], PD[11], d);
	Scalarmultiplication17to24(PD[6], PD[12], PD[13], d);
	Scalarmultiplication25to32(PD[7], PD[14], PD[15], d);

	Scalarmultiplication1to4(PD[8], PD[16], PD[17], d);
	Scalarmultiplication5to8(PD[9], PD[18], PD[19], d);
	Scalarmultiplication9to12(PD[10], PD[20], PD[21], d);
	Scalarmultiplication13to16(PD[11], PD[22], PD[23], d);
	Scalarmultiplication17to20(PD[12], PD[24], PD[25], d);
	Scalarmultiplication21to24(PD[13], PD[26], PD[27], d);
	Scalarmultiplication25to28(PD[14], PD[28], PD[29], d);
	Scalarmultiplication29to32(PD[15], PD[30], PD[31], d);

	Scalarmultiplication1to2(PD[16], PD[32], PD[33], d);
	Scalarmultiplication3to4(PD[17], PD[34], PD[35], d);
	Scalarmultiplication5to6(PD[18], PD[36], PD[37], d);
	Scalarmultiplication7to8(PD[19], PD[38], PD[39], d);
	Scalarmultiplication9to10(PD[20], PD[40], PD[41], d);
	Scalarmultiplication11to12(PD[21], PD[42], PD[43], d);
	Scalarmultiplication13to14(PD[22], PD[44], PD[45], d);
	Scalarmultiplication15to16(PD[23], PD[46], PD[47], d);
	Scalarmultiplication17to18(PD[24], PD[48], PD[49], d);
	Scalarmultiplication19to20(PD[25], PD[50], PD[51], d);
	Scalarmultiplication21to22(PD[26], PD[52], PD[53], d);
	Scalarmultiplication23to24(PD[27], PD[54], PD[55], d);
	Scalarmultiplication25to26(PD[28], PD[56], PD[57], d);
	Scalarmultiplication27to28(PD[29], PD[58], PD[59], d);
	Scalarmultiplication29to30(PD[30], PD[60], PD[61], d);
	Scalarmultiplication31to32(PD[31], PD[62], PD[63], d);

	Big number = 1;

	for (int i = 32; i < 64; i++)
	{
		number = PD[i].Z;

		if (number == 0)
		{
			IsSupersingularity = 0;
			return -1;
		}
	}
	number = 179;
	copy(PD[32].X, P1[1].X);
	copy(PD[32].Y, P1[1].Y);
	copy(PD[32].Z, P1[1].Z);
	MontgomeryLadder(P1[1], number.getbig(), P1[0], d);

	number = PD[0].Z;
	if (number != 0)
	{
		IsSupersingularity = 0;
		return -1;
	}

	return 1;
}

int preComputationEx1(affinePoint P1, extProjectPoint *P, long long &CPUcycles, int Cmax, big d)
{
	big one = mirvar(1);
	copy(P1.x, P[1].X);
	copy(P1.y, P[1].Y);
	copy(one, P[1].Z);
	nres_modmult(P[1].X, P[1].Y, P[1].T);                       

	
	
	return 0;
}


int cDBL(MontgomeryPoint P, MontgomeryPoint P2, big A)
{
	big A24 = mirvar(1);
	big two = mirvar(2);
	big four = mirvar(4);
	big t0 = mirvar(2);
	big t1 = mirvar(4);
	big t2 = mirvar(4);
	big middle0 = mirvar(2);
	big middle1 = mirvar(4);
	big middle2 = mirvar(4);
	
	nres_modadd(A, two, A24);
	nres_moddiv(A24, four, A24);//A24=(A + 2) / 4;

	nres_modadd(P.X, P.Z, middle0);//XP+ZP
	nres_modsub(P.X, P.Z, middle1);//XP-ZP
	nres_modmult(middle0, middle0, t0);//t0=(XP+ZP)^2
	nres_modmult(middle1, middle1, t1);//t1=(XP-ZP)^2
	nres_modmult(t0, t1, P2.X);//P2.X=t0t1
	nres_modsub(t0, t1, t2);//t2=t0-t1
	nres_modmult(A24, t2,t0);//t0=A24 t2
	nres_modadd(t0, t1, middle2);//middle2=t0+t1
	nres_modmult(t2, middle2, P2.Z);//P2.Z=t2(t0+t1)
	return 0;
}

int cADD(MontgomeryPoint PminusQ, MontgomeryPoint P, MontgomeryPoint Q, MontgomeryPoint PaddQ, big A)
{//PminusQ(X_(p-q)^{-1},Z_{P-Q}^{-1});P(XP,ZP):Q(XQ,ZQ) 
	big A24 = mirvar(1);
	big two = mirvar(2);
	big four = mirvar(4);
	big t0 = mirvar(2);
	big t1 = mirvar(4);
	big t2 = mirvar(4);
	big middle0 = mirvar(2);
	big middle1 = mirvar(4);
	big middle2 = mirvar(4);

 

	nres_modadd(Q.X, Q.Z, middle0);//XQ+ZQ
	nres_modsub(P.X, P.Z, middle1);//XP-ZP
	nres_modmult(middle0, middle1, t0);//t0=(XP-ZP)(XQ+ZQ)

	nres_modsub(Q.X, Q.Z, middle0);//XQ-ZQ
	nres_modadd(P.X, P.Z, middle1);//XP+ZP
	nres_modmult(middle0, middle1, t1);//t1=(XP+ZP)(XQ-ZQ)


	nres_modadd(t0, t1, middle0);//t0+t1	
	nres_modmult(middle0, middle0, middle1);//middle1=(t0+t1)^2
	nres_modmult(PminusQ.X, middle1, PaddQ.X);//PaddQ.X=X_(p-q)^{-1}(t0+t1)^2

	nres_modsub(t0, t1, middle0);//t0-t1	
	nres_modmult(middle0, middle0, middle1);//middle1=(t0-t1)^2
	nres_modmult(PminusQ.Z, middle1, PaddQ.Z);//PaddQ.Z=Z_(p-q)^{-1}(t0+t1)^2

	nres_moddiv(PaddQ.X, four, PaddQ.X);//PaddQ.X=PaddQ.X / 4;
	nres_moddiv(PaddQ.Z, four, PaddQ.Z);//PaddQ.X=PaddQ.Z / 4;
	return 0;
}

int copy(MontgomeryPoint S, MontgomeryPoint T)
{
	copy(S.X, T.X);
	copy(S.Z, T.Z);
	return 0;
}

int cLadder(Big nScalar, MontgomeryPoint PminusQ, MontgomeryPoint P, MontgomeryPoint Q, MontgomeryPoint nP, MontgomeryPoint nPaddQ, big A)
{
	MontgomeryPoint S0;
	MontgomeryPoint S1;
	MontgomeryPoint T;
	MontgomeryPoint R;
	MontgomeryPoint invQ;
	S0.X = mirvar(0);
	S0.Z = mirvar(1);
	S1.X = mirvar(0);
	S1.Z = mirvar(1);
	T.X = mirvar(0);
	T.Z = mirvar(1);
	R.X = mirvar(0);
	R.Z = mirvar(1);

	big one = mirvar(0);
	invQ.X = mirvar(0);
	invQ.Z = mirvar(1);

	copy(P.X, S1.X);
	copy(P.Z, S1.Z);

	copy(Q.X, T.X);
	copy(Q.Z, T.Z);

	nres_moddiv(one, Q.X, invQ.X);
	nres_moddiv(one, Q.Z, invQ.Z);

	int b = 1024;
	int n[1024];
	binary(nScalar, n, b);

	for (int i = b - 1; i > -1; i--)
	{

		cADD(PminusQ, S0, S1, R, A);
		if (n[i] == n[i + 1])
		{
			swap(S0, S1); 
			swap(PminusQ, invQ);
		}

		cADD(invQ, T, S0, T, A);

		cDBL(S0, S0, A);
		copy(R, S1);
	}

	copy(S0, nP);
	copy(T, nPaddQ);
	return 0;
}


int traslate(MontgomeryPoint P, MontgomeryPoint T, MontgomeryPoint Q, big A)
{


	big r = mirvar(1);
	big s = mirvar(2);
	big four = mirvar(4);
	big t0 = mirvar(2);
	big t1 = mirvar(4);
	big t2 = mirvar(4);
	big middle0 = mirvar(2);
	big middle1 = mirvar(4);
	big middle2 = mirvar(4);


	copy(T.X, r);
	copy(T.Z, s);

	nres_modmult(P.X, r, middle1);//rXP 
	nres_modmult(P.Z, s, middle0);//sZP	
	nres_modsub(middle1, middle0, Q.X);//rXP-sZP

	nres_modmult(P.X, s, middle1);//sXP 
	nres_modmult(P.Z, r, middle0);//rZP	
	nres_modsub(middle1, middle0, Q.Z);//sXP-rZP

	Big zero = 0;
	zero = s;
	if (zero == 0)
	{
		copy(P.X, Q.X);
		copy(P.Z, Q.Z);
	}
	zero = r;
	if (zero == 0)
	{
		copy(P.Z, Q.X);
		copy(P.X, Q.Z);
	}
	return 0;
}

int TatePairing(big l, MontgomeryPoint P, MontgomeryPoint Q, big eTl, big A)
{
	MontgomeryPoint invPminusQ;
	MontgomeryPoint negativeQ;
	invPminusQ.X = mirvar(0);
	invPminusQ.Z = mirvar(1);
	negativeQ.X = mirvar(0);
	negativeQ.Z = mirvar(1);
	big zero = mirvar(0);
	big one = mirvar(1);
	affinePoint P0;
	affinePoint Q0;

	P0.x = mirvar(0);
	P0.y = mirvar(1);
	Q0.x = mirvar(0);
	Q0.y = mirvar(1);
	ConstructMontgomeryPoint(A, P0);
	ConstructMontgomeryPoint(A, Q0);
	MontgomeryAddition(P0,Q0, invPminusQ, A);

	copy(P0.x, P.X);
	copy(Q0.x, Q.X);
	copy(one, P.Z);
	copy(one, Q.Z);

	MontgomeryPoint nP;
	MontgomeryPoint nPaddQ;
	MontgomeryPoint lP;
	MontgomeryPoint lPaddQ;

	nP.X = mirvar(0);
	nP.Z = mirvar(1);
	nPaddQ.X = mirvar(0);
	nPaddQ.Z = mirvar(1);

	lP.X = mirvar(0);
	lP.Z = mirvar(1);
	lPaddQ.X = mirvar(0);
	lPaddQ.Z = mirvar(1);
	Big n = l;
	int flag = 0;
	if (n % 2 == 0)
	{
		n = n / 2;
		flag = 1;
	}
	//int cLadder(Big nScalar, MontgomeryPoint PminusQ, MontgomeryPoint P, MontgomeryPoint Q, MontgomeryPoint nP, MontgomeryPoint nPaddQ, big A)
	cLadder(n, invPminusQ, P,Q,nP, nPaddQ,A);

	traslate(nPaddQ, nP, lPaddQ, A);
	traslate(nP, nP, lP, A);

	big lambdaQ = mirvar(0);
	big lambdaP = mirvar(1);
	copy(lPaddQ.Z, lambdaQ);
	copy(lP.X, lambdaP);


	nres_moddiv(lambdaQ, lambdaP, eTl);
	return 0;

}
int MontgomeryLadder(big lTate, big number, big &Q, big ConstantA)
{
	nres_powmod(lTate, number, Q);
	return 0;
}
int checkOrder(big lTate, int &isSupersingular)
{//lTate is FirstPart order
	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };	 
	Big FirstPart = 1;
	isSupersingular = 1;
	for (int i = 0; i < 33; i++)
	{
		FirstPart = FirstPart*l[i];
	}
	

	big P[65];
	for (int i = 0; i < 65; i++)
	{
		P[i]=mirvar(1);
	}
	big constantA = mirvar(0);
	Big number = 1;
	
	//Level 1
	number = 1;
	for (int i = 16; i < 34; i++)
	{
		number = number*l[i];
	}
	
	MontgomeryLadder(P[1], number.getbig(), P[2], constantA);
	//P[2]=n/(l[0]l[1]\ldotsl[115])P

	number = 1;
	for (int i = 0; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[1], number.getbig(), P[3], constantA);
	//P[3]=n/(l[16]l[17]\ldotsl[33])P



	//Level 2
	number = 1;
	for (int i = 8; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[2], number.getbig(), P[4], constantA);
	//P[4]=n/(l[0]l[1]\ldotsl[7])P

	number = 1;
	for (int i = 0; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[2], number.getbig(), P[5], constantA);
	//P[5]=n/(l[8]l[9]\ldotsl[15])P



	number = 1;
	for (int i = 24; i < 33; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[3], number.getbig(), P[6], constantA);
	//P[6]=n/(l[16]l[17]\ldotsl[23])P

	number = 1;
	for (int i = 16; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[3], number.getbig(), P[7], constantA);
	//P[7]=n/(l[24]l[25]\ldotsl[32])P

	//Level 3
	number = 1;
	for (int i = 4; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[4], number.getbig(), P[8], constantA);
	//P[8]=n/(l[0]l[1]\ldotsl[3])P

	number = 1;
	for (int i = 0; i < 4; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[4], number.getbig(), P[9], constantA);
	//P[9]=n/(l[4]l[5]\ldotsl[7])P

	number = 1;
	for (int i = 12; i < 15; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[5], number.getbig(), P[10], constantA);
	//P[10]=n/(l[8]l[9]\ldotsl[11])P
	number = 1;
	for (int i = 8; i < 12; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[5], number.getbig(), P[11], constantA);
	//P[11]=n/(l[12]l[13]\ldotsl[15])P


	number = 1;
	for (int i = 20; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[6], number.getbig(), P[12], constantA);
	//P[12]=n/(l[16]l[17]\ldotsl[19])P
	number = 1;
	for (int i = 16; i < 20; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[6], number.getbig(), P[13], constantA);
	//P[13]=n/(l[20]l[21]\ldotsl[23])P

	number = 1;
	for (int i = 28; i < 32; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[7], number.getbig(), P[14], constantA);
	//P[14]=n/(l[24]l[25]\ldotsl[27])P

	number = 1;
	for (int i = 24; i < 28; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[7], number.getbig(), P[15], constantA);
	//P[15]=n/(l[28]l[29]\ldotsl[32])P

	//////////////////////
	//Level 4
	number = 1;
	for (int i = 2; i < 4; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[8], number.getbig(), P[16], constantA);
	//P[16]=n/(l[0]l[1])P

	number = 1;
	for (int i = 0; i < 2; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[8], number.getbig(), P[17], constantA);
	//P[17]=n/(l[2]l[3])P



	number = 1;
	for (int i = 6; i < 8; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[9], number.getbig(), P[18], constantA);
	//P[18]=n/(l[4]l[5])P

	number = 1;
	for (int i = 4; i < 6; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[9], number.getbig(), P[19], constantA);
	//P[19]=n/(l[6]l[7])P

	number = 1;
	for (int i = 10; i < 12; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[10], number.getbig(), P[20], constantA);
	//P[20]=n/(l[8]l[9])P
	number = 1;
	for (int i = 8; i < 10; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[10], number.getbig(), P[21], constantA);
	//P[21]=n/(l[10]l[11])P

	number = 1;
	for (int i = 14; i < 16; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[11], number.getbig(), P[22], constantA);
	//P[22]=n/(l[12]l[13]\ldotsl[15])P
	number = 1;
	for (int i = 12; i < 14; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[11], number.getbig(), P[23], constantA);
	//P[23]=n/(l[14]l[15]\ldotsl[15])P


	number = 1;
	for (int i = 18; i < 20; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[12], number.getbig(), P[24], constantA);
	//P[24]=n/(l[16]l[17])P
	number = 1;
	for (int i = 16; i < 18; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[12], number.getbig(), P[25], constantA);
	//P[25]=n/(l[18]l[19])P

	number = 1;
	for (int i = 22; i < 24; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[13], number.getbig(), P[26], constantA);
	//P[26]=n/(l[20]l[21])P
	number = 1;
	for (int i = 20; i < 22; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[13], number.getbig(), P[27], constantA);
	//P[27]=n/(l[22]l[23])P

	number = 1;
	for (int i = 26; i < 28; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[14], number.getbig(), P[28], constantA);
	//P[28]=n/(l[24]l[25])P

	number = 1;
	for (int i = 24; i < 26; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[14], number.getbig(), P[29], constantA);
	//P[29]=n/(l[26]l[27])P

	number = 1;
	for (int i = 30; i < 33; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[15], number.getbig(), P[30], constantA);
	//P[30]=n/(l[28]l[29]\ldotsl[32])P

	number = 1;
	for (int i = 28; i < 30; i++)
	{
		number = number*l[i];
	}
	MontgomeryLadder(P[15], number.getbig(), P[31], constantA);
	//P[31]=n/(l[30]l[31][32])P

	//////////////////////
	//Level 5
	number = 1;
	number = l[1];
	MontgomeryLadder(P[16], number.getbig(), P[32], constantA);
	//P[32]=n/(l[0])P

	number = 1;
	number = number*l[0];
	MontgomeryLadder(P[16], number.getbig(), P[33], constantA);
	//P[33]=n/(l[1])P

	number = 1;
	number = l[3];
	MontgomeryLadder(P[17], number.getbig(), P[34], constantA);
	//P[34]=n/(l[2])P

	number = 1;
	number = number*l[2];
	MontgomeryLadder(P[17], number.getbig(), P[35], constantA);
	//P[35]=n/(l[3])P



	number = 1;
	number = l[5];
	MontgomeryLadder(P[18], number.getbig(), P[36], constantA);
	//P[36]=n/(l[4])P

	number = 1;
	number = number*l[4];
	MontgomeryLadder(P[18], number.getbig(), P[37], constantA);
	//P[37]=n/(l[5])P

	number = 1;
	number = l[7];
	MontgomeryLadder(P[19], number.getbig(), P[38], constantA);
	//P[38]=n/(l[6])P

	number = 1;
	number = number*l[6];
	MontgomeryLadder(P[19], number.getbig(), P[39], constantA);
	//P[39]=n/(l[7])P

	number = 1;
	number = l[9];
	MontgomeryLadder(P[20], number.getbig(), P[40], constantA);
	//P[40]=n/(l[8])P

	number = 1;
	number = number*l[8];
	MontgomeryLadder(P[20], number.getbig(), P[41], constantA);
	//P[41]=n/(l[9])P

	number = 1;
	number = l[11];
	MontgomeryLadder(P[21], number.getbig(), P[42], constantA);
	//P[42]=n/(l[10])P

	number = 1;
	number = number*l[10];
	MontgomeryLadder(P[21], number.getbig(), P[43], constantA);
	//P[43]=n/(l[11])P

	number = 1;
	number = l[13];
	MontgomeryLadder(P[22], number.getbig(), P[44], constantA);
	//P[42]=n/(l[12])P

	number = 1;
	number = number*l[12];
	MontgomeryLadder(P[22], number.getbig(), P[45], constantA);
	//P[43]=n/(l[13])P

	number = 1;
	number = l[15];
	MontgomeryLadder(P[23], number.getbig(), P[46], constantA);
	//P[44]=n/(l[14])P

	number = 1;
	number = number*l[14];
	MontgomeryLadder(P[23], number.getbig(), P[47], constantA);
	//P[45]=n/(l[15])P

	number = 1;
	number = l[17];
	MontgomeryLadder(P[24], number.getbig(), P[48], constantA);
	//P[48]=n/(l[16])P

	number = 1;
	number = number*l[16];
	MontgomeryLadder(P[24], number.getbig(), P[49], constantA);
	//P[49]=n/(l[17])P

	number = 1;
	number = l[19];
	MontgomeryLadder(P[25], number.getbig(), P[50], constantA);
	//P[50]=n/(l[18])P

	number = 1;
	number = number*l[18];
	MontgomeryLadder(P[25], number.getbig(), P[51], constantA);
	//P[51]=n/(l[19])P


	number = 1;
	number = l[21];
	MontgomeryLadder(P[26], number.getbig(), P[52], constantA);
	//P[52]=n/(l[20])P

	number = 1;
	number = number*l[20];
	MontgomeryLadder(P[26], number.getbig(), P[53], constantA);
	//P[53]=n/(l[21])P

	number = 1;
	number = l[23];
	MontgomeryLadder(P[27], number.getbig(), P[54], constantA);
	//P[54]=n/(l[22])P

	number = 1;
	number = number*l[22];
	MontgomeryLadder(P[27], number.getbig(), P[55], constantA);
	//P[55]=n/(l[23])P

	number = 1;
	number = l[25];
	MontgomeryLadder(P[28], number.getbig(), P[56], constantA);
	//P[56]=n/(l[24])P

	number = 1;
	number = number*l[24];
	MontgomeryLadder(P[28], number.getbig(), P[57], constantA);
	//P[57]=n/(l[25])P

	number = 1;
	number = l[27];
	MontgomeryLadder(P[29], number.getbig(), P[58], constantA);
	//P[58]=n/(l[26])P

	number = 1;
	number = number*l[26];
	MontgomeryLadder(P[29], number.getbig(), P[59], constantA);
	//P[59]=n/(l[27])P


	number = 1;
	number = l[29];
	MontgomeryLadder(P[30], number.getbig(), P[60], constantA);
	//P[60]=n/(l[28])P

	number = 1;
	number = number*l[28];
	MontgomeryLadder(P[30], number.getbig(), P[61], constantA);
	//P[61]=n/(l[29])P

	number = 1;
	number = l[31] * l[32];
	MontgomeryLadder(P[31], number.getbig(), P[62], constantA);
	//P[62]=n/(l[30])P

	number = 1;
	number = l[30] * l[32];
	MontgomeryLadder(P[31], number.getbig(), P[63], constantA);
	//P[63]=n/(l[31])P

	number = 1;
	number = l[30] * l[31];
	MontgomeryLadder(P[31], number.getbig(), P[64], constantA);
	//P[64]=n/(l[32])P
	 
	Big x = P[0];
	for (int i = 32; i < 64; i++)
	{
		x = P[i];
		if (x == 1)
		{
			isSupersingular = 0;
			return -1;
		}
	}
	number = l[32];
	MontgomeryLadder(P[64], number.getbig(), P[0], constantA);
	x = P[0];
	if (x != 1)
	{
		isSupersingular = 0;
		return -1;
	}
	return 0;
}
int pairingMontgomeryLadder(big A,  int &IsSupersingularity)
{

	int l[74] = { 179, 239, 313, 337, 197, 251, 263, 347, 227, 257, 283, 293, 241, 233, 277, 311, 269, 317, 193, 223,
		181, 211, 281, 353, 191, 229, 359, 349, 199, 271, 307, 331,	367, 373, 587, 173, 167, 163, 157, 151,
		149, 139, 137, 131, 127, 113, 109, 107, 103, 101, 97, 89, 83, 79, 73, 71, 67, 61, 59, 53,
		47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3 };
	big eT1 = mirvar(0);
	Big FirstPart = 1;
	Big SecondPart = 1;


	for (int i = 0; i < 33; i++)
	{
		FirstPart = FirstPart*l[i];
	}

	for (int i = 33; i < 74; i++)
	{
		SecondPart = SecondPart*l[i];
	}

	SecondPart = SecondPart * 4;
	big lTate = mirvar(1);
	copy(SecondPart.getbig(), lTate);
	MontgomeryPoint P;
	MontgomeryPoint Q;
	P.X = mirvar(0);
	P.Z = mirvar(1);
	Q.X = mirvar(0);
	Q.Z = mirvar(1);
	TatePairing(lTate, P, Q, lTate, A);
	//Check the order of lTate
	checkOrder(lTate, IsSupersingularity);
	return 0;
}