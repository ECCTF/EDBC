#include "preComputation.h"

//int MontgomeryTrick(projectPoint *P, int Cmax)
//{//montgomery Trick
//
//	big a[33];
//	big b[33];
//	big c[33];
//	int n = Cmax;
//	big one = mirvar(1), d = mirvar(1), R1 = mirvar(0);
//	for (int i = 0; i < 33; i++)
//	{
//		a[i] = mirvar(0);
//		b[i] = mirvar(0);
//		c[i] = mirvar(1);
//
//	}
//	for (int i = 1; i < n + 1; i++)
//	{
//		copy(P[2 * i + 1].Z, a[i]);
//	}
//	//a[i]=3P,5P,7P,9P,....
//	copy(a[1], c[1]);
//	//c1=a1
//
//	for (int i = 2; i < n + 1; i++)
//	{
//		modmult(a[i], c[i - 1], c[i]);
//	}
//	//ci=ci-1*ai
//
//	inverse(c[n], d);
//
//	for (int i = n; i > 1; i--)
//	{
//		modmult(d, c[i - 1], b[i]);
//		modmult(a[i], d, d);
//	}
//	copy(d, b[1]);
//
//
//	for (int i = 1; i < n + 1; i++)
//	{
//		modmult2(b[i], P[2 * i + 1].X, P[2 * i + 1].X);
//		modsquare2(b[i], R1);
//		modmult2(R1, P[2 * i + 1].Y, P[2 * i + 1].Y);
//		copy(one, P[2 * i + 1].Z);
//	}
//
//	//Big out=0;
//	//for(int i=1;i<n+1;i++)
//	//{	
//	//	cout <<"Q_"<<2*i+1<<endl;
//	//	out =P[2*i+1].X;
//	//	cout <<"X:"<<out <<endl;
//	//	out =P[2*i+1].Y;
//	//	cout <<"Y:"<<out <<endl;
//	//	out =P[2*i+1].Z;
//	//	cout <<"Z:"<<out <<endl;
//	//	 
//	//}
//
//	return 0;
//}



int preComputationW(affinePoint P1, projectPoint *P, long long &CPUcycles, int Cmax)
{//P[i]=i*P, a=-3
	
	big one = mirvar(1);
	copy(P1.x, P[1].X);
	copy(P1.y, P[1].Y);
	copy(one, P[1].Z);
	switch (Cmax)
	{
	case 3:
		JacobiTriple(P[1], P[3]);
		break;
	case 5:
		JacobiDouble(P[1], P[2]);
		JacobiDouble(P[2], P[4]);
		JacobiAdditionSubtraction(P[4], P[1], P[5], P[3]);
		break;
	case 7:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		break;
	case 9:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiTriple(P[3], P[9]);
		break;
	case 11:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[5], P[10]);
		JacobiAdditionSubtraction(P[10], P[1], P[11], P[9]);
		break;
	case 13:
		JacobiDouble(P[1], P[2]);
		JacobiDouble(P[2], P[4]);
		JacobiAdditionSubtraction(P[4], P[1], P[5], P[3]);
		JacobiDouble(P[5], P[10]);
		JacobiAdditionSubtraction(P[10], P[1], P[11], P[9]);
		JacobiAdditionSubtraction(P[10], P[3], P[13], P[7]);
		break;
	case 15:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[6], P[12]);
		JacobiAdditionSubtraction(P[12], P[1], P[13], P[11]);
		JacobiAdditionSubtraction(P[12], P[3], P[15], P[9]);
		break;
	case 17:
		JacobiTriple(P[1], P[3]);
		JacobiTriple(P[3], P[9]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[7], P[14]);
		JacobiAdditionSubtraction(P[14], P[1], P[15], P[13]);
		JacobiAdditionSubtraction(P[14], P[3], P[17], P[11]);
		break;
	case 19:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[7], P[14]);
		JacobiAdditionSubtraction(P[14], P[1], P[15], P[13]);
		JacobiAdditionSubtraction(P[14], P[3], P[17], P[11]);
		JacobiAdditionSubtraction(P[14], P[5], P[19], P[9]);
		break;
	case 21:
		JacobiDouble(P[1], P[2]);
		JacobiDouble(P[2], P[4]);
		JacobiAdditionSubtraction(P[4], P[1], P[5], P[3]);
		JacobiDouble(P[5], P[10]);
		JacobiAdditionSubtraction(P[10], P[1], P[11], P[9]);
		JacobiAdditionSubtraction(P[10], P[3], P[13], P[7]);
		JacobiDouble(P[9], P[18]);
		JacobiAdditionSubtraction(P[18], P[1], P[19], P[17]);
		JacobiAdditionSubtraction(P[18], P[3], P[21], P[15]);
		break;
	case 23:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[5], P[10]);
		JacobiAdditionSubtraction(P[10], P[1], P[11], P[9]);
		JacobiDouble(P[9], P[18]);
		JacobiAdditionSubtraction(P[18], P[1], P[19], P[17]);
		JacobiAdditionSubtraction(P[18], P[3], P[21], P[15]);
		JacobiAdditionSubtraction(P[18], P[5], P[23], P[13]);
		break;
	case 25:
		JacobiDouble(P[1], P[2]);
		JacobiDouble(P[2], P[4]);
		JacobiAdditionSubtraction(P[4], P[1], P[5], P[3]);
		JacobiDouble(P[5], P[10]);
		JacobiAdditionSubtraction(P[10], P[1], P[11], P[9]);
		JacobiAdditionSubtraction(P[10], P[3], P[13], P[7]);
		JacobiDouble(P[10], P[20]);
		JacobiAdditionSubtraction(P[20], P[1], P[21], P[19]);
		JacobiAdditionSubtraction(P[20], P[3], P[23], P[17]);
		JacobiAdditionSubtraction(P[20], P[5], P[25], P[15]);
		break;
	case 27:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[6], P[12]);
		JacobiAdditionSubtraction(P[12], P[1], P[13], P[11]);
		JacobiAdditionSubtraction(P[12], P[3], P[15], P[9]);
		JacobiDouble(P[11], P[22]);
		JacobiAdditionSubtraction(P[22], P[1], P[23], P[21]);
		JacobiAdditionSubtraction(P[22], P[3], P[25], P[19]);
		JacobiAdditionSubtraction(P[22], P[5], P[27], P[17]);
		break;
	case 29:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[6], P[12]);
		JacobiAdditionSubtraction(P[12], P[1], P[13], P[11]);
		JacobiAdditionSubtraction(P[12], P[3], P[15], P[9]);
		JacobiAdditionSubtraction(P[12], P[5], P[17], P[7]);
		JacobiDouble(P[12], P[24]);
		JacobiAdditionSubtraction(P[24], P[1], P[25], P[23]);
		JacobiAdditionSubtraction(P[24], P[3], P[27], P[21]);
		JacobiAdditionSubtraction(P[24], P[5], P[29], P[19]);
		break;
	case 31:
		JacobiTriple(P[1], P[3]);
		JacobiDouble(P[3], P[6]);
		JacobiAdditionSubtraction(P[6], P[1], P[7], P[5]);
		JacobiDouble(P[6], P[12]);
		JacobiAdditionSubtraction(P[12], P[1], P[13], P[11]);
		JacobiAdditionSubtraction(P[12], P[3], P[15], P[9]);
		JacobiDouble(P[12], P[24]);
		JacobiAdditionSubtraction(P[24], P[1], P[25], P[23]);
		JacobiAdditionSubtraction(P[24], P[3], P[27], P[21]);
		JacobiAdditionSubtraction(P[24], P[5], P[29], P[19]);
		JacobiAdditionSubtraction(P[24], P[7], P[31], P[17]);
		break;
	default:
		cout << "Can not deal!" << endl;//不能处理的情况
	}
	return 0;
}



int preComputationE(affinePoint P1, projectPoint *P, long long &CPUcycles, int Cmax, big d,big a)
{
	big one = mirvar(1);
	copy(P1.x, P[1].X);
	copy(P1.y, P[1].Y);
	copy(one, P[1].Z);
	switch (Cmax)
	{
	case 3:
		InvEdwardsTriple(P[1], P[3],d);
		break;
	case 5:
		InvEdwardsDouble(P[1], P[2], d);
		InvEdwardsAddition(P[2], P[1], P[3], d, a);
		InvEdwardsAddition(P[2], P[3], P[5], d, a);
		break;
	case 7:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		break;
	case 9:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsAddition(P[6], P[3], P[9], d, a);
		break;
	case 11:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[5], P[10], d);
		InvEdwardsAdditionSubtraction(P[10], P[1], P[11], P[9], d, a);
		break;
	case 13:
		InvEdwardsDouble(P[1], P[2], d);
		InvEdwardsAddition(P[2], P[1], P[3], d, a);
		InvEdwardsAddition(P[2], P[3], P[5], d, a);
		InvEdwardsDouble(P[5], P[10], d);
		InvEdwardsAdditionSubtraction(P[10], P[1], P[11], P[9], d, a);
		InvEdwardsAdditionSubtraction(P[10], P[3], P[13], P[7], d, a);
		break;
	case 15:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[6], P[12], d);
		InvEdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d, a);
		break;
	case 17:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[7], P[14], d);
		InvEdwardsAdditionSubtraction(P[14], P[1], P[15], P[13], d, a);
		InvEdwardsTriple(P[3], P[9], d);
		InvEdwardsAdditionSubtraction(P[14], P[3], P[17], P[11], d, a);
		break;
	case 19:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[7], P[14], d);
		InvEdwardsAdditionSubtraction(P[14], P[1], P[15], P[13], d, a);
		InvEdwardsAdditionSubtraction(P[14], P[3], P[17], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[14], P[5], P[19], P[9], d, a);
		break;
	case 21:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[7], P[14], d);
		InvEdwardsAdditionSubtraction(P[14], P[1], P[15], P[13], d, a);
		InvEdwardsAdditionSubtraction(P[14], P[3], P[17], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[14], P[5], P[19], P[9], d, a);
		InvEdwardsAddition(P[14], P[7], P[21], d, a);
		break;
	case 23:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[5], P[10], d);
		InvEdwardsAdditionSubtraction(P[10], P[1], P[11], P[9], d, a);
		InvEdwardsDouble(P[9], P[18], d);
		InvEdwardsAdditionSubtraction(P[18], P[1], P[19], P[17], d, a);
		InvEdwardsAdditionSubtraction(P[18], P[3], P[21], P[15], d, a);
		InvEdwardsAdditionSubtraction(P[18], P[5], P[23], P[13], d, a);
		break;
	case 25:
		InvEdwardsDouble(P[1], P[2], d);
		InvEdwardsAddition(P[2], P[1], P[3], d, a);
		InvEdwardsAddition(P[2], P[3], P[5], d, a);
		InvEdwardsDouble(P[5], P[10], d);
		InvEdwardsAdditionSubtraction(P[10], P[1], P[11], P[9], d, a);
		InvEdwardsAdditionSubtraction(P[10], P[3], P[13], P[7], d, a);
		InvEdwardsDouble(P[10], P[20], d);
		InvEdwardsAdditionSubtraction(P[20], P[1], P[21], P[19], d, a);
		InvEdwardsAdditionSubtraction(P[20], P[3], P[23], P[17], d, a);
		InvEdwardsAdditionSubtraction(P[20], P[5], P[25], P[15], d, a);
		break;
	case 27:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[6], P[12], d);
		InvEdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d, a);
		InvEdwardsDouble(P[11], P[22], d);
		InvEdwardsAdditionSubtraction(P[22], P[1], P[23], P[21], d, a);
		InvEdwardsAdditionSubtraction(P[22], P[3], P[25], P[19], d, a);
		InvEdwardsAdditionSubtraction(P[22], P[5], P[27], P[17], d, a);
		break;
	case 29:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[6], P[12], d);
		InvEdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d, a);
		InvEdwardsAddition(P[12], P[5], P[17], d, a);
		InvEdwardsDouble(P[12], P[24], d);
		InvEdwardsAdditionSubtraction(P[24], P[1], P[25], P[23], d, a);
		InvEdwardsAdditionSubtraction(P[24], P[3], P[27], P[21], d, a);
		InvEdwardsAdditionSubtraction(P[24], P[5], P[29], P[19], d, a);
		break;
	case 31:
		InvEdwardsTriple(P[1], P[3], d);
		InvEdwardsDouble(P[3], P[6], d);
		InvEdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d, a);
		InvEdwardsDouble(P[6], P[12], d);
		InvEdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d, a);
		InvEdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d, a);
		InvEdwardsDouble(P[12], P[24], d);
		InvEdwardsAdditionSubtraction(P[24], P[1], P[25], P[23], d, a);
		InvEdwardsAdditionSubtraction(P[24], P[3], P[27], P[21], d, a);
		InvEdwardsAdditionSubtraction(P[24], P[5], P[29], P[19], d, a);
		InvEdwardsAdditionSubtraction(P[24], P[7], P[31], P[17], d, a);
		break;
	default:
		cout << "Can not deal!" << endl;//不能处理的情况
	}
	return 0;
}



int preComputationD(affinePoint P1, extProjectPoint *P, long long &CPUcycles, int Cmax)
{
	big one = mirvar(1);
	copy(P1.x, P[1].X);
	copy(P1.y, P[1].Y);
	copy(one, P[1].Z);
	copy(one, P[1].T);

	switch(Cmax)
	{
	case 3:
		DIKTripling(P[1], P[3]);
		break;
	case 5:
		DIKDoubling(P[1], P[2]);
		DIKDoubling(P[2], P[4]);
		DIKAdditionSubtraction(P[4], P[1], P[5], P[3]);
		break;
	case 7:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		break;
	case 9:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKTripling(P[3], P[9]);
		break;
	case 11:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[5], P[10]);
		DIKAdditionSubtraction(P[10], P[1], P[11], P[9]);
		break;
	case 13:
		DIKDoubling(P[1], P[2]);
		DIKDoubling(P[2], P[4]);
		DIKAdditionSubtraction(P[4], P[1], P[5], P[3]);
		DIKDoubling(P[5], P[10]);
		DIKAdditionSubtraction(P[10], P[1], P[11], P[9]);
		DIKAdditionSubtraction(P[10], P[3], P[13], P[7]);
		break;
	case 15:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[6], P[12]);
		DIKAdditionSubtraction(P[12], P[1], P[13], P[11]);
		DIKAdditionSubtraction(P[12], P[3], P[15], P[9]);
		break;
	case 17:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[7], P[14]);
		DIKAdditionSubtraction(P[14], P[1], P[15], P[13]);
		DIKAdditionSubtraction(P[14], P[3], P[17], P[11]);
		DIKTripling(P[3], P[9]);
		break;
	case 19:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[7], P[14]);
		DIKAdditionSubtraction(P[14], P[1], P[15], P[13]);
		DIKAdditionSubtraction(P[14], P[3], P[17], P[11]);
		DIKAdditionSubtraction(P[14], P[5], P[19], P[9]);
		break;
	case 21:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[7], P[14]);
		DIKAdditionSubtraction(P[14], P[1], P[15], P[13]);
		DIKAdditionSubtraction(P[14], P[3], P[17], P[11]);
		DIKAdditionSubtraction(P[14], P[5], P[19], P[9]);
		DIKAddition(P[14], P[7], P[21]);
		break;
	case 23:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[5], P[10]);
		DIKAdditionSubtraction(P[10], P[1], P[11], P[9]);
		DIKDoubling(P[9], P[18]);
		DIKAdditionSubtraction(P[18], P[1], P[19], P[17]);
		DIKAdditionSubtraction(P[18], P[3], P[21], P[15]);
		DIKAdditionSubtraction(P[18], P[5], P[23], P[13]);
		break;
	case 25:
		DIKDoubling(P[1], P[2]);
		DIKDoubling(P[2], P[4]);
		DIKAdditionSubtraction(P[4], P[1], P[5], P[3]);
		DIKDoubling(P[5], P[10]);
		DIKAdditionSubtraction(P[10], P[1], P[11], P[9]);
		DIKAdditionSubtraction(P[10], P[3], P[13], P[7]);
		DIKDoubling(P[10], P[20]);
		DIKAdditionSubtraction(P[20], P[1], P[21], P[19]);
		DIKAdditionSubtraction(P[20], P[3], P[23], P[17]);
		DIKAdditionSubtraction(P[20], P[5], P[25], P[15]);
		break;
	case 27:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[6], P[12]);
		DIKAdditionSubtraction(P[12], P[1], P[13], P[11]);
		DIKAdditionSubtraction(P[12], P[3], P[15], P[9]);
		DIKDoubling(P[11], P[22]);
		DIKAdditionSubtraction(P[22], P[1], P[23], P[21]);
		DIKAdditionSubtraction(P[22], P[3], P[25], P[19]);
		DIKAdditionSubtraction(P[22], P[5], P[27], P[17]);
		break;
	case 29:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[6], P[12]);
		DIKAdditionSubtraction(P[12], P[1], P[13], P[11]);
		DIKAdditionSubtraction(P[12], P[3], P[15], P[9]);
		DIKAddition(P[12], P[5], P[17]);
		DIKDoubling(P[12], P[24]);
		DIKAdditionSubtraction(P[24], P[1], P[25], P[23]);
		DIKAdditionSubtraction(P[24], P[3], P[27], P[21]);
		DIKAdditionSubtraction(P[24], P[5], P[29], P[19]);
		break;
	case 31:
		DIKTripling(P[1], P[3]);
		DIKDoubling(P[3], P[6]);
		DIKAdditionSubtraction(P[6], P[1], P[7], P[5]);
		DIKDoubling(P[6], P[12]);
		DIKAdditionSubtraction(P[12], P[1], P[13], P[11]);
		DIKAdditionSubtraction(P[12], P[3], P[15], P[9]);
		DIKDoubling(P[12], P[24]);
		DIKAdditionSubtraction(P[24], P[1], P[25], P[23]);
		DIKAdditionSubtraction(P[24], P[3], P[27], P[21]);
		DIKAdditionSubtraction(P[24], P[5], P[29], P[19]);
		DIKAdditionSubtraction(P[24], P[7], P[31], P[17]);
		break;
	default:
		cout << "Can not deal!" << endl;//不能处理的情况
	}
	
	return 0;
}



int preComputationEx(affinePoint P1, extProjectPoint *P, long long &CPUcycles, int Cmax, big d)
{
	big one = mirvar(1);
	copy(P1.x, P[1].X);
	copy(P1.y, P[1].Y);
	copy(one, P[1].Z);
	nres_modmult(P[1].X, P[1].Y, P[1].T);

	switch (Cmax)
	{
	case 3:
		EdwardsTriple(P[1], P[3], d);
		break;
	case 5:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[2], P[3], P[5], d);
		break;
	case 7:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		break;
	case 9:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		EdwardsAddition(P[7], P[2], P[9], d);
		break;
	case 11:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		EdwardsAddition(P[7], P[2], P[9], d);
		EdwardsAddition(P[9], P[2], P[11], d);
		break;
	case 13:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		EdwardsAddition(P[7], P[2], P[9], d);
		EdwardsAddition(P[9], P[2], P[11], d);
		EdwardsAddition(P[11], P[2], P[13], d);
		break;
	case 15:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		EdwardsAddition(P[7], P[2], P[9], d);
		EdwardsAddition(P[9], P[2], P[11], d);
		EdwardsAddition(P[11], P[2], P[13], d);
		EdwardsAddition(P[13], P[2], P[15], d);
		break;
	case 17:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsAddition(P[5], P[2], P[7], d);
		EdwardsAddition(P[7], P[2], P[9], d);
		EdwardsAddition(P[9], P[2], P[11], d);
		EdwardsAddition(P[11], P[2], P[13], d);
		EdwardsAddition(P[13], P[2], P[15], d);
		EdwardsAddition(P[15], P[2], P[17], d);
		break;
	case 19:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsDouble(P[7], P[14], d);
		EdwardsAdditionSubtraction(P[14], P[1], P[15], P[13], d);
		EdwardsAdditionSubtraction(P[14], P[3], P[17], P[11], d);
		EdwardsAdditionSubtraction(P[14], P[5], P[19], P[9], d);
		break;
	case 21:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsDouble(P[7], P[14], d);
		EdwardsAdditionSubtraction(P[14], P[1], P[15], P[13], d);
		EdwardsAdditionSubtraction(P[14], P[3], P[17], P[11], d);
		EdwardsAdditionSubtraction(P[14], P[5], P[19], P[9], d);
		EdwardsAddition(P[14], P[7], P[21], d);
		break;
	case 23:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsAddition(P[6], P[3], P[9], d);
		EdwardsDouble(P[9], P[18], d);
		EdwardsNegative(P[7], P[7], d);
		EdwardsAddition(P[18], P[7], P[11], d);
		EdwardsAdditionSubtraction(P[18], P[1], P[19], P[17], d);
		EdwardsAdditionSubtraction(P[18], P[3], P[21], P[15], d);
		EdwardsAdditionSubtraction(P[18], P[5], P[23], P[13], d);
		break;
	case 25:
		EdwardsDouble(P[1], P[2], d);
		EdwardsAddition(P[1], P[2], P[3], d);
		EdwardsAddition(P[3], P[2], P[5], d);
		EdwardsDouble(P[5], P[10], d);
		EdwardsAdditionSubtraction(P[10], P[1], P[11], P[9], d);
		EdwardsAdditionSubtraction(P[10], P[3], P[13], P[7], d);
		EdwardsDouble(P[10], P[20], d);
		EdwardsAdditionSubtraction(P[20], P[1], P[21], P[19], d);
		EdwardsAdditionSubtraction(P[20], P[3], P[23], P[17], d);
		EdwardsAdditionSubtraction(P[20], P[5], P[25], P[15], d);
		break;
	case 27:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsDouble(P[6], P[12], d);
		EdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d);
		EdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d);
		EdwardsDouble(P[11], P[22], d);
		EdwardsAdditionSubtraction(P[22], P[1], P[23], P[21], d);
		EdwardsAdditionSubtraction(P[22], P[3], P[25], P[19], d);
		EdwardsAdditionSubtraction(P[22], P[5], P[27], P[17], d);
		break;
	case 29:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsDouble(P[6], P[12], d);
		EdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d);
		EdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d);
		EdwardsAddition(P[12], P[5], P[17], d);
		EdwardsDouble(P[12], P[24], d);
		EdwardsAdditionSubtraction(P[24], P[1], P[25], P[23], d);
		EdwardsAdditionSubtraction(P[24], P[3], P[27], P[21], d);
		EdwardsAdditionSubtraction(P[24], P[5], P[29], P[19], d);
		break;
	case 31:
		EdwardsTriple(P[1], P[3], d);
		EdwardsDouble(P[3], P[6], d);
		EdwardsAdditionSubtraction(P[6], P[1], P[7], P[5], d);
		EdwardsDouble(P[6], P[12], d);
		EdwardsAdditionSubtraction(P[12], P[1], P[13], P[11], d);
		EdwardsAdditionSubtraction(P[12], P[3], P[15], P[9], d);
		EdwardsDouble(P[12], P[24], d);
		EdwardsAdditionSubtraction(P[24], P[1], P[25], P[23], d);
		EdwardsAdditionSubtraction(P[24], P[3], P[27], P[21], d);
		EdwardsAdditionSubtraction(P[24], P[5], P[29], P[19], d);
		EdwardsAdditionSubtraction(P[24], P[7], P[31], P[17], d);
		break;
	default:
		cout << "Can not deal!" << endl;//不能处理的情况
	}
	
	return 0;
}