
#include "AdaptiveDBC.h"

long long cpucycles(void)
{
	return __rdtsc();
}

Big pow(int b, int t, int q)
{//return a Big =2^b3^t5^q
	Big result = 1;
	Big pow3 = 3;
	Big pow5 = 5;
	while (t > 0)
	{
		if ((t % 2) == 1)
		{
			result *= pow3;
		}
		pow3 = pow3*pow3;
		t = t / 2;
	}
	while (q > 0)
	{
		if ((q % 2) == 1)
		{
			result *= pow5;
		}
		pow5 = pow5*pow5;
		t = q / 2;
	}
	result = result << b;
	return result;
}

Big pow1(int b, int t)
{//return a Big =2^b3^t
	Big result = 1;
	Big pow3 = 3;
	while (t > 0)
	{
		if ((t % 2) == 1)
		{
			result *= pow3;
		}
		pow3 = pow3*pow3;
		t = t / 2;
	}
	result = result << b;
	return result;
}


//long long power(int b, int t, int q)
//{//return a Big =2^b3^t
//	long long result = 1;
//	result = result << b;
//	for (int i = 0; i < t; i++)
//		result = result * 3;
//	for (int i = 0; i < q; i++)
//		result = result * 5;
//	return result;
//}


long long power(int b, int t, int q)
{//return a Big =2^b3^t5^q
	long long result = 1;
	long long pow3 = 3;
	long long pow5 = 5;
	while (t > 0)
	{
		if ((t % 2) == 1)
		{
			result *= pow3;
		}
		pow3 = pow3*pow3;
		t = t / 2;
	}
	while (q > 0)
	{
		if ((q % 2) == 1)
		{
			result *= pow5;
		}
		pow5 = pow5*pow5;
		t = q / 2;
	}
	result = result << b;
	return result;
}


int positive(int i)
{
	if (i < 0)
	{
		return i + 1;
	}
	return i;
}


Big rand(int n)
{
	Big value = 1;
	for (int i = 1; i < n; i++)
	{
		value *= 2;
		value = value + (rand() % 2);
	}
	return value;
}
int weishu(Big n)
{
	int weishu = 0;
	while (n != 0)
	{
		n = n / 10;
		weishu++;
	}
	return weishu;
}

double value(int bmax, int tmax, int length, double C1, double C2, double C3)
{//return the value of the DBC, which is the value of the value function in our dynamic programming
	return bmax*C2 + tmax*C3 + (length - 1)*C1;
}

int Cmap(int c)
{
	int realC;
	if (c % 2 == 0)
	{
		realC = c / 2;
	}
	else
	{
		realC = -(c + 1) / 2;
	}
	return realC;
}


int CmapInverse(int c)
{//return to the Item[i][j][k][realC]
	int realC;
	if (c >= 0)
	{
		realC = c * 2;
	}
	else
	{
		realC = -2 * c - 1;
	}
	return realC;
}

int intervalValue(Big n, int i, int j, int k)
{
	int interval = -1;
	Big x = pow(i, j, k);
	if (x > n)
	{
		interval = 0;
		return interval;
	}
	else
	{
		Big y = n / x;
		if (y > 1024)
		{
			cout << "This work can not deal! Cmax is large than 512" << endl;
			return 1024;
		}
		interval = y % 1024;
		return interval;
	}
	return -1;
}



int copy(DBC *target, DBC *source)
{
	target->b = source->b;
	target->t = source->t;
	//target->q = source->q;
	target->coefficient = source->coefficient;
	target->additionCost = source->additionCost;
	target->HammingWeight = source->HammingWeight;
	//target->direction = source->direction;
	target->flag = source->flag;
	return 0;
}
int assignment(DBC *target, DBC *source, int b, int t, int coefficient, double additionCost, int flag)
{
	target->b = b;
	target->t = t;
	//target->q = q;
	target->coefficient = coefficient;
	target->additionCost = source->additionCost + additionCost;
	target->HammingWeight = source->HammingWeight + 1;

	target->flag = flag;
	return 0;
}
DBC DynamicProgrammingDBC(Big n, DBC Item[MatrixB][MatrixT][2 * Cmax], 
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles)
{
	const int Infinite = 100000000;
	DBC DBCfinal;
	DBCfinal.b = 0;
	DBCfinal.t = 0;
	DBCfinal.HammingWeight = Infinite;
	DBCfinal.additionCost = Infinite;
	DBCfinal.coefficient = 0;
	DBCfinal.flag = -1;
	long long CPUcycleStart = cpucycles();
	for (int i = 0; i < MatrixB; i++)
	{
		for (int j = 0; j < MatrixT; j++)
		{
			for (int k = 0; k < 2 * Cmax; k++)
			{
				Item[i][j][k].additionCost = Infinite;
				Item[i][j][k].HammingWeight = Infinite;
				Item[i][j][k].b = 0;
				Item[i][j][k].t = 0;
				Item[i][j][k].coefficient = 0;
				Item[i][j][k].flag = -1;
			}
		}
	}


	
	Item[0][0][0].additionCost = 0;
	Item[0][0][0].HammingWeight = 0;
	Big rho = (2 * Cmax + 1)*n;
	const int alpha = 8;
	const int beta = 8;
	int powValue[alpha][beta];
	int finalChainLength = 0;
	const int infiniteCost = Infinite;
	const int infinite = Infinite;
	DBC DBCmiddle;
	 
	double finalValue = Infinite;
	int finalB = 0;
	int finalT = 0;
	int finalSign = 0;
	double log32 = log(2.0) / log(3.0);
	//cout <<"log 3:"<<log32<<endl;

	Big nBound = n * 2 * Cmax;//adaptive	
	double valuemiddle = 0;//value
	int lengthmiddle = 0;//the value of the middle length, used to compare with others			 

	int bBound = (logb2(nBound.getbig()));///int bBound=logb2(nBound.getbig())-1;
	int tBound = 0;



	Big bBoundLarge = 1;
	//int reduceTBound[bLength];
	//int reduceBBound[tLength];

	bBoundLarge = bBoundLarge << bBound;

	if (bBoundLarge > nBound)
	{//bBound
		bBound -= 1;
	}
	Big middlenegative = 1;
	Big residue = 0;
	Big module = 1;
	Big quotient = n;
	Big reducedResidue = 0;
	Big reducedResidueMiddle = 0;
	Big reducedModule = 0;
	Big realModlue = 0;
	Big n6 = 0;

	tBound = floor(log32*(logb2(nBound.getbig())));
	if (pow(0, tBound) > nBound)
	{//bBound
		tBound -= 1;
	}
	int reduceBBound[MatrixB];


	Big BoundSpecialMiddle = n * 2 * Cmax;//adaptive
	int tBoundMiddle = 0;
	for (int j = 0; j <= MatrixT; j++)
	{//reduceBBound[j]=log B/3^j
	 //2^equivalentTBound[j]  3^{j} exactly less than 2n	
		tBoundMiddle = logb2(BoundSpecialMiddle.getbig());//tBoundLarge =pow( bBound,j);														
		bBoundLarge = pow(tBoundMiddle, j);
		if (bBoundLarge > nBound)
		{//bBound
			tBoundMiddle -= 1;
		}
		BoundSpecialMiddle = BoundSpecialMiddle / 3;
		reduceBBound[j] = tBoundMiddle;
	}

	const double log3 = 1.5849625007211561814537389439478;
	int bLength = bits(rho);
	int tLength = bLength / log3;
	if (bLength > bBound)
	{
		bLength = bBound;
	}
	if (tLength > tBound)
	{
		tLength = tBound;
	}

	int alphaSquare = alpha*alpha;
	int betaSquare = beta*beta;
	Big Module1Square = pow(alphaSquare + 1, betaSquare + 1);
	Big n1 = 0;
	Big n1Previous = 0;
	Big n1Prime = 0;
	Big n1PrimeMidlle = 0;
	Big n1PrimePrevious = 0;

	int  Module1 = power(alpha + 1, beta + 1);
	int reducedMiddleResidue = 0;
	int block = power(alpha, beta);
	Big middle = 0;
	int middleTestn1 = 0;
	int blockjj = power(0, beta);
	int blockii = power(alpha, 0);
	Big blockjjSquare = pow(0, betaSquare);
	Big blockiiSquare = pow(alphaSquare, 0);
	Big Moduleii = 1;
	Moduleii = pow(alpha, 0);// 
	Big Modulejj = 1;
	Big newModulejj = 1;

	n6 = 6 * n;

	powValue[0][0] = 1;
	for (int j = 0; j < beta; j++)
	{
		if (j > 0)
		{
			powValue[0][j] = powValue[0][j - 1] * 3;
		}

		for (int i = 1; i < alpha; i++)
		{
			powValue[i][j] = (powValue[i - 1][j]) << 1;
			//cout << powValue[i][j] << endl;
		}
	}
	reducedResidue = n6;
	int ii1Index = 0;
	int jj1Index = 0;
	int i1Index = 0;
	int j1Index = 0;
	int jIndex = 0;
	int iIndex = 0;
	int i = 0;
	int j = 0;
	int xmiddle;
	int ymiddle;
	int reducedXMiddle, reducedYMiddle;
	

	for (int jj1 = 0; jj1 < tBound / betaSquare + 1; jj1++)
	{
		n1PrimePrevious = reducedResidue;
		reducedResidue = reducedResidue / blockjjSquare;
		ii1Index = 0;
		for (int ii1 = 0; ii1 < reduceBBound[jj1*betaSquare] / alphaSquare + 1; ii1++)
		{
			//n1PrimeMidlle = n1PrimePrevious;
			//n1Prime = n1PrimeMidlle%Module1Square;
			n1Prime = n1PrimePrevious%Module1Square;
			//n1PrimePrevious = n1PrimePrevious / blockiiSquare;
			n1PrimePrevious = n1PrimePrevious >> alphaSquare;
			j1Index = jj1Index;
			for (int j1 = 0; j1 < beta; j1++)
			{
				n1Previous = n1Prime;
				n1Prime = n1Prime / blockjj;
				i1Index = ii1Index;
				for (int i1 = 0; i1 < alpha; i1++)
				{
					//n1 = n1Previous;
					reducedMiddleResidue = n1Previous%Module1;
					//reducedMiddleResidue = n1;
					n1Previous = n1Previous / blockii;
					//middle = n6 / pow(ii1*alphaSquare + i1*alpha, jj1*betaSquare + j1*beta);
					//middleTestn1 = middle%Module1;
					//if (reducedMiddleResidue != middleTestn1)
					//{
					//	cout << "6n:" << n6 << endl;
					//	cout << "n1PrimePrevious:" << n1PrimePrevious << " n1Prime:" << n1Prime << " n1Previous:" << n1Previous << " n1: "<<n1 << endl;
					//	cout << "alpha1Square" << blockiiSquare<<endl;
					//	cout << "n: " << n << " ii1: " << ii1 << " jj1: " << jj1 << " i1: " << i1 << " j1: " << j1 << endl;
					//  cout << "n1: " << reducedMiddleResidue<< endl;
					//	cout << "Real n1: " << middleTestn1<< endl;
					//	cout << "wrong of the equation!" << endl;
					//	
					//}
					//iIndex = i1Index;
					//jIndex = j1Index;

					for (int l = 0; l < beta; l++)
					{

						j = j1Index + l;
						//int j = jj*beta + l;
						//int j = jj1*betaSquare+j1*beta + l;
						//jIndex = j1Index + l;
						//if (jIndex != j)
						//{
						//	cout << "Index error!" << jIndex <<" :" <<j << endl;
						//}
						if (j > tBound)
						{
							break;
						}
						for (int k = 0; k < alpha; k++)
						{
							i = i1Index + k;////
											//int i = ii1*alphaSquare + i1*alpha + k;
											//iIndex = i1Index + k;
											//if (iIndex != i)
											//{
											//	cout << "Index error!" << iIndex << " :" << i << endl;
											//}
							if (i > reduceBBound[j])
							{
								break;
							}


							//int nowResidue = reducedMiddleResidue / nowModule;
							int nowResidue = reducedMiddleResidue / powValue[k][l];
							int reducedMiddle = nowResidue % 6;
							//reducedXMiddle = reducedMiddle / 3;//modified on Dec. 30 2018
							//xmiddle = reducedXMiddle % 2;
							xmiddle = reducedMiddle / 3;
							int nowIBit = xmiddle;
							ymiddle = reducedMiddle / 2;
							int nowJBit = ymiddle;
							//reducedYMiddle = reducedMiddle / 2;//modified on Dec. 30 2018
							//ymiddle = reducedYMiddle % 3;
							//Lemma 6 start
							//通过2的方向计算Item[i][j]和ItemNegative[i][j]		 

							if (i > 0)
							{
								if (nowIBit == 0 || nowIBit == 1)
								{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
								 //In this case n_{b,t,q}=n_{b-1,t,q}+2^{b-1}3^{t}5^{q}=n_{b-1,t,q,1}
									for (int c = 0; c < 2 * Cmax; c++)
									{
										int realC = Cmap(c);
										int alphaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b-1,t,q,possibleI} + alphaLemma*2^{b-1}3^{t}5^{q}
														   //alpha==0
										if (alphaLemma == 0)
										{//In this case n_{b,t,q}=n_{b-1,t,q}
											int possibleI = 2 * realC + nowIBit;
											if (possibleI<Cmax && possibleI>-Cmax - 1)
											{//It is legal
											 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
											 //i.e. n_{b,t,q,c}=n_{b-1,t,q,possibleI}

												int sourceI = CmapInverse(possibleI);
												if (Item[i - 1][j][sourceI].additionCost < Infinite)
												{//alpha==0
													if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost)
													{
														copy(&Item[i][j][c], &Item[i - 1][j][sourceI]);
													}
												}
											}
										}
										for (alphaLemma = -Cmax; alphaLemma <= Cmax; alphaLemma = alphaLemma + 2)
										{
											int possibleI = 2 * realC + nowIBit - alphaLemma;
											if (possibleI<Cmax && possibleI>-Cmax - 1)
											{//It is legal
											 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
											 //i.e. n_{b,t,q,c}=n_{b-1,t,q,sourceI}
												int sourceI = CmapInverse(possibleI);
												if (Item[i - 1][j][sourceI].additionCost < Infinite)
												{
													double additionCost = A;
													if (alphaLemma == 1 || alphaLemma == -1)
													{
														additionCost = mA;
													}
													if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost + additionCost)
													{
														assignment(&Item[i][j][c], &Item[i - 1][j][sourceI], i - 1, j, alphaLemma, additionCost, sourceI);
													}
												}
											}
										}

									}
								}
								else
								{
									cout << "direction from binary is wrong!" << endl;
								}
							}


							if (j > 0)
							{
								if (nowJBit > -1 && nowJBit < 3)
								{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
								 //In this case n_{b,t,q}=n_{b,t-1,q}+2^{b}3^{t-1}5^{q}
									for (int c = 0; c < 2 * Cmax; c++)
									{
										int realC = Cmap(c);
										int betaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b,t-1,q,possibleJ} + betaLemma*2^{b}3^{t-1}5^{q}

										if (betaLemma == 0)
										{//In this case n_{b,t,q}=n_{b,t-1,q}
											int possibleJ = 3 * realC + nowJBit;
											if (possibleJ<Cmax && possibleJ>-Cmax - 1)
											{//It is legal
											 //n_{b,t,q}=n_{b,t-1,q,possibleJ}-2*realC*2^{b}3^{t-1}5^{q}
											 //i.e. n_{b,t,q,c}=n_{B,t-1,q,sourceJ}
												int sourceJ = CmapInverse(possibleJ);
												if (Item[i][j - 1][sourceJ].additionCost < Infinite)
												{//alpha==0

													if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost)
													{
														copy(&Item[i][j][c], &Item[i][j - 1][sourceJ]);
													}
												}
											}
										}
										for (betaLemma = -Cmax; betaLemma <= Cmax; betaLemma = betaLemma + 2)
										{
											int possibleJ = 3 * realC + nowJBit - betaLemma;
											if (possibleJ<Cmax && possibleJ>-Cmax - 1)
											{//It is legal
											 //n_{b,t,q}=n_{b,t-1,q,sourceJ}-2*realC*2^{b}3^{t-1}5^{q}
											 //i.e. n_{b,t,q,c}=n_{b,t-1,q,sourceJ}
												int sourceJ = CmapInverse(possibleJ);
												if (Item[i][j - 1][sourceJ].additionCost < Infinite)
												{
													double additionCost = A;
													if (betaLemma == 1 || betaLemma == -1)
													{
														additionCost = mA;
													}
													if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost + additionCost)
													{
														assignment(&Item[i][j][c], &Item[i][j - 1][sourceJ], i, j - 1, betaLemma, additionCost, sourceJ);
													}
												}
											}
										}

									}
								}
								else
								{
									cout << "direction from ternary is wrong!" << endl;
								}
							}
						}//cout << endl;						
					}
					i1Index += alpha;
				}
				j1Index += beta;
			}
			ii1Index += alphaSquare;
		}
		jj1Index += betaSquare;
		//Modulejj = Modulejj*blockjj;
	}
	long long CPUcycleEnd = cpucycles();
	CPUcycles = CPUcycleEnd - CPUcycleStart;
	CPUcycles = CPUcycles / 100;

	int coefficientSign = 1;
	double lowerBound = bits(n / (2 * Cmax)) - 1;
	double upperBound = bits(n * (2 * Cmax));
	//Calculate TBCfinal
	double finalCost = Infinite;
	//FinalCost
	for (int i = 0; i <= bLength; i++)
	{
		for (int j = 0; j <= tLength; j++)
		{
			double determineInterval = i + j*log3;

			if (determineInterval > upperBound)
			{
				break;
			}
			if (determineInterval < lowerBound)
			{
				continue;
			}
			int interval = -1;//determine the interval  pow(i, j, k)>n with interval=0;
							  // n/(interval+1)<pow(i, j, k)le n/interval
			interval = intervalValue(n, i, j, 0);
			int CoefficientAlpha = 0;
			int realC = interval;
			if (realC<Cmax &&realC>-Cmax - 1)
			{
				int cIndex = CmapInverse(realC);
				//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost)
				double firstItemCost = A;
				if (Item[i][j][cIndex].coefficient == 1 || Item[i][j][cIndex].coefficient == -1)
				{
					firstItemCost = mA;
				}
				double finalCostMiddle = Item[i][j][cIndex].additionCost + D*Item[i][j][cIndex].b + T*Item[i][j][cIndex].t
					- firstItemCost;
				if (finalCost > finalCostMiddle)
				{
					finalCost = finalCostMiddle;
					copy(&DBCfinal, &Item[i][j][cIndex]);
				}
			}

			for (int c = -Cmax; c < Cmax + 1; c = c + 2)
			{

				int AddtionCost = A;
				if (c == 1 || c == -1)
				{
					AddtionCost = mA;
				}
				realC = interval - c;
				if (realC<Cmax &&realC>-Cmax - 1)
				{
					int cIndex = CmapInverse(realC);
					//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost+ AddtionCost)
					double finalCostMiddle = Item[i][j][cIndex].additionCost + D*i
						+ T*j;
					if (finalCost > finalCostMiddle)
					{
						finalCost = finalCostMiddle;						
						assignment(&DBCfinal, &Item[i][j][cIndex], i, j, c, AddtionCost, cIndex);
					}
				}
			}
		}
	}
	//cout << "finalCost:" << finalCost << endl;
	//cout << "Hamming weight:" << DBCfinal.HammingWeight << endl;
	//cout << DBCfinal.coefficient << " *2^" << DBCfinal.b << " 3^" << DBCfinal.t << endl;

	DBC NowTBC;

	copy(&NowTBC, &Item[DBCfinal.b][DBCfinal.t][DBCfinal.flag]);
	//out the result chain
	finalChain[DBCfinal.HammingWeight - 1].sign = DBCfinal.coefficient;
	finalChain[DBCfinal.HammingWeight - 1].b = DBCfinal.b;
	finalChain[DBCfinal.HammingWeight - 1].t = DBCfinal.t;

	int presentSign = 1;
	for (int i = DBCfinal.HammingWeight - 1; i > 0; i--)
	{
		finalChain[i - 1].sign = NowTBC.coefficient;
		finalChain[i - 1].b = NowTBC.b;
		finalChain[i - 1].t = NowTBC.t;
		//cout << "Hamming weight:" << NowTBC.HammingWeight << endl;
		//cout << NowTBC.coefficient << " *2^" << NowTBC.b << " 3^" << NowTBC.t << endl;
		int middleB;
		int middleT;
		int middleFlag;
		int middleHammingWeight;
		int middleCoefficient;
		middleB = Item[NowTBC.b][NowTBC.t][NowTBC.flag].b;
		middleT = Item[NowTBC.b][NowTBC.t][NowTBC.flag].t;
		middleFlag = Item[NowTBC.b][NowTBC.t][NowTBC.flag].flag;
		middleHammingWeight = Item[NowTBC.b][NowTBC.t][NowTBC.flag].HammingWeight;
		middleCoefficient = Item[NowTBC.b][NowTBC.t][NowTBC.flag].coefficient;
		NowTBC.b = middleB;
		NowTBC.t = middleT;
		NowTBC.flag = middleFlag;
		NowTBC.HammingWeight = middleHammingWeight;
		NowTBC.coefficient = middleCoefficient;	
	}


	//for (int i = DBCfinal.HammingWeight - 1; i > -1; i--)
	//{
	//	cout << i << ":" << finalChain[i].sign << "*2^" <<
	//		finalChain[i].b << " 3^" <<
	//		finalChain[i].t << endl;
	//}
	return DBCfinal;
}


//DBC DynamicProgrammingDBCBernstein(Big n, DBC Item[MatrixB][MatrixT][2 * Cmax],
//	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles)
//{
//	const int Infinite = 100000000;
//	DBC DBCfinal;
//	DBCfinal.b = 0;
//	DBCfinal.t = 0;
//	DBCfinal.HammingWeight = Infinite;
//	DBCfinal.additionCost = Infinite;
//	DBCfinal.coefficient = 0;
//	DBCfinal.flag = -1;
//	long long CPUcycleStart = cpucycles();
//	for (int i = 0; i < MatrixB; i++)
//	{
//		for (int j = 0; j < MatrixT; j++)
//		{
//			for (int k = 0; k < 2 * Cmax; k++)
//			{
//				Item[i][j][k].additionCost = Infinite;
//				Item[i][j][k].HammingWeight = Infinite;
//				Item[i][j][k].b = 0;
//				Item[i][j][k].t = 0;
//				Item[i][j][k].coefficient = 0;
//				Item[i][j][k].flag = -1;
//			}
//		}
//	}
//
//
//
//	Item[0][0][0].additionCost = 0;
//	Item[0][0][0].HammingWeight = 0;
//	Big rho = (2 * Cmax + 1)*n;
//	const int alpha = 8;
//	const int beta = 8;
//	int powValue[alpha][beta];
//	int finalChainLength = 0;
//	const int infiniteCost = Infinite;
//	const int infinite = Infinite;
//	DBC DBCmiddle;
//
//	double finalValue = Infinite;
//	int finalB = 0;
//	int finalT = 0;
//	int finalSign = 0;
//	double log32 = log(2.0) / log(3.0);
//	//cout <<"log 3:"<<log32<<endl;
//
//	Big nBound = n * 2 * Cmax;//adaptive	
//	double valuemiddle = 0;//value
//	int lengthmiddle = 0;//the value of the middle length, used to compare with others			 
//
//	int bBound = (logb2(nBound.getbig()));///int bBound=logb2(nBound.getbig())-1;
//	int tBound = 0;
//
//
//
//	Big bBoundLarge = 1;
//	//int reduceTBound[bLength];
//	//int reduceBBound[tLength];
//
//	bBoundLarge = bBoundLarge << bBound;
//
//	if (bBoundLarge > nBound)
//	{//bBound
//		bBound -= 1;
//	}
//	Big middlenegative = 1;
//	Big residue = 0;
//	Big module = 1;
//	Big quotient = n;
//	Big reducedResidue = 0;
//	Big reducedResidueMiddle = 0;
//	Big reducedModule = 0;
//	Big realModlue = 0;
//	Big n6 = 0;
//
//	tBound = floor(log32*(logb2(nBound.getbig())));
//	if (pow(0, tBound) > nBound)
//	{//bBound
//		tBound -= 1;
//	}
//	int reduceBBound[MatrixB];
//
//
//	Big BoundSpecialMiddle = n * 2 * Cmax;//adaptive
//	int tBoundMiddle = 0;
//	for (int j = 0; j <= MatrixT; j++)
//	{//reduceBBound[j]=log B/3^j
//	 //2^equivalentTBound[j]  3^{j} exactly less than 2n	
//		tBoundMiddle = logb2(BoundSpecialMiddle.getbig());//tBoundLarge =pow( bBound,j);														
//		bBoundLarge = pow(tBoundMiddle, j);
//		if (bBoundLarge > nBound)
//		{//bBound
//			tBoundMiddle -= 1;
//		}
//		BoundSpecialMiddle = BoundSpecialMiddle / 3;
//		reduceBBound[j] = tBoundMiddle;
//	}
//
//	const double log3 = 1.5849625007211561814537389439478;
//	int bLength = bits(rho);
//	int tLength = bLength / log3;
//	if (bLength > bBound)
//	{
//		bLength = bBound;
//	}
//	if (tLength > tBound)
//	{
//		tLength = tBound;
//	}
//
//	int alphaSquare = alpha*alpha;
//	int betaSquare = beta*beta;
//	Big Module1Square = pow(alphaSquare + 1, betaSquare + 1);
//	Big n1 = 0;
//	Big n1Previous = 0;
//	Big n1Prime = 0;
//	Big n1PrimeMidlle = 0;
//	Big n1PrimePrevious = 0;
//
//	int Module1 = power(alpha + 1, beta + 1);
//	int reducedMiddleResidue = 0;
//	int block = power(alpha, beta);
//	Big middle = 0;
//	int middleTestn1 = 0;
//	int blockjj = power(0, beta);
//	int blockii = power(alpha, 0);
//	Big blockjjSquare = pow(0, betaSquare);
//	Big blockiiSquare = pow(alphaSquare, 0);
//	Big Moduleii = 1;
//	Moduleii = pow(alpha, 0);// 
//	Big Modulejj = 1;
//	Big newModulejj = 1;
//
//	n6 = 6 * n;
//
//	powValue[0][0] = 1;
//	for (int j = 0; j < beta; j++)
//	{
//		if (j > 0)
//		{
//			powValue[0][j] = powValue[0][j - 1] * 3;
//		}
//
//		for (int i = 1; i < alpha; i++)
//		{
//			powValue[i][j] = (powValue[i - 1][j]) << 1;
//			//cout << powValue[i][j] << endl;
//		}
//	}
//	reducedResidue = n6;
//	int ii1Index = 0;
//	int jj1Index = 0;
//	int i1Index = 0;
//	int j1Index = 0;
//	int jIndex = 0;
//	int iIndex = 0;
//	int i = 0;
//	int j = 0;
//	int xmiddle;
//	int ymiddle;
//	int reducedXMiddle, reducedYMiddle;
//
//
//	for (int jj1 = 0; jj1 < tBound / betaSquare + 1; jj1++)
//	{
//		n1PrimePrevious = reducedResidue;
//		reducedResidue = reducedResidue / blockjjSquare;
//		ii1Index = 0;
//		for (int ii1 = 0; ii1 < reduceBBound[jj1*betaSquare] / alphaSquare + 1; ii1++)
//		{
//			//n1PrimeMidlle = n1PrimePrevious;
//			//n1Prime = n1PrimeMidlle%Module1Square;
//			n1Prime = n1PrimePrevious%Module1Square;
//			//n1PrimePrevious = n1PrimePrevious / blockiiSquare;
//			n1PrimePrevious = n1PrimePrevious >> alphaSquare;
//			j1Index = jj1Index;
//
//			
//			for (int j1 = 0; j1 < beta; j1++)
//			{
//				n1Previous = n1Prime;
//				n1Prime = n1Prime / blockjj;
//				i1Index = ii1Index;
//				for (int i1 = 0; i1 < alpha; i1++)
//				{
//					//n1 = n1Previous;
//					reducedMiddleResidue = n1Previous%Module1;
//					//reducedMiddleResidue = n1;
//					n1Previous = n1Previous / blockii;
//					//middle = n6 / pow(ii1*alphaSquare + i1*alpha, jj1*betaSquare + j1*beta);
//					//middleTestn1 = middle%Module1;
//					//if (reducedMiddleResidue != middleTestn1)
//					//{
//					//	cout << "6n:" << n6 << endl;
//					//	cout << "n1PrimePrevious:" << n1PrimePrevious << " n1Prime:" << n1Prime << " n1Previous:" << n1Previous << " n1: "<<n1 << endl;
//					//	cout << "alpha1Square" << blockiiSquare<<endl;
//					//	cout << "n: " << n << " ii1: " << ii1 << " jj1: " << jj1 << " i1: " << i1 << " j1: " << j1 << endl;
//					//  cout << "n1: " << reducedMiddleResidue<< endl;
//					//	cout << "Real n1: " << middleTestn1<< endl;
//					//	cout << "wrong of the equation!" << endl;
//					//	
//					//}
//					//iIndex = i1Index;
//					//jIndex = j1Index;
//
//					for (int l = 0; l < beta; l++)
//					{
//
//						j = j1Index + l;
//						//int j = jj*beta + l;
//						//int j = jj1*betaSquare+j1*beta + l;
//						//jIndex = j1Index + l;
//						//if (jIndex != j)
//						//{
//						//	cout << "Index error!" << jIndex <<" :" <<j << endl;
//						//}
//						Big nowResidue = n6 / pow(i, j);
//						int reducedMiddle1 = nowResidue % 6;
//						if (j > tBound)
//						{
//							break;
//						}
//						for (int k = 0; k < alpha; k++)
//						{
//							i = i1Index + k;////
//											//int i = ii1*alphaSquare + i1*alpha + k;
//											//iIndex = i1Index + k;
//											//if (iIndex != i)
//											//{
//											//	cout << "Index error!" << iIndex << " :" << i << endl;
//											//}
//							if (i > reduceBBound[j])
//							{
//								break;
//							}
//
//
//							//int nowResidue = reducedMiddleResidue / nowModule;
//							int nowResidue = reducedMiddleResidue / powValue[k][l];
//							int reducedMiddle = nowResidue % 6;
//							//reducedXMiddle = reducedMiddle / 3;//modified on Dec. 30 2018
//							//xmiddle = reducedXMiddle % 2;
//							xmiddle = reducedMiddle / 3;
//							int nowIBit = xmiddle;
//							ymiddle = reducedMiddle / 2;
//							int nowJBit = ymiddle;
//							//reducedYMiddle = reducedMiddle / 2;//modified on Dec. 30 2018
//							//ymiddle = reducedYMiddle % 3;
//							//Lemma 6 start
//							//通过2的方向计算Item[i][j]和ItemNegative[i][j]		 
//
//							if (i > 0)
//							{
//								if (nowIBit == 0 || nowIBit == 1)
//								{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
//								 //In this case n_{b,t,q}=n_{b-1,t,q}+2^{b-1}3^{t}5^{q}=n_{b-1,t,q,1}
//									for (int c = 0; c < 2 * Cmax; c++)
//									{
//										int realC = Cmap(c);
//										int alphaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b-1,t,q,possibleI} + alphaLemma*2^{b-1}3^{t}5^{q}
//														   //alpha==0
//										if (alphaLemma == 0)
//										{//In this case n_{b,t,q}=n_{b-1,t,q}
//											int possibleI = 2 * realC + nowIBit;
//											if (possibleI<Cmax && possibleI>-Cmax - 1)
//											{//It is legal
//											 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
//											 //i.e. n_{b,t,q,c}=n_{b-1,t,q,possibleI}
//
//												int sourceI = CmapInverse(possibleI);
//												if (Item[i - 1][j][sourceI].additionCost < Infinite)
//												{//alpha==0
//													if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost)
//													{
//														copy(&Item[i][j][c], &Item[i - 1][j][sourceI]);
//													}
//												}
//											}
//										}
//										for (alphaLemma = -Cmax; alphaLemma <= Cmax; alphaLemma = alphaLemma + 2)
//										{
//											int possibleI = 2 * realC + nowIBit - alphaLemma;
//											if (possibleI<Cmax && possibleI>-Cmax - 1)
//											{//It is legal
//											 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
//											 //i.e. n_{b,t,q,c}=n_{b-1,t,q,sourceI}
//												int sourceI = CmapInverse(possibleI);
//												if (Item[i - 1][j][sourceI].additionCost < Infinite)
//												{
//													double additionCost = A;
//													if (alphaLemma == 1 || alphaLemma == -1)
//													{
//														additionCost = mA;
//													}
//													if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost + additionCost)
//													{
//														assignment(&Item[i][j][c], &Item[i - 1][j][sourceI], i - 1, j, alphaLemma, additionCost, sourceI);
//													}
//												}
//											}
//										}
//
//									}
//								}
//								else
//								{
//									cout << "direction from binary is wrong!" << endl;
//								}
//							}
//
//
//							if (j > 0)
//							{
//								if (nowJBit > -1 && nowJBit < 3)
//								{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
//								 //In this case n_{b,t,q}=n_{b,t-1,q}+2^{b}3^{t-1}5^{q}
//									for (int c = 0; c < 2 * Cmax; c++)
//									{
//										int realC = Cmap(c);
//										int betaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b,t-1,q,possibleJ} + betaLemma*2^{b}3^{t-1}5^{q}
//
//										if (betaLemma == 0)
//										{//In this case n_{b,t,q}=n_{b,t-1,q}
//											int possibleJ = 3 * realC + nowJBit;
//											if (possibleJ<Cmax && possibleJ>-Cmax - 1)
//											{//It is legal
//											 //n_{b,t,q}=n_{b,t-1,q,possibleJ}-2*realC*2^{b}3^{t-1}5^{q}
//											 //i.e. n_{b,t,q,c}=n_{B,t-1,q,sourceJ}
//												int sourceJ = CmapInverse(possibleJ);
//												if (Item[i][j - 1][sourceJ].additionCost < Infinite)
//												{//alpha==0
//
//													if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost)
//													{
//														copy(&Item[i][j][c], &Item[i][j - 1][sourceJ]);
//													}
//												}
//											}
//										}
//										for (betaLemma = -Cmax; betaLemma <= Cmax; betaLemma = betaLemma + 2)
//										{
//											int possibleJ = 3 * realC + nowJBit - betaLemma;
//											if (possibleJ<Cmax && possibleJ>-Cmax - 1)
//											{//It is legal
//											 //n_{b,t,q}=n_{b,t-1,q,sourceJ}-2*realC*2^{b}3^{t-1}5^{q}
//											 //i.e. n_{b,t,q,c}=n_{b,t-1,q,sourceJ}
//												int sourceJ = CmapInverse(possibleJ);
//												if (Item[i][j - 1][sourceJ].additionCost < Infinite)
//												{
//													double additionCost = A;
//													if (betaLemma == 1 || betaLemma == -1)
//													{
//														additionCost = mA;
//													}
//													if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost + additionCost)
//													{
//														assignment(&Item[i][j][c], &Item[i][j - 1][sourceJ], i, j - 1, betaLemma, additionCost, sourceJ);
//													}
//												}
//											}
//										}
//
//									}
//								}
//								else
//								{
//									cout << "direction from ternary is wrong!" << endl;
//								}
//							}
//						}//cout << endl;						
//					}
//					i1Index += alpha;
//				}
//				j1Index += beta;
//			}
//			ii1Index += alphaSquare;
//		}
//		jj1Index += betaSquare;
//		//Modulejj = Modulejj*blockjj;
//	}
//	long long CPUcycleEnd = cpucycles();
//	CPUcycles = CPUcycleEnd - CPUcycleStart;
//	CPUcycles = CPUcycles / 10;
//
//	int coefficientSign = 1;
//	double lowerBound = bits(n / (2 * Cmax)) - 1;
//	double upperBound = bits(n * (2 * Cmax));
//	//Calculate TBCfinal
//	double finalCost = Infinite;
//	//FinalCost
//	for (int i = 0; i <= bLength; i++)
//	{
//		for (int j = 0; j <= tLength; j++)
//		{
//			double determineInterval = i + j*log3;
//
//			if (determineInterval > upperBound)
//			{
//				break;
//			}
//			if (determineInterval < lowerBound)
//			{
//				continue;
//			}
//			int interval = -1;//determine the interval  pow(i, j, k)>n with interval=0;
//							  // n/(interval+1)<pow(i, j, k)le n/interval
//			interval = intervalValue(n, i, j, 0);
//			int CoefficientAlpha = 0;
//			int realC = interval;
//			if (realC<Cmax &&realC>-Cmax - 1)
//			{
//				int cIndex = CmapInverse(realC);
//				//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost)
//				double firstItemCost = A;
//				if (Item[i][j][cIndex].coefficient == 1 || Item[i][j][cIndex].coefficient == -1)
//				{
//					firstItemCost = mA;
//				}
//				double finalCostMiddle = Item[i][j][cIndex].additionCost + D*Item[i][j][cIndex].b + T*Item[i][j][cIndex].t
//					- firstItemCost;
//				if (finalCost > finalCostMiddle)
//				{
//					finalCost = finalCostMiddle;
//					copy(&DBCfinal, &Item[i][j][cIndex]);
//				}
//			}
//
//			for (int c = -Cmax; c < Cmax + 1; c = c + 2)
//			{
//
//				int AddtionCost = A;
//				if (c == 1 || c == -1)
//				{
//					AddtionCost = mA;
//				}
//				realC = interval - c;
//				if (realC<Cmax &&realC>-Cmax - 1)
//				{
//					int cIndex = CmapInverse(realC);
//					//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost+ AddtionCost)
//					double finalCostMiddle = Item[i][j][cIndex].additionCost + D*i
//						+ T*j;
//					if (finalCost > finalCostMiddle)
//					{
//						finalCost = finalCostMiddle;
//						assignment(&DBCfinal, &Item[i][j][cIndex], i, j, c, AddtionCost, cIndex);
//					}
//				}
//			}
//		}
//	}
//	//cout << "finalCost:" << finalCost << endl;
//	//cout << "Hamming weight:" << DBCfinal.HammingWeight << endl;
//	//cout << DBCfinal.coefficient << " *2^" << DBCfinal.b << " 3^" << DBCfinal.t << endl;
//
//	DBC NowTBC;
//
//	copy(&NowTBC, &Item[DBCfinal.b][DBCfinal.t][DBCfinal.flag]);
//	//out the result chain
//	finalChain[DBCfinal.HammingWeight - 1].sign = DBCfinal.coefficient;
//	finalChain[DBCfinal.HammingWeight - 1].b = DBCfinal.b;
//	finalChain[DBCfinal.HammingWeight - 1].t = DBCfinal.t;
//
//	int presentSign = 1;
//	for (int i = DBCfinal.HammingWeight - 1; i > 0; i--)
//	{
//		finalChain[i - 1].sign = NowTBC.coefficient;
//		finalChain[i - 1].b = NowTBC.b;
//		finalChain[i - 1].t = NowTBC.t;
//		//cout << "Hamming weight:" << NowTBC.HammingWeight << endl;
//		//cout << NowTBC.coefficient << " *2^" << NowTBC.b << " 3^" << NowTBC.t << endl;
//		int middleB;
//		int middleT;
//		int middleFlag;
//		int middleHammingWeight;
//		int middleCoefficient;
//		middleB = Item[NowTBC.b][NowTBC.t][NowTBC.flag].b;
//		middleT = Item[NowTBC.b][NowTBC.t][NowTBC.flag].t;
//		middleFlag = Item[NowTBC.b][NowTBC.t][NowTBC.flag].flag;
//		middleHammingWeight = Item[NowTBC.b][NowTBC.t][NowTBC.flag].HammingWeight;
//		middleCoefficient = Item[NowTBC.b][NowTBC.t][NowTBC.flag].coefficient;
//		NowTBC.b = middleB;
//		NowTBC.t = middleT;
//		NowTBC.flag = middleFlag;
//		NowTBC.HammingWeight = middleHammingWeight;
//		NowTBC.coefficient = middleCoefficient;
//	}
//
//
//	//for (int i = DBCfinal.HammingWeight - 1; i > -1; i--)
//	//{
//	//	cout << i << ":" << finalChain[i].sign << "*2^" <<
//	//		finalChain[i].b << " 3^" <<
//	//		finalChain[i].t << endl;
//	//}
//	return DBCfinal;
//}


DBC DynamicProgrammingDBCBernstein(Big n, DBC Item[MatrixB][MatrixT][2 * Cmax], 
	DBCChainItem finalChain[MatrixB], double mA, double A, double D, double T, long long &CPUcycles)
{
	const int Infinite = 100000000;
	DBC DBCfinal;
	DBCfinal.b = 0;
	DBCfinal.t = 0;
	DBCfinal.HammingWeight = Infinite;
	DBCfinal.additionCost = Infinite;
	DBCfinal.coefficient = 0;
	DBCfinal.flag = -1;
	for (int i = 0; i < MatrixB; i++)
	{
		for (int j = 0; j < MatrixT; j++)
		{
			for (int k = 0; k < 2 * Cmax; k++)
			{
				Item[i][j][k].additionCost = Infinite;
				Item[i][j][k].HammingWeight = Infinite;
				Item[i][j][k].b = 0;
				Item[i][j][k].t = 0;
				Item[i][j][k].coefficient = 0;
				Item[i][j][k].flag = -1;
			}
		}
	}
	Item[0][0][0].additionCost = 0;
	Item[0][0][0].HammingWeight = 0;
	Big rho = (2 * Cmax + 1)*n;

	int finalChainLength = 0;
	const int infiniteCost = Infinite;
	const int infinite = Infinite;
	DBC DBCmiddle;
	
	double finalValue = Infinite;

	double log32 = log(2.0) / log(3.0);
	//cout <<"log 3:"<<log32<<endl;

	Big nBound = n * 2 * Cmax;//adaptive	
	

	int bBound = (logb2(nBound.getbig()));///int bBound=logb2(nBound.getbig())-1;
	int tBound = 0;

	Big bBoundLarge = 1;

	bBoundLarge = bBoundLarge << bBound;

	if (bBoundLarge > nBound)
	{//bBound
		bBound -= 1;
	}

	Big reducedResidue = 0; 
	Big n6 = 0;
	tBound = floor(log32*(logb2(nBound.getbig())));
	if (pow(0, tBound) > nBound)
	{//bBound
		tBound -= 1;
	}

	const double log3 = 1.5849625007211561814537389439478;
	int bLength = bits(rho);
	int tLength = bLength / log3;
	if (bLength > bBound)
	{
		bLength = bBound;
	}
	if (tLength > tBound)
	{
		tLength = tBound;
	}
 
	n6 = 6 * n; 
	reducedResidue = n6; 
	int xmiddle;
	int ymiddle;
	int reducedXMiddle, reducedYMiddle;

	long long CPUcycleStart = cpucycles();
	for (int j = 0; j <= tBound; j++)
	{
		for (int i = 0; i <= bBound; i++)
		{
			Big nowResidue = n6 / pow(i, j);
			int reducedMiddle = nowResidue % 6;
			//reducedXMiddle = reducedMiddle / 3;//modified on Dec. 30 2018
			//xmiddle = reducedXMiddle % 2;
			xmiddle = reducedMiddle / 3;
			int nowIBit = xmiddle;
			ymiddle = reducedMiddle / 2;
			int nowJBit = ymiddle;
			//reducedYMiddle = reducedMiddle / 2;//modified on Dec. 30 2018
			//ymiddle = reducedYMiddle % 3;
			//Lemma 6 start
			//通过2的方向计算Item[i][j]和ItemNegative[i][j]		 

			if (i > 0)
			{
				if (nowIBit == 0 || nowIBit == 1)
				{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
				 //In this case n_{b,t,q}=n_{b-1,t,q}+2^{b-1}3^{t}5^{q}=n_{b-1,t,q,1}
					for (int c = 0; c < 2 * Cmax; c++)
					{
						int realC = Cmap(c);
						int alphaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b-1,t,q,possibleI} + alphaLemma*2^{b-1}3^{t}5^{q}
										   //alpha==0
						if (alphaLemma == 0)
						{//In this case n_{b,t,q}=n_{b-1,t,q}
							int possibleI = 2 * realC + nowIBit;
							if (possibleI<Cmax && possibleI>-Cmax - 1)
							{//It is legal
							 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
							 //i.e. n_{b,t,q,c}=n_{b-1,t,q,possibleI}

								int sourceI = CmapInverse(possibleI);
								if (Item[i - 1][j][sourceI].additionCost < Infinite)
								{//alpha==0
									if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost)
									{
										copy(&Item[i][j][c], &Item[i - 1][j][sourceI]);
									}
								}
							}
						}
						for (alphaLemma = -Cmax; alphaLemma <= Cmax; alphaLemma = alphaLemma + 2)
						{
							int possibleI = 2 * realC + nowIBit - alphaLemma;
							if (possibleI<Cmax && possibleI>-Cmax - 1)
							{//It is legal
							 //n_{b,t,q}=n_{b-1,t,q,sourceI}-2*realC*2^{b-1}3^{t}5^{q}
							 //i.e. n_{b,t,q,c}=n_{b-1,t,q,sourceI}
								int sourceI = CmapInverse(possibleI);
								if (Item[i - 1][j][sourceI].additionCost < Infinite)
								{
									double additionCost = A;
									if (alphaLemma == 1 || alphaLemma == -1)
									{
										additionCost = mA;
									}
									if (Item[i][j][c].additionCost > Item[i - 1][j][sourceI].additionCost + additionCost)
									{
										assignment(&Item[i][j][c], &Item[i - 1][j][sourceI], i - 1, j, alphaLemma, additionCost, sourceI);
									}
								}
							}
						}

					}
				}
				else
				{
					cout << "direction from binary is wrong!" << endl;
				}
			}


			if (j > 0)
			{
				if (nowJBit > -1 && nowJBit < 3)
				{//The last [2*i] positive Item[i].  [2*i-1] negative  Item[-i]
				 //In this case n_{b,t,q}=n_{b,t-1,q}+2^{b}3^{t-1}5^{q}
					for (int c = 0; c < 2 * Cmax; c++)
					{
						int realC = Cmap(c);
						int betaLemma = 0;//The coefficient n_{b,t,q,c}=n_{b,t-1,q,possibleJ} + betaLemma*2^{b}3^{t-1}5^{q}

						if (betaLemma == 0)
						{//In this case n_{b,t,q}=n_{b,t-1,q}
							int possibleJ = 3 * realC + nowJBit;
							if (possibleJ<Cmax && possibleJ>-Cmax - 1)
							{//It is legal
							 //n_{b,t,q}=n_{b,t-1,q,possibleJ}-2*realC*2^{b}3^{t-1}5^{q}
							 //i.e. n_{b,t,q,c}=n_{B,t-1,q,sourceJ}
								int sourceJ = CmapInverse(possibleJ);
								if (Item[i][j - 1][sourceJ].additionCost < Infinite)
								{//alpha==0

									if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost)
									{
										copy(&Item[i][j][c], &Item[i][j - 1][sourceJ]);
									}
								}
							}
						}
						for (betaLemma = -Cmax; betaLemma <= Cmax; betaLemma = betaLemma + 2)
						{
							int possibleJ = 3 * realC + nowJBit - betaLemma;
							if (possibleJ<Cmax && possibleJ>-Cmax - 1)
							{//It is legal
							 //n_{b,t,q}=n_{b,t-1,q,sourceJ}-2*realC*2^{b}3^{t-1}5^{q}
							 //i.e. n_{b,t,q,c}=n_{b,t-1,q,sourceJ}
								int sourceJ = CmapInverse(possibleJ);
								if (Item[i][j - 1][sourceJ].additionCost < Infinite)
								{
									double additionCost = A;
									if (betaLemma == 1 || betaLemma == -1)
									{
										additionCost = mA;
									}
									if (Item[i][j][c].additionCost > Item[i][j - 1][sourceJ].additionCost + additionCost)
									{
										assignment(&Item[i][j][c], &Item[i][j - 1][sourceJ], i, j - 1, betaLemma, additionCost, sourceJ);
									}
								}
							}
						}

					}
				}
				else
				{
					cout << "direction from ternary is wrong!" << endl;
				}
			}
		}
	}

	long long CPUcycleEnd = cpucycles();
	CPUcycles = CPUcycleEnd - CPUcycleStart;

	CPUcycles = CPUcycles / 100;

	int coefficientSign = 1;
	double lowerBound = bits(n / (2 * Cmax)) - 1;
	double upperBound = bits(n * (2 * Cmax));
	//Calculate TBCfinal
	double finalCost = Infinite;
	//FinalCost
	for (int i = 0; i <= bLength; i++)
	{
		for (int j = 0; j <= tLength; j++)
		{
			double determineInterval = i + j*log3;

			if (determineInterval > upperBound)
			{
				break;
			}
			if (determineInterval < lowerBound)
			{
				continue;
			}
			int interval = -1;//determine the interval  pow(i, j, k)>n with interval=0;
							  // n/(interval+1)<pow(i, j, k)le n/interval
			interval = intervalValue(n, i, j, 0);
			int CoefficientAlpha = 0;
			int realC = interval;
			if (realC<Cmax &&realC>-Cmax - 1)
			{
				int cIndex = CmapInverse(realC);
				//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost)
				double firstItemCost = A;
				if (Item[i][j][cIndex].coefficient == 1 || Item[i][j][cIndex].coefficient == -1)
				{
					firstItemCost = mA;
				}
				double finalCostMiddle = Item[i][j][cIndex].additionCost + D*Item[i][j][cIndex].b + T*Item[i][j][cIndex].t
					- firstItemCost;
				if (finalCost > finalCostMiddle)
				{
					finalCost = finalCostMiddle;
					copy(&DBCfinal, &Item[i][j][cIndex]);
				}
			}

			for (int c = -Cmax; c < Cmax + 1; c = c + 2)
			{
				int AddtionCost = A;
				if (c == 1 || c == -1)
				{
					AddtionCost = mA;
				}
				realC = interval - c;
				if (realC<Cmax &&realC>-Cmax - 1)
				{
					int cIndex = CmapInverse(realC);
					//if (TBCfinal.additionCost>Item[i][j][k][cIndex].additionCost+ AddtionCost)
					double finalCostMiddle = Item[i][j][cIndex].additionCost + D*i
						+ T*j;
					if (finalCost > finalCostMiddle)
					{
						finalCost = finalCostMiddle; 
						assignment(&DBCfinal, &Item[i][j][cIndex], i, j, c, AddtionCost, cIndex);
					}
				}
			}
		}
	}
	//cout << "finalCost:" << finalCost << endl;
	//cout << "Hamming weight:" << DBCfinal.HammingWeight << endl;
	//cout << DBCfinal.coefficient << " *2^" << DBCfinal.b << " 3^" << DBCfinal.t << endl;

	DBC NowTBC;
	copy(&NowTBC, &Item[DBCfinal.b][DBCfinal.t][DBCfinal.flag]);	 


	finalChain[DBCfinal.HammingWeight - 1].sign = DBCfinal.coefficient;
	finalChain[DBCfinal.HammingWeight - 1].b = DBCfinal.b;
	finalChain[DBCfinal.HammingWeight - 1].t = DBCfinal.t;

	int presentSign = 1;
	for (int i = DBCfinal.HammingWeight - 1; i > 0; i--)
	{
		finalChain[i - 1].sign = NowTBC.coefficient;
		finalChain[i - 1].b = NowTBC.b;
		finalChain[i - 1].t = NowTBC.t;
		//cout << "Hamming weight:" << NowTBC.HammingWeight << endl;
		//cout << NowTBC.coefficient << " *2^" << NowTBC.b << " 3^" << NowTBC.t << endl;
		int middleB;
		int middleT;
		int middleFlag;
		int middleHammingWeight;
		int middleCoefficient;
		middleB = Item[NowTBC.b][NowTBC.t][NowTBC.flag].b;
		middleT = Item[NowTBC.b][NowTBC.t][NowTBC.flag].t;
		middleFlag = Item[NowTBC.b][NowTBC.t][NowTBC.flag].flag;
		middleHammingWeight = Item[NowTBC.b][NowTBC.t][NowTBC.flag].HammingWeight;
		middleCoefficient = Item[NowTBC.b][NowTBC.t][NowTBC.flag].coefficient;
		NowTBC.b = middleB;
		NowTBC.t = middleT;
		NowTBC.flag = middleFlag;
		NowTBC.HammingWeight = middleHammingWeight;
		NowTBC.coefficient = middleCoefficient;		
	}

	//for (int i = DBCfinal.HammingWeight - 1; i > -1; i--)
	//{
	//	cout <<i <<":" <<finalChain[i].sign << "*2^"<<
	//	finalChain[i].b << " 3^"<<
	//	finalChain[i].t <<endl;
	//}
	return DBCfinal;
}