

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include "big.h"

using namespace std;
long long cpucycles(void)
{
	return __rdtsc();
}


miracl *mip;

const int infinite = 10000000;
const int bLength = 700;
const int tLength = 400;




typedef struct {
public:
	int sign;//1 represents positive, -1 represents negative
	int  b;
	int  t;
	Big numberN; //The number of DBCs for n
	Big number2N;  //The number of DBCs for n+2^b3^t
	Big number3N; //The number of DBCs for n+2\cdot 2^b3^t
	Big numberNegativeN;//The number of DBCs for n-2^b3^t
	Big numberNegative2N;//The number of DBCs for n-2\cdot 2^b3^t
	Big numberNegative3N;//The number of DBCs for n-3\cdot 2^b3^t
	int length;
	int flag;
	//Big n;
}DBC;


Big pow(int b, int t)
{//return a Big =2^b3^t
	Big result = 1;

	result = result << b;
	for (int i = 0; i < t; i++)
		result = result * 3;
	return result;
}


int CalculateDBCNumber(Big n, DBC(*Item)[tLength], int bMax, int tMax)// 
{
	Item[0][0].numberN = 1;
	Big threeInversion = 1;//2^a3^{b-1}
	Big twoInversion = 1;//2^{a-1}3^b
	Big threeInversionTwo = threeInversion * 2;//2^{a+1}3^{b-1}
	Big module = 1;//2^a3^b
	Big module3 = 3;//2^a3^{b+1}
	Big residue = 0;

	//Big moduleEx = 3;//3*module
	//value Item[0][j]
	//Calculte Item[0][j]
	for (int j = 1; j <= tMax; j++)
	{	
		module = module * 3;		
		//cout << "0," << j << ":" << "2^{b}3^{t+1}=" << module3 << endl;
		residue = n%module;
		//cout << "n_{b,t}=" << residue << endl;
		Big indexTable = 0;
		indexTable = residue * 2 * 3 / module;//n_{b,t}/2^{b-1}3^{t-1}
		int indexSwitch = indexTable % 6;
		//cout << "swithcNumber:" << indexSwitch << endl;

		switch (indexSwitch)
		{
		case 0:
		case 1:
			Item[0][j].numberNegative3N = 0;
			Item[0][j].numberNegative2N = Item[0][j - 1].numberNegative3N;
			Item[0][j].numberNegativeN = Item[0][j - 1].numberNegative3N
				+ Item[0][j - 1].numberNegative2N + Item[0][j - 1].numberN;
			Item[0][j].numberN = Item[0][j - 1].numberNegative3N
				+ Item[0][j - 1].numberNegativeN
				+ Item[0][j - 1].numberN + Item[0][j - 1].number2N;
			Item[0][j].number2N = Item[0][j - 1].numberN + Item[0][j - 1].number3N;
			Item[0][j].number3N = 0;
			break;
		case 2:
		case 3:
			Item[0][j].numberNegative3N = 0;
			Item[0][j].numberNegative2N = Item[0][j - 1].numberNegative2N;
			Item[0][j].numberNegativeN = Item[0][j - 1].numberNegative3N
				+ Item[0][j - 1].numberNegative2N + Item[0][j - 1].numberNegativeN
				+ Item[0][j - 1].number2N;
			Item[0][j].numberN = Item[0][j - 1].numberNegative2N
				+ Item[0][j - 1].numberN
				+ Item[0][j - 1].number2N + Item[0][j - 1].number3N;
			Item[0][j].number2N = Item[0][j - 1].number2N;
			Item[0][j].number3N = 0;
			break;
		case 4:
		case 5:
			Item[0][j].numberNegative3N = 0;
			Item[0][j].numberNegative2N = Item[0][j - 1].numberNegative3N
				+ Item[0][j - 1].numberNegativeN;
			Item[0][j].numberNegativeN = Item[0][j - 1].numberNegative2N
				+ Item[0][j - 1].numberNegativeN + Item[0][j - 1].numberN
				+ Item[0][j - 1].number3N;
			Item[0][j].numberN = Item[0][j - 1].numberNegativeN
				+ Item[0][j - 1].number2N + Item[0][j - 1].number3N;
			Item[0][j].number2N = Item[0][j - 1].number3N;
			Item[0][j].number3N = 0;
			break;
		default:cout << "error" << endl;
		}
	}

	module = 1;
	//value Item[i][0]
	//Calculte Item[i][0]
	for (int i = 1; i <= bMax; i++)
	{		
		module = module * 2;
		//cout << i << ",0:" << "2^{b}3^{t}=" << module << endl;
		residue = n%module;
		//cout << "n_{b,t}=" << residue << endl;

		Big indexTable = 0;
		indexTable = residue * 2 * 3 / module;//n_{b,t}/2^{b-1}3^{t-1}
		int indexSwitch = indexTable % 6;

		//cout << "swithcNumber:" << indexSwitch << endl;

		switch (indexSwitch)
		{
		case 0:
		case 1:
		case 2:
			Item[i][0].numberNegative3N = Item[i - 1][0].numberNegative3N;
			Item[i][0].numberNegative2N = Item[i - 1][0].numberNegative3N
				+ Item[i - 1][0].numberNegativeN;
			Item[i][0].numberNegativeN = Item[i - 1][0].numberNegative3N
				+ Item[i - 1][0].numberNegative2N + Item[i - 1][0].numberNegativeN
				+ Item[i - 1][0].number2N;
			Item[i][0].numberN = Item[i - 1][0].numberNegative3N
				+ Item[i - 1][0].numberNegativeN + Item[i - 1][0].numberN
				+ Item[i - 1][0].number2N;
			Item[i][0].number2N = Item[i - 1][0].numberNegativeN
				+ Item[i - 1][0].number2N + Item[i - 1][0].number3N;
			Item[i][0].number3N = Item[i - 1][0].number2N;
			break;
		case 3:
		case 4:
		case 5:
			Item[i][0].numberNegative3N = Item[i - 1][0].numberNegative2N;
			Item[i][0].numberNegative2N = Item[i - 1][0].numberNegative3N
				+ Item[i - 1][0].numberNegative2N + Item[i - 1][0].numberN;
			Item[i][0].numberNegativeN = Item[i - 1][0].numberNegative2N
				+ Item[i - 1][0].numberNegativeN + Item[i - 1][0].numberN
				+ Item[i - 1][0].number3N;
			Item[i][0].numberN = Item[i - 1][0].numberNegative2N
				+ Item[i - 1][0].numberN + Item[i - 1][0].number2N
				+ Item[i - 1][0].number3N;
			Item[i][0].number2N = Item[i - 1][0].numberN
				+ Item[i - 1][0].number3N;
			Item[i][0].number3N = Item[i - 1][0].number3N;
			break;
		default:cout << "error" << endl;
		}
	}

	for (int j = 1; j <= tMax; j++)
	{
		module = pow(0, j);
		for (int i = 1; i <= bMax; i++)
		{		
			module = module * 2;		
			//cout << i << "," << j << ":" << "2^{b}3^{t}=" << module << endl;
			residue = n%module;
			//cout << "n_{b,t}=" << residue << endl;
			Big indexTable = 0;
			indexTable = residue * 2 * 3 / module;//n_{b,t}/2^{b-1}3^{t-1}
			int indexSwitch = indexTable % 6;
			//cout << "swithcNumber:" << indexSwitch << endl;

			switch (indexSwitch)
			{
			case 0:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative3N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN
					+ Item[i][j - 1].numberNegative3N;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegative2N
					+ Item[i][j - 1].numberN
					- Item[i - 1][j - 1].numberNegative3N;
				Item[i][j].numberN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].numberN + Item[i][j - 1].number2N
					- Item[i - 1][j - 1].numberNegative3N - Item[i - 1][j - 1].numberNegativeN
					- Item[i - 1][j - 1].numberN - Item[i - 1][j - 1].number2N;
				Item[i][j].number2N = Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N + Item[i - 1][j].number3N
					+ Item[i][j - 1].numberN + Item[i][j - 1].number3N;
				Item[i][j].number3N = Item[i - 1][j].number2N;
				break;
			case 1:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative3N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN
					+ Item[i][j - 1].numberNegative3N;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegative2N
					+ Item[i][j - 1].numberN
					- Item[i - 1][j - 1].numberNegative2N;
				Item[i][j].numberN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].numberN + Item[i][j - 1].number2N
					- Item[i - 1][j - 1].numberNegative2N - Item[i - 1][j - 1].numberN
					- Item[i - 1][j - 1].number2N - Item[i - 1][j - 1].number3N;
				Item[i][j].number2N = Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N + Item[i - 1][j].number3N
					+ Item[i][j - 1].numberN + Item[i][j - 1].number3N;
				Item[i][j].number3N = Item[i - 1][j].number2N;
				break;
			case 2:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative3N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN
					+ Item[i][j - 1].numberNegative2N;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegative2N
					+ Item[i][j - 1].numberNegativeN + Item[i][j - 1].number2N
					- Item[i - 1][j - 1].numberNegative3N
					- Item[i - 1][j - 1].numberNegativeN;
				Item[i][j].numberN = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number2N
					+ Item[i][j - 1].numberNegative2N + Item[i][j - 1].numberN
					+ Item[i][j - 1].number2N + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].numberNegativeN
					- Item[i - 1][j - 1].number2N - Item[i - 1][j - 1].number3N;
				Item[i][j].number2N = Item[i - 1][j].numberNegativeN
					+ Item[i - 1][j].number2N + Item[i - 1][j].number3N
					+ Item[i][j - 1].number2N;
				Item[i][j].number3N = Item[i - 1][j].number2N;
				break;
			case 3:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative2N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberN
					+ Item[i][j - 1].numberNegative2N;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegative3N + Item[i][j - 1].numberNegative2N
					+ Item[i][j - 1].numberNegativeN + Item[i][j - 1].number2N
					- Item[i - 1][j - 1].numberNegative3N
					- Item[i - 1][j - 1].numberNegative2N
					- Item[i - 1][j - 1].numberN;
				Item[i][j].numberN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberN + Item[i - 1][j].number2N
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegative2N + Item[i][j - 1].numberN
					+ Item[i][j - 1].number2N + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].numberN - Item[i - 1][j - 1].number3N;
				Item[i][j].number2N = Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].number2N;
				Item[i][j].number3N = Item[i - 1][j].number3N;
				break;
			case 4:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative2N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberN
					+ Item[i][j - 1].numberNegative3N
					+ Item[i][j - 1].numberNegativeN;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegative2N + Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].numberN + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].numberNegative3N
					- Item[i - 1][j - 1].numberNegative2N
					- Item[i - 1][j - 1].numberNegativeN
					- Item[i - 1][j - 1].number2N;
				Item[i][j].numberN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberN + Item[i - 1][j].number2N
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].number2N + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].number2N;
				Item[i][j].number2N = Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].number3N;
				Item[i][j].number3N = Item[i - 1][j].number3N;
				break;
			case 5:
				Item[i][j].numberNegative3N = Item[i - 1][j].numberNegative2N;
				Item[i][j].numberNegative2N = Item[i - 1][j].numberNegative3N
					+ Item[i - 1][j].numberNegative2N + Item[i - 1][j].numberN
					+ Item[i][j - 1].numberNegative3N
					+ Item[i][j - 1].numberNegativeN;
				Item[i][j].numberNegativeN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberNegativeN + Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegative2N + Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].numberN + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].numberNegative2N
					- Item[i - 1][j - 1].numberNegativeN
					- Item[i - 1][j - 1].numberN
					- Item[i - 1][j - 1].number3N;
				Item[i][j].numberN = Item[i - 1][j].numberNegative2N
					+ Item[i - 1][j].numberN + Item[i - 1][j].number2N
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].numberNegativeN
					+ Item[i][j - 1].number2N + Item[i][j - 1].number3N
					- Item[i - 1][j - 1].number3N;
				Item[i][j].number2N = Item[i - 1][j].numberN
					+ Item[i - 1][j].number3N
					+ Item[i][j - 1].number3N;
				Item[i][j].number3N = Item[i - 1][j].number3N;
				break;
			default:cout << "error" << endl;
			}
		}
	}
	return 0;
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


int main()
{
	//srand(time(NULL));
	mip = mirsys(2048, 2);//all numbers are less than 2^{1024}
	mip->IOBASE = 10;//The base of number is decimal


	Big n = "314159265358979323846264338327950288419716939937510582097494459230781640628620";//pi the ratio of circle multiply 10^{}

	n = "31415926535897932384626433832795028841971693993751058209749445923078164062862";

	//Big e="271828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774";
	//n = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024";
	//n = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165";

	//n = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559";
	n = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647";

	//n = "31415926535897932384626433832795028841971693993751058209749445923078164062862";
	//n = "3141592653589793238462643383279502884197169399375105820974944";

	 

	cout << bits(n) << endl;


	DBC  Item[bLength][tLength];

	for (int i = 0; i < bLength; i++)
	{
		for (int j = 0; j < tLength; j++)
		{

			Item[i][j].numberN = 0;//Initial value is NULL
			Item[i][j].numberNegativeN = 0;//Initial value is NULL
			Item[i][j].number2N = 0;//Initial value is NULL
			Item[i][j].numberNegative2N = 0;//Initial value is NULL
			Item[i][j].number3N = 0;//Initial value is NULL
			Item[i][j].numberNegative3N = 0;//Initial value is NULL
		}
	}

	cout << "An example of Algorithm 1: pi*10^{120}=3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647" << endl;

	int bMax=240;
	int tMax=120;

	CalculateDBCNumber(n, Item, bMax, tMax);


	cout << bits(n) << endl;
	Big nCompare = 0;
	nCompare = pow(bMax, tMax);
	if (nCompare>n)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].numberNegative3N
			+ Item[bMax][tMax].numberNegativeN
			+ Item[bMax][tMax].numberN
			+ Item[bMax][tMax].number2N
			<< "." << endl;
	}
	else if (nCompare>n / 2)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].numberNegative2N
			+ Item[bMax][tMax].numberN
			+ Item[bMax][tMax].number2N
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n / 3)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].numberNegativeN
			+ Item[bMax][tMax].number2N
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n / 4)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].numberN
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n / 5)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].number2N << "." << endl;
	}
	else if (nCompare>n / 6)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is " << endl
			<< Item[bMax][tMax].number3N << "." << endl;
	}
	else
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for pi*10^{120} is  0"
			<< "." << endl;
	}
	 



	cout << "please input a positive integer and the power of 2 bMax (less than 700) and the power of 3 tMax(less than 400)" << endl;

	cout << "If one wants to use  bigger bMax and tMax,Please change the value of bLength and tLength in the number.cpp" << endl;
	cout << " Please input a positive integer(less than 700 bit)." << endl;
	cin >> n;
	cout << "bMax<700. Please input" << endl;
	cin >> bMax;
	//bMax = 7;

	while (bMax > 699)
	{
		cout << "bMax<700. Please input again" << endl;
		cin >> bMax;
	}

	cout << "tMax<400. Please input" << endl;
	cin >> tMax;
	//tMax = 0;


	while (tMax > 399)
	{
		cout << "tMax<400. Please input again" << endl;
		cin >> tMax;
	}
	cout << "The bits of your input is " << bits(n) << "." << endl;
	//Big nCompare = 0;
	nCompare = pow(bMax, tMax);
	//if (n > 2 * nCompare)
	//{
	//	cout << "The number of DBCs with a leading term dividing  2^{" << bMax << "}3^{" << tMax << "} for " << n << " is 0." << endl;
	//	return 0;
	//}
	//n=randbits(256);
	double start = cpucycles();

	int indexT = 300;
	//n = 0;
	CalculateDBCNumber(n, Item, bMax, tMax);
	double end = cpucycles();
	//cout << Item[500][indexT].numberN+ Item[500][indexT].numberNegativeN <<endl<<endl;
	/*cout << (end - start) / 1000000 << endl<<endl;


	cout << "The number of DBCs with a leading term dividing  2^{" << bMax << "}3^{" 
		<< tMax << "} for " << n
		<< " is " << endl
		<<"\\begin{equation*}"<<endl
		<<"\\begin{aligned}"<<endl
		<< "|\\mathit{\\operatorname{\\Phi}}(" << bMax << "," << tMax << ","
		<< nCompare*(-3) << ")=&"
		<< Item[bMax][tMax].numberNegative3N << ",&" << endl
		<< "|\\mathit{\\operatorname{\\Phi}}("<< bMax<<","<< tMax<<"," 
		<< nCompare*(-2) <<")=&"
		<< Item[bMax][tMax].numberNegative2N << ",&" << endl
		<< "|\\mathit{\\operatorname{\\Phi}}(" << bMax << "," << tMax << ","
		<< nCompare*(-1) << ")=&"
		<< Item[bMax][tMax].numberNegativeN << ",\\\\" << endl
		<< "|\\mathit{\\operatorname{\\Phi}}(" << bMax << "," << tMax << ","
		<< "0" << ")=&"
		<< Item[bMax][tMax].numberN << ",&" << endl
		<< "|\\mathit{\\operatorname{\\Phi}}(" << bMax << "," << tMax << ","
		<< nCompare << ")=&"
		<< Item[bMax][tMax].number2N << ",&" << endl
		<< "|\\mathit{\\operatorname{\\Phi}}(" << bMax << "," << tMax << ","
		<< nCompare*(2) << ")=&"
		<< Item[bMax][tMax].number3N << ".\\\\" << endl
		<< "\\end{aligned}" << endl
		<< "\\end{equation*}" << endl;*/
		

	if (nCompare>n  )
	{
		cout << "The number of DBCs with a leading factor dividing  2^{" 
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl 
			<< Item[bMax][tMax].numberNegative3N
			+ Item[bMax][tMax].numberNegativeN
			+ Item[bMax][tMax].numberN
			+ Item[bMax][tMax].number2N
			<< "." << endl;		
	}
	else if (  nCompare>n / 2)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{" 
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl 
			<< Item[bMax][tMax].numberNegative2N
			+ Item[bMax][tMax].numberN
			+ Item[bMax][tMax].number2N
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n / 3)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl
			<< Item[bMax][tMax].numberNegativeN
			+ Item[bMax][tMax].number2N
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n / 4)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl
			<< Item[bMax][tMax].numberN		
			+ Item[bMax][tMax].number3N
			<< "." << endl;
	}
	else if (nCompare>n /5)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl
			<< Item[bMax][tMax].number2N<< "." << endl;
	}
	else if (nCompare>n / 6)
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for " << n << " is " << endl
			<< Item[bMax][tMax].number3N << "." << endl;
	}
	else
	{
		cout << "The number of DBCs with a leading factor dividing  2^{"
			<< bMax << "}3^{" << tMax << "} for " << n << " is 0" 
			<<  "." << endl;
	}

	/*for (int i = 0; i < 200; i++)
	{
		cout <<"i:"<<i<<endl			
			<< Item[0][i].numberNegative3N
			+ Item[0][i].numberNegativeN
			+ Item[0][i].numberN
			+ Item[0][i].number2N
			<< ":" << endl
			<< 2*(Item[0][i].numberNegative3N
			+ Item[0][i].numberNegativeN
			+ Item[0][i].numberN
			+ Item[0][i].number2N)
			<< "." << endl;
	}*/

	return 0;
}
