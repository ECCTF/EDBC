#include "scalarMultiplication.h"

int scalarMultiplicationJ(affinePoint P, affinePoint &nP, int *item, int number, int Atomicity)
{
	projectPoint P2;
	P2.X=mirvar(0); P2.Y=mirvar(0);P2.Z=mirvar(1);
	//projectPoint P3;
	//P3.X=mirvar(0); P3.Y=mirvar(0);P3.Z=mirvar(1);
	nres(P2.Z,P2.Z);
	//Big X(0),Y(0),Z(0);
	copy(P.x,P2.X); //P2.X=P.x; 
	copy(P.y,P2.Y); //P2.Y=P.y;
	//P2.Z=1;
	affinePoint negativeP;
	negativeP.x=mirvar(0);
	negativeP.y=mirvar(0);
	copy(P.x,negativeP.x);
	nres_negate(P.y, negativeP.y);
	
	switch (Atomicity)
	{
	case 0:
		for (int i = number - 2; i > -1; i--)
		{
			JacobiDouble(P2, P2);

			if (item[i] == 1)
			{
				JacobiMixedAddition(P2, P, P2);
			}
			else if (item[i] == -1)
			{
				JacobiMixedAddition(P2, negativeP, P2);
			}
		}
		break;
	case 1:
		for (int i = number - 2; i > -1; i--)
		{
			JacobiDoubleAtomicity(P2, P2);

			if (item[i] == 1)
			{
				JacobiMixedAdditionAtomicity(P2, P, P2);
			}
			else if (item[i] == -1)
			{
				JacobiMixedAdditionAtomicity(P2, negativeP, P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}


 

	big mid=mirvar(0);
	big one=mirvar(1);

	//Big x=one; 
	
	nres_modmult(P2.Z,P2.Z,mid);
	nres_modmult(mid,P2.Z,mid);
	//redc(mid, mid);
	nres(one,one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if(compare(mid,one)!=0)
	{		
		//nres(mid,mid);
		nres_moddiv(one,mid,mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(mid,P2.Y,nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid,P2.X,mid);
	nres_modmult(mid,P2.Z,nP.x);



	return 0;
}


int scalarMultiplicationJ(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, int Atomicity)
{

	projectPoint P2;
	P2.X=mirvar(0); P2.Y=mirvar(0);P2.Z=mirvar(1);
	
	nres(P2.Z,P2.Z);
	
	copy(P.x,P2.X); //P2.X=P.x; 
	copy(P.y,P2.Y); //P2.Y=P.y;
	//projectPoint P2;
	
	//P2.X=P.x;
	//P2.Y=P.y;
	//P2.Z=1;


	affinePoint negativeP;
	negativeP.x=mirvar(0);
	negativeP.y=mirvar(0);
	copy(P.x,negativeP.x);
	nres_negate(P.y, negativeP.y);
	////affinePoint negativeP;
	////negativeP.x=P.x;
	////negativeP.y=(-P.y)%p;
	switch (Atomicity)
	{
	case 0:
		if (number == 0)
		{
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				JacobiDouble(P2, P2);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				JacobiTriple(P2, P2);
			}
		}
		else
		{
			for (int j = 0; j < item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				JacobiDouble(P2, P2);
			}
			for (int j = 0; j < item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				JacobiTriple(P2, P2);
			}

			for (int i = number - 1; i > 0; i--)
			{
				if (item[i].sign == 1)
				{
					JacobiMixedAddition(P2, P, P2);
				}
				else if (item[i].sign == -1)
				{
					JacobiMixedAddition(P2, negativeP, P2);
				}
				for (int j = 0; j < item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					JacobiDouble(P2, P2);
				}
				for (int j = 0; j < item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					JacobiTriple(P2, P2);
				}
			}
			if (item[0].sign == 1)
			{
				JacobiMixedAddition(P2, P, P2);
			}
			else if (item[0].sign == -1)
			{
				JacobiMixedAddition(P2, negativeP, P2);
			}
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				JacobiDouble(P2, P2);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				JacobiTriple(P2, P2);
			}
		}
		break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}
		}
		else
		{
			for (int j = 0; j < item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j < item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}

			for (int i = number - 1; i > 0; i--)
			{
				if (item[i].sign == 1)
				{
					JacobiMixedAdditionAtomicity(P2, P, P2);
				}
				else if (item[i].sign == -1)
				{
					JacobiMixedAdditionAtomicity(P2, negativeP, P2);
				}
				for (int j = 0; j < item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					JacobiDoubleAtomicity(P2, P2);
				}
				for (int j = 0; j < item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					JacobiTripleAtomicity(P2, P2);
				}
			}
			if (item[0].sign == 1)
			{
				JacobiMixedAdditionAtomicity(P2, P, P2);
			}
			else if (item[0].sign == -1)
			{
				JacobiMixedAdditionAtomicity(P2, negativeP, P2);
			}
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid=mirvar(0);
	big one=mirvar(1);

	//Big x=one; 
	
	nres_modmult(P2.Z,P2.Z,mid);
	nres_modmult(mid,P2.Z,mid);
	//redc(mid, mid);
	nres(one,one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if(compare(mid,one)!=0)
	{		
		//nres(mid,mid);
		nres_moddiv(one,mid,mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(mid,P2.Y,nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid,P2.X,mid);
	nres_modmult(mid,P2.Z,nP.x);
	//Big mid(0);
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p,p);
	////	cout <<":" <<P2.X <<":"<<P2.Y <<":"<<P2.Z<<endl;
	//nP.x=(mid*P2.X*P2.Z)%p;
	//nP.y=(mid*P2.Y)%p;

	
	return 0;
}

int scalarMultiplicationAd(affinePoint P, affinePoint &nP, int *item, int number, int Atomicity)
{
	
	projectPoint P2;
	P2.X=mirvar(0); P2.Y=mirvar(0);P2.Z=mirvar(1);
	
	nres(P2.Z,P2.Z);
	
	copy(P.x,P2.X); //P2.X=P.x; 
	copy(P.y,P2.Y); //P2.Y=P.y;
	//projectPoint P2;
	////Big X(0),Y(0),Z(0);
	//P2.X=P.x;
	//P2.Y=P.y;
	//P2.Z=1;
	
	big mid=mirvar(0);
	big mid1=mirvar(0);
	big one=mirvar(1);
	nres(one,one);

	big d=mirvar(0);
	nres_modmult(P.x,P.x,d);
	nres_modmult(P.y,P.y,mid);
	nres_modmult(d,mid,mid1);

	nres_modadd(d,mid,d);
	nres_modsub(d,mid1,d);
	//  d=(P.x*P.x+P.y*P.y-(P.x*P.x*P.y*P.y))%p;


	affinePoint negativeP;
	negativeP.x=mirvar(0);
	negativeP.y=mirvar(0);
	copy(P.y,negativeP.y);
	nres_negate(P.x, negativeP.x);

	//affinePoint negativeP;
	//negativeP.x=(-P.x)%p;
	//negativeP.y=P.y;
	switch (Atomicity)
	{
	case 0:
	for(int i=number-2; i>-1; i--)
	{
		InvEdwardsDouble(P2, P2,d);
		if(item[i]==1)
		{
			InvEdwardsMixedAddition(P2,P, P2,d);
		}
		else if(item[i]==-1)
		{
			InvEdwardsMixedAddition(P2,negativeP, P2, d);
		}		
	}
	break;
	case 1:
		for (int i = number - 2; i>-1; i--)
		{
			InvEdwardsDoubleAtomicity(P2, P2, d);
			if (item[i] == 1)
			{
				InvEdwardsMixedAdditionAtomicity(P2, P, P2, d);
			}
			else if (item[i] == -1)
			{
				InvEdwardsMixedAdditionAtomicity(P2, negativeP, P2, d);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	
	if(compare(P2.X,one)!=0)
	{			
		nres_moddiv(one,P2.X,mid);
	}
	else
	{
		copy(one,mid);
	}
	//mid=inverse(P2.X,p);
	nres_modmult(mid,P2.Z,nP.x);
	//nP.x=(mid*P2.Z)%p;

	if(compare(P2.Y,one)!=0)
	{			
		nres_moddiv(one,P2.Y,mid);
	}
	else
	{
		copy(one,mid);
	}
	
	nres_modmult(mid,P2.Z,nP.y);
	
	//mid=inverse(P2.Y,p);
	//nP.y=(mid*P2.Z)%p;

	return 0;
}

 
int scalarMultiplicationAd(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, int Atomicity)
{

	
	projectPoint P2;
	P2.X=mirvar(0); P2.Y=mirvar(0);P2.Z=mirvar(1);
	
	nres(P2.Z,P2.Z);
	
	copy(P.x,P2.X); //P2.X=P.x; 
	copy(P.y,P2.Y); //P2.Y=P.y;
	////projectPoint P2;
	//Big X(0),Y(0),Z(0);
	//P2.X=P.x;
	//P2.Y=P.y;
	//P2.Z=1;

	big mid=mirvar(0);
	big mid1=mirvar(0);
	big one=mirvar(1);
	nres(one,one);

	big d=mirvar(0);
	nres_modmult(P.x,P.x,d);
	nres_modmult(P.y,P.y,mid);
	nres_modmult(d,mid,mid1);

	nres_modadd(d,mid,d);
	nres_modsub(d,mid1,d);
	//  d=(P.x*P.x+P.y*P.y-(P.x*P.x*P.y*P.y))%p;


	affinePoint negativeP;
	negativeP.x=mirvar(0);
	negativeP.y=mirvar(0);
	copy(P.y,negativeP.y);
	nres_negate(P.x, negativeP.x);

	////affinePoint negativeP;
	////negativeP.x=(-P.x)%p;
	////negativeP.y=P.y;

	////Big d=(P.x*P.x+P.y*P.y-(P.x*P.x*P.y*P.y))%p;
	switch (Atomicity)
	{
	case 0:
		if (number == 0)
		{
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				InvEdwardsDouble(P2, P2, d);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				InvEdwardsTriple(P2, P2, d);
			}
		}
		else
		{
			for (int j = 0; j < item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				InvEdwardsDouble(P2, P2, d);
			}
			for (int j = 0; j < item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				InvEdwardsTriple(P2, P2, d);
			}

			for (int i = number - 1; i > 0; i--)
			{
				if (item[i].sign == 1)
				{
					InvEdwardsMixedAddition(P2, P, P2, d);
				}
				else if (item[i].sign == -1)
				{
					InvEdwardsMixedAddition(P2, negativeP, P2, d);
				}
				for (int j = 0; j < item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					InvEdwardsDouble(P2, P2, d);
				}
				for (int j = 0; j < item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					InvEdwardsTriple(P2, P2, d);
				}
			}
			if (item[0].sign == 1)
			{
				InvEdwardsMixedAddition(P2, P, P2, d);
			}
			else if (item[0].sign == -1)
			{
				InvEdwardsMixedAddition(P2, negativeP, P2, d);
			}
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				InvEdwardsDouble(P2, P2, d);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				InvEdwardsTriple(P2, P2, d);
			}
		}


		break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}
		}
		else
		{
			for (int j = 0; j < item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j < item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}

			for (int i = number - 1; i > 0; i--)
			{
				if (item[i].sign == 1)
				{
					InvEdwardsMixedAdditionAtomicity(P2, P, P2, d);
				}
				else if (item[i].sign == -1)
				{
					InvEdwardsMixedAdditionAtomicity(P2, negativeP, P2, d);
				}
				for (int j = 0; j < item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					InvEdwardsDoubleAtomicity(P2, P2, d);
				}
				for (int j = 0; j < item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					InvEdwardsTripleAtomicity(P2, P2, d);
				}
			}
			if (item[0].sign == 1)
			{
				InvEdwardsMixedAdditionAtomicity(P2, P, P2, d);
			}
			else if (item[0].sign == -1)
			{
				InvEdwardsMixedAdditionAtomicity(P2, negativeP, P2, d);
			}
			for (int j = 0; j < item[0].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j < item[0].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}
		}
		break;
	default: cout << "error!" << endl;
	}

		if(compare(P2.X,one)!=0)
	{			
		nres_moddiv(one,P2.X,mid);
	}
	else
	{
		copy(one,mid);
	}
	//mid=inverse(P2.X,p);
	nres_modmult(mid,P2.Z,nP.x);
	//nP.x=(mid*P2.Z)%p;

	if(compare(P2.Y,one)!=0)
	{			
		nres_moddiv(one,P2.Y,mid);
	}
	else
	{
		copy(one,mid);
	}
	
	nres_modmult(mid,P2.Z,nP.y);
	//Big mid(0);
	//mid=inverse(P2.X,p);
	//nP.x=(mid*P2.Z)%p;
	//mid=inverse(P2.Y,p);
	//nP.y=(mid*P2.Z)%p;

	
	return 0;
}


int windowScalarMultiplicationJ(affinePoint P, affinePoint &nP, int *item, int number, projectPoint *Pre, int Cmax, int Atomicity)
{
	//*Pre[i]=[i]P,i=1,3,5,7,....Cmax

	projectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1);
	//projectPoint P3;
	//P3.X=mirvar(0); P3.Y=mirvar(0);P3.Z=mirvar(1);
	nres(P2.Z, P2.Z);
	//Big X(0),Y(0),Z(0);
	copy(P.x, P2.X); //P2.X=P.x; 
	copy(P.y, P2.Y); //P2.Y=P.y;
					 //P2.Z=1;


	switch (Atomicity)
	{
	case 0:
	for (int i = number - 2; i>-1; i--)
	{
		JacobiDouble(P2, P2);

		if (item[i] > 0)
		{
			JacobiAddition(P2, Pre[item[i]], P2);
		}
		else if (item[i] <0)
		{
			int index = item[i] * (-1);
			JacobiNegative(Pre[index], Pre[0]);
			JacobiAddition(P2, Pre[0], P2);
		}
	}
	break;
	case 1:
		for (int i = number - 2; i>-1; i--)
		{
			JacobiDoubleAtomicity(P2, P2);

			if (item[i] > 0)
			{
				JacobiAdditionAtomicity(P2, Pre[item[i]], P2);
			}
			else if (item[i] <0)
			{
				int index = item[i] * (-1);
				JacobiNegativeAtomicity(Pre[index], Pre[0]);
				JacobiAdditionAtomicity(P2, Pre[0], P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}

	big mid = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	nres_modmult(P2.Z, P2.Z, mid);
	nres_modmult(mid, P2.Z, mid);
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(mid, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, mid, mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(mid, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid, P2.X, mid);
	nres_modmult(mid, P2.Z, nP.x);

	return 0;
}


int windowScalarMultiplicationJ(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, projectPoint *Pre, int Cmax, int Atomicity)
{

	projectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1);

	nres(P2.Z, P2.Z);

	int index = item[number - 1].sign;
	//P2=Pre[index]
	copy(Pre[index].X, P2.X); //P2.X=P.x; 
	copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
	copy(Pre[index].Z, P2.Z); 				 //projectPoint P2;


	affinePoint negativeP;
	negativeP.x = mirvar(0);
	negativeP.y = mirvar(0);
	copy(P.x, negativeP.x);
	nres_negate(P.y, negativeP.y);
	////affinePoint negativeP;
	////negativeP.x=P.x;
	////negativeP.y=(-P.y)%p;
	switch (Atomicity)
	{
	case 0:
	if (number == 0)
	{
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			JacobiDouble(P2, P2);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			JacobiTriple(P2, P2);
		}
	}
	else
	{
		for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
		{
			JacobiDouble(P2, P2);
		}
		for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
		{
			JacobiTriple(P2, P2);
		}

		for (int i = number - 1; i>0; i--)
		{
			if (item[i].sign > 0)
			{
				JacobiMixedAddition(P2, P, P2);
			}
			else if (item[i].sign < 0)
			{
				JacobiMixedAddition(P2, negativeP, P2);
			}
			for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
			{
				JacobiDouble(P2, P2);
			}
			for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
			{
				JacobiTriple(P2, P2);
			}
		}
		if (item[0].sign > 0)
		{
			JacobiMixedAddition(P2, P, P2);
		}
		else if (item[0].sign < 0)
		{
			JacobiMixedAddition(P2, negativeP, P2);
		}
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			JacobiDouble(P2, P2);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			JacobiTriple(P2, P2);
		}
	}
	break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}
		}
		else
		{
			for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}

			for (int i = number - 1; i>0; i--)
			{
				if (item[i].sign > 0)
				{
					JacobiMixedAdditionAtomicity(P2, P, P2);
				}
				else if (item[i].sign < 0)
				{
					JacobiMixedAdditionAtomicity(P2, negativeP, P2);
				}
				for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					JacobiDoubleAtomicity(P2, P2);
				}
				for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					JacobiTripleAtomicity(P2, P2);
				}
			}
			if (item[0].sign > 0)
			{
				JacobiMixedAdditionAtomicity(P2, P, P2);
			}
			else if (item[0].sign < 0)
			{
				JacobiMixedAdditionAtomicity(P2, negativeP, P2);
			}
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				JacobiDoubleAtomicity(P2, P2);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				JacobiTripleAtomicity(P2, P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	nres_modmult(P2.Z, P2.Z, mid);
	nres_modmult(mid, P2.Z, mid);
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(mid, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, mid, mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(mid, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid, P2.X, mid);
	nres_modmult(mid, P2.Z, nP.x);
	//Big mid(0);
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p,p);
	////	cout <<":" <<P2.X <<":"<<P2.Y <<":"<<P2.Z<<endl;
	//nP.x=(mid*P2.X*P2.Z)%p;
	//nP.y=(mid*P2.Y)%p;


	return 0;
}


int windowScalarMultiplicationEx(affinePoint P, affinePoint &nP, int *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity)
{
	//*Pre[i]=[i]P,i=1,3,5,7,....Cmax

	extProjectPoint P2;
	big X3, Y3, Z3,T3;
	X3 = mirvar(0), Y3 = mirvar(0), Z3 = mirvar(0), T3 = mirvar(0);
	P2.X = X3; P2.Y = Y3; P2.Z =Z3,P2.T=T3;
	//projectPoint P3;
	//P3.X=mirvar(0); P3.Y=mirvar(0);P3.Z=mirvar(1);
	nres(P2.Z, P2.Z);
	//Big X(0),Y(0),Z(0);
	copy(P.x, P2.X); //P2.X=P.x; 
	copy(P.y, P2.Y); //P2.Y=P.y;
					 //P2.Z=1;
	nres_modmult(P2.X, P2.Y, P2.T);
	//P2.T=P2.X*P2.Y

	big d = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	for (int i = number - 2; i>-1; i--)
	{
		EdwardsDouble(P2, P2, d);

		if (item[i] > 0)
		{
			EdwardsAddition(P2, Pre[item[i]], P2, d);
		}
		else if (item[i] <0)
		{
			int index = item[i] * (-1);
			EdwardsNegative(Pre[index], Pre[0], d);
			EdwardsAddition(P2, Pre[0], P2, d);
		}
	}
	break;
	case 1:
		for (int i = number - 2; i>-1; i--)
		{
			EdwardsDoubleAtomicity(P2, P2, d);

			if (item[i] > 0)
			{
				EdwardsAdditionAtomicity(P2, Pre[item[i]], P2, d);
			}
			else if (item[i] <0)
			{
				int index = item[i] * (-1);
				EdwardsNegativeAtomicity(Pre[index], Pre[0], d);
				EdwardsAdditionAtomicity(P2, Pre[0], P2, d);
			}
		}
		break;
	default: cout << "error!" << endl;
	}


	big mid = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	//nres_modmult(P2.Z, P2.Z, mid);
	//nres_modmult(mid, P2.Z, mid);
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	//if (compare(mid, one) != 0)
	//{
	//	//nres(mid,mid);
	//	nres_moddiv(one, mid, mid);
	//}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	if (compare(P2.Z, one) != 0)
	{
		nres_moddiv(one, P2.Z, mid);
	}
	nres_modmult(mid, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid, P2.X, nP.x);
	//nP.y=Y/Z
	//nP.X=X/Z

	return 0;
}


int windowScalarMultiplicationEx(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity)
{
	extProjectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1), P2.T = mirvar(1);

	nres(P2.Z, P2.Z);

	int index = item[number - 1].sign;
	//P2=Pre[index]
	if (index < 0)
	{
	
		index = -index;
		nres_negate(Pre[index].X, P2.X); //P2.X=P.x; 
		copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
		nres_negate(Pre[index].Z, P2.Z); 				 //projectPoint P2;
		copy(Pre[index].T, P2.T);

	}
	else
	{
		copy(Pre[index].X, P2.X); //P2.X=P.x; 
		copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
		copy(Pre[index].Z, P2.Z); 				 //projectPoint P2;
		copy(Pre[index].T, P2.T);
	}

	extProjectPoint P1;
	P1.X = mirvar(0); P1.Y = mirvar(0); P1.Z = mirvar(1), P1.T = mirvar(1);
	copy(P.x, P1.X);
	copy(P.y, P1.Y);
	copy(P1.Z, P1.Z);
	nres_modmult(P1.X, P1.Y, P1.T);
	nres_moddiv(P1.T, P1.Z, P1.T);
	//P1=affine(P)

	extProjectPoint negativeP;
	negativeP.X = mirvar(0);
	negativeP.Y = mirvar(0);
	negativeP.T = mirvar(0);
	negativeP.Z = mirvar(1);
	copy(P.x, negativeP.X);
	copy(P.y, negativeP.Y);
	nres_negate(negativeP.X, negativeP.X);
	copy(negativeP.Z, negativeP.Z);
	nres_modmult(negativeP.X, negativeP.Y, negativeP.T);
	nres_moddiv(negativeP.T, negativeP.Z, negativeP.T);
	nres_negate(negativeP.T, negativeP.T);
	////affinePoint negativeP;
	////negativeP.x=P.x;
	////negativeP.y=(-P.y)%p;
	big d = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	if (number == 0)
	{
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			EdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			EdwardsTriple(P2, P2, d);
		}
	}
	else
	{
		for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
		{
			EdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
		{
			EdwardsTriple(P2, P2, d);
		}

		for (int i = number - 1; i>0; i--)
		{
			if (item[i].sign > 0) //Cmax>=item[i].sign
			{
				EdwardsAddition(P2, P1, P2, d);
			}
			else if (item[i].sign < 0) //item[i].sign>=-Cmax
			{
				EdwardsAddition(P2, negativeP, P2, d);
			}
			for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
			{
				EdwardsDouble(P2, P2, d);
			}
			for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
			{
				EdwardsTriple(P2, P2, d);
			}
		}
		if (item[0].sign >0)
		{
			EdwardsAddition(P2, P1, P2, d);
		}
		else if (item[0].sign <0)
		{
			EdwardsAddition(P2, negativeP, P2, d);
		}
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			EdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			EdwardsTriple(P2, P2, d);
		}
	}


	break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				EdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				EdwardsTripleAtomicity(P2, P2, d);
			}
		}
		else
		{
			for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				EdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				EdwardsTripleAtomicity(P2, P2, d);
			}

			for (int i = number - 1; i>0; i--)
			{
				if (item[i].sign > 0) //Cmax>=item[i].sign
				{
					EdwardsAdditionAtomicity(P2, P1, P2, d);
				}
				else if (item[i].sign < 0) //item[i].sign>=-Cmax
				{
					EdwardsAdditionAtomicity(P2, negativeP, P2, d);
				}
				for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					EdwardsDoubleAtomicity(P2, P2, d);
				}
				for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					EdwardsTripleAtomicity(P2, P2, d);
				}
			}
			if (item[0].sign >0)
			{
				EdwardsAdditionAtomicity(P2, P1, P2, d);
			}
			else if (item[0].sign <0)
			{
				EdwardsAdditionAtomicity(P2, negativeP, P2, d);
			}
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				EdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				EdwardsTripleAtomicity(P2, P2, d);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	/*nres_modmult(P2.Z, P2.Z, mid);
	nres_modmult(mid, P2.Z, mid);*/
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(P2.Z, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, P2.Z, mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(mid, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid, P2.X, nP.x);
	//Big mid(0);
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p,p);
	////	cout <<":" <<P2.X <<":"<<P2.Y <<":"<<P2.Z<<endl;
	//nP.x=(mid*P2.X*P2.Z)%p;
	//nP.y=(mid*P2.Y)%p;

	return 0;
}


int windowScalarMultiplicationInvE(affinePoint P, affinePoint &nP, int *item, int number, projectPoint *Pre, int Cmax, int Atomicity)
{
	//*Pre[i]=[i]P,i=1,3,5,7,....Cmax

	projectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1);
	//projectPoint P3;
	//P3.X=mirvar(0); P3.Y=mirvar(0);P3.Z=mirvar(1);
	nres(P2.Z, P2.Z);
	//Big X(0),Y(0),Z(0);
	copy(P.x, P2.X); //P2.X=P.x; 
	copy(P.y, P2.Y); //P2.Y=P.y;
					 //P2.Z=1;
	//nres_modmult(P2.X, P2.Y, P2.T);
	//P2.T=P2.X*P2.Y

	big d = mirvar(0), a = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	for (int i = number - 2; i>-1; i--)
	{
		InvEdwardsDouble(P2, P2, d);

		if (item[i] > 0)
		{
			InvEdwardsAddition(P2, Pre[item[i]], P2, d, a);
		}
		else if (item[i] < 0)
		{
			int index = item[i] * (-1);
			InvEdwardsNegative(Pre[index], Pre[0]);
			InvEdwardsAddition(P2, Pre[0], P2, d, a);
		}
	}
	break;
	case 1:
		for (int i = number - 2; i>-1; i--)
		{
			InvEdwardsDoubleAtomicity(P2, P2, d);

			if (item[i] > 0)
			{
				InvEdwardsAdditionAtomicity(P2, Pre[item[i]], P2, d, a);
			}
			else if (item[i] < 0)
			{
				int index = item[i] * (-1);
				InvEdwardsNegativeAtomicity(Pre[index], Pre[0]);
				InvEdwardsAdditionAtomicity(P2, Pre[0], P2, d, a);
			}
		}
		break;
	default: cout << "error!" << endl;
	}

	big mid1 = mirvar(0), mid2 = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	//nres_modmult(P2.Z, P2.Z, mid);
	//nres_modmult(mid, P2.Z, mid);
	////redc(mid, mid);
	nres(one, one);
	////x=mid;
	////cout<<"Z^3:"<<x<<endl;
	//if (compare(mid, one) != 0)
	//{
	//	//nres(mid,mid);
	//	nres_moddiv(one, mid, mid);
	//}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	if (compare(P2.Y, one) != 0)
	{
		nres_moddiv(one, P2.Y, mid1);
	}
	if (compare(P2.X, one) != 0)
	{
		nres_moddiv(one, P2.X, mid2);
	}
	nres_modmult(mid1, P2.Z, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid2, P2.Z, nP.x);
	//nP.y=Z/Y
	//nP.X=Z/X

	return 0;
}


int windowScalarMultiplicationInvE(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, projectPoint *Pre, int Cmax, int Atomicity)
{
	projectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1);

	nres(P2.Z, P2.Z);

	int index = item[number - 1].sign;
	//P2=Pre[index]
	if (index > 0)
	{
		copy(Pre[index].X, P2.X); //P2.X=P.x; 
		copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
		copy(Pre[index].Z, P2.Z); 				 //projectPoint P2;
	}
	else
	{
		index = -index;
		nres_negate(Pre[index].X, P2.X); //P2.X=P.x; 
		copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
		copy(Pre[index].Z, P2.Z); 				 //projectPoint P2;
	}

	projectPoint P1;
	P1.X = mirvar(0); P1.Y = mirvar(0); P1.Z = mirvar(1);
	copy(P1.Z, P1.Z);
	nres_moddiv(P1.Z, P.x, P1.X);
	nres_moddiv(P1.Z, P.y, P1.Y);
	//P1=affine(P)

	projectPoint negativeP;
	negativeP.X = mirvar(0);
	negativeP.Y = mirvar(0);
	negativeP.Z = mirvar(1);
	copy(negativeP.Z, negativeP.Z);
	nres_moddiv(negativeP.Z, P.y, negativeP.Y);
	nres_moddiv(negativeP.Z, P.x, negativeP.X);
	nres_negate(negativeP.X, negativeP.X);
	
	////affinePoint negativeP;
	////negativeP.x=P.x;
	////negativeP.y=(-P.y)%p;
	big d = mirvar(0),a = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	if (number == 0)
	{
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			InvEdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			InvEdwardsTriple(P2, P2, d);
		}
	}
	else
	{
		for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
		{
			InvEdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
		{
			InvEdwardsTriple(P2, P2, d);
		}

		for (int i = number - 1; i>0; i--)
		{
			if (item[i].sign > 0) //Cmax>=item[i].sign
			{
				InvEdwardsAddition(P2, P1, P2, d, a);
			}
			else if (item[i].sign < 0) //item[i].sign>=-Cmax
			{
				InvEdwardsAddition(P2, negativeP, P2, d, a);
			}
			for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
			{
				InvEdwardsDouble(P2, P2, d);
			}
			for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
			{
				InvEdwardsTriple(P2, P2, d);
			}
		}
		if (item[0].sign >0 )
		{
			InvEdwardsAddition(P2, P1, P2, d, a);
		}
		else if (item[0].sign <0)
		{
			InvEdwardsAddition(P2, negativeP, P2, d, a);
		}
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			InvEdwardsDouble(P2, P2, d);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			InvEdwardsTriple(P2, P2, d);
		}
	}
	break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}
		}
		else
		{
			for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}

			for (int i = number - 1; i>0; i--)
			{
				if (item[i].sign > 0) //Cmax>=item[i].sign
				{
					InvEdwardsAdditionAtomicity(P2, P1, P2, d, a);
				}
				else if (item[i].sign < 0) //item[i].sign>=-Cmax
				{
					InvEdwardsAdditionAtomicity(P2, negativeP, P2, d, a);
				}
				for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					InvEdwardsDoubleAtomicity(P2, P2, d);
				}
				for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					InvEdwardsTripleAtomicity(P2, P2, d);
				}
			}
			if (item[0].sign >0)
			{
				InvEdwardsAdditionAtomicity(P2, P1, P2, d, a);
			}
			else if (item[0].sign <0)
			{
				InvEdwardsAdditionAtomicity(P2, negativeP, P2, d, a);
			}
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				InvEdwardsDoubleAtomicity(P2, P2, d);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				InvEdwardsTripleAtomicity(P2, P2, d);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid = mirvar(0), mid1 = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	/*nres_modmult(P2.Z, P2.Z, mid);
	nres_modmult(mid, P2.Z, mid);*/
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(P2.X, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, P2.X, mid);
	}
	if (compare(P2.Y, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, P2.Y, mid1);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	nres_modmult(P2.Z, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(P2.Z, P2.X, nP.x);
	//Big mid(0);
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p,p);
	////	cout <<":" <<P2.X <<":"<<P2.Y <<":"<<P2.Z<<endl;
	//nP.x=(mid*P2.X*P2.Z)%p;
	//nP.y=(mid*P2.Y)%p;

	return 0;
}


int MontgomeryTrick(projectPoint *P, int w, int Atomicity)
{//montgomery Trick Jacobi coordinates to affine coordinates

	big a[33];
	big b[33];
	big c[33];
	int n = (1 << (w - 2)) - 1;
	big one = mirvar(1), d = mirvar(1);
	for (int i = 0; i < 33; i++)
	{
		a[i] = mirvar(0);
		b[i] = mirvar(0);
		c[i] = mirvar(1);

	}
	for (int i = 1; i < n + 1; i++)
	{
		copy(P[2 * i + 1].Z, a[i]);
	}
	copy(a[1], c[1]);

	for (int i = 2; i < n + 1; i++)
	{
		modmult(a[i], c[i - 1], c[i]);
	}

	inverse(c[n], d);

	for (int i = n; i > 1; i--)
	{
		modmult(d, c[i - 1], b[i]);
		modmult(a[i], d, d);
	}
	copy(d, b[1]);//b[i]=1/a[i]=1/P[2 * i + 1].Z


	for (int i = 1; i < n + 1; i++)
	{
		modmult(b[i], b[i], c[i]);//c[i]=b[i]^2
		modmult(c[i], P[2 * i + 1].X, P[2 * i + 1].X); //P[2 * i + 1].X=P[2 * i + 1].X/P[2 * i + 1].Z^2
		modmult(b[i], c[i], c[i]);//c[i]=b[i]^3
		modmult(c[i], P[2 * i + 1].Y, P[2 * i + 1].Y);//P[2 * i + 1].Y=P[2 * i + 1].Y/P[2 * i + 1].Z^3
		copy(one, P[2 * i + 1].Z);
	}

	//Big out=0;
	//for(int i=1;i<n+1;i++)
	//{	
	//	cout <<"Q_"<<2*i+1<<endl;
	//	out =P[2*i+1].X;
	//	cout <<"X:"<<out <<endl;
	//	out =P[2*i+1].Y;
	//	cout <<"Y:"<<out <<endl;
	//	out =P[2*i+1].Z;
	//	cout <<"Z:"<<out <<endl;
	//	 
	//}

	return 0;
}


int MontgomeryTrick(extProjectPoint *P, int w, int Atomicity)
{//montgomery Trick DIK coordinates to affine coordinates

	big a[33];
	big b[33];
	big c[33];
	int n = (1 << (w - 2)) - 1;
	big one = mirvar(1), d = mirvar(1);
	for (int i = 0; i < 33; i++)
	{
		a[i] = mirvar(0);
		b[i] = mirvar(0);
		c[i] = mirvar(1);

	}
	for (int i = 1; i < n + 1; i++)
	{
		copy(P[2 * i + 1].Z, a[i]);
	}
	copy(a[1], c[1]);

	for (int i = 2; i < n + 1; i++)
	{
		modmult(a[i], c[i - 1], c[i]);
	}

	inverse(c[n], d);

	for (int i = n; i > 1; i--)
	{
		modmult(d, c[i - 1], b[i]);
		modmult(a[i], d, d);
	}
	copy(d, b[1]);//b[i]=1/a[i]=1/P[2 * i + 1].Z


	for (int i = 1; i < n + 1; i++)
	{
		//modmult(b[i], b[i], c[i]);//c[i]=b[i]^2
		modmult(b[i], P[2 * i + 1].X, P[2 * i + 1].X); //P[2 * i + 1].X=P[2 * i + 1].X/P[2 * i + 1].Z^2
		//modmult(b[i], c[i], c[i]);//c[i]=b[i]^3
		modmult(b[i], P[2 * i + 1].Y, P[2 * i + 1].Y);//P[2 * i + 1].Y=P[2 * i + 1].Y/P[2 * i + 1].Z^3
		copy(one, P[2 * i + 1].Z);
		copy(one, P[2 * i + 1].T);
	}

	//Big out=0;
	//for(int i=1;i<n+1;i++)
	//{	
	//	cout <<"Q_"<<2*i+1<<endl;
	//	out =P[2*i+1].X;
	//	cout <<"X:"<<out <<endl;
	//	out =P[2*i+1].Y;
	//	cout <<"Y:"<<out <<endl;
	//	out =P[2*i+1].Z;
	//	cout <<"Z:"<<out <<endl;
	//	 
	//}

	return 0;
}


int windowScalarMultiplicationD(affinePoint P, affinePoint &nP, int *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity)
{
	//*Pre[i]=[i]P,i=1,3,5,7,....Cmax

	extProjectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1);
	//projectPoint P3;
	//P3.X=mirvar(0); P3.Y=mirvar(0);P3.Z=mirvar(1);
	nres(P2.Z, P2.Z);
	//Big X(0),Y(0),Z(0);
	copy(P.x, P2.X); //P2.X=P.x; 
	copy(P.y, P2.Y); //P2.Y=P.y;
					 //P2.Z=1;
	nres_modmult(P2.X, P2.Y, P2.T);
	//P2.T=P2.X*P2.Y

	big d = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	for (int i = number - 2; i>-1; i--)
	{
		DIKDoubling(P2, P2);

		if (item[i] > 0)
		{
			DIKAddition(P2, Pre[item[i]], P2);
		}
		else if (item[i] <0)
		{
			int index = item[i] * (-1);
	
			DIKNegative(Pre[index], Pre[0]);
			DIKAddition(P2, Pre[0], P2);
		}
	}
	break;
	case 1:
		for (int i = number - 2; i>-1; i--)
		{
			DIKDoublingAtomicity(P2, P2);

			if (item[i] > 0)
			{
				DIKAdditionAtomicity(P2, Pre[item[i]], P2);
			}
			else if (item[i] <0)
			{
				int index = item[i] * (-1);

				DIKNegativeAtomicity(Pre[index], Pre[0]);
				DIKAdditionAtomicity(P2, Pre[0], P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid1 = mirvar(0), mid2 = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	nres_modmult(P2.Z, P2.Z, mid1);
	//nres_modmult(mid, P2.Z, mid);
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(mid1, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, mid1, mid1);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	if (compare(P2.Z, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, P2.Z, mid2);
	}
	nres_modmult(mid1, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_modmult(mid2, P2.X, nP.x);
	//nP.x=X/Z
	//nP.y=Y/ZZ
	//ZZ=Z^2
	return 0;
}


int windowScalarMultiplicationD(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity)
{
	extProjectPoint P2;
	P2.X = mirvar(0); P2.Y = mirvar(0); P2.Z = mirvar(1), P2.T = mirvar(1);

	nres(P2.Z, P2.Z);

	int index = item[number - 1].sign;
	//P2=Pre[index]
	copy(Pre[index].X, P2.X); //P2.X=P.x; 
	copy(Pre[index].Y, P2.Y); //P2.Y=P.y;
	copy(Pre[index].Z, P2.Z); 				 //projectPoint P2;
	copy(Pre[index].T, P2.T);

	extProjectPoint P1;
	P1.X = mirvar(0); P1.Y = mirvar(0); P1.Z = mirvar(1), P1.T = mirvar(1);
	copy(P1.Z, P1.Z);
	nres_modmult(P.x, P1.Z, P1.X);
	nres_modmult(P.y, P1.Z, P1.Y);
	nres_modmult(P1.Y, P1.Z, P1.Y);

	nres_modmult(P1.Z, P1.Z, P1.T);
	//P1=affine(P)

	extProjectPoint negativeP;
	negativeP.X = mirvar(0);
	negativeP.Y = mirvar(0);
	negativeP.T = mirvar(0);
	negativeP.Z = mirvar(1);
	
	nres_modmult(P.x, negativeP.Z, negativeP.X);
	nres_modmult(negativeP.Z, negativeP.Z, negativeP.T);
	nres_modmult(P.y, negativeP.T, negativeP.Y);
	
	////affinePoint negativeP;
	////negativeP.x=P.x;
	////negativeP.y=(-P.y)%p;
	big d = mirvar(0);
	switch (Atomicity)
	{
	case 0:
	if (number == 0)
	{
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			DIKDoubling(P2, P2);

		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			DIKTripling(P2, P2);
		}
	}
	else
	{
		for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
		{
			DIKDoubling(P2, P2);
		}
		for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
		{
			DIKTripling(P2, P2);
		}

		for (int i = number - 1; i>0; i--)
		{
			if (item[i].sign > 0) //Cmax>=item[i].sign
			{
				DIKAddition(P2, P1, P2);
			}
			else if (item[i].sign < 0) //item[i].sign>=-Cmax
			{
				DIKAddition(P2, negativeP, P2);
			}
			for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
			{
				DIKDoubling(P2, P2);
			}
			for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
			{
				DIKTripling(P2, P2);
			}
		}
		if (item[0].sign >0)
		{
			DIKAddition(P2, P1, P2);
		}
		else if (item[0].sign <0)
		{
			DIKAddition(P2, negativeP, P2);
		}
		for (int j = 0; j<item[0].binaryExponent; j++)
		{
			DIKDoubling(P2, P2);
		}
		for (int j = 0; j<item[0].ternaryExponent; j++)
		{
			DIKTripling(P2, P2);
		}
	}
	break;
	case 1:
		if (number == 0)
		{
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				DIKDoublingAtomicity(P2, P2);

			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				DIKTriplingAtomicity(P2, P2);
			}
		}
		else
		{
			for (int j = 0; j<item[number].binaryExponent - item[number - 1].binaryExponent; j++)
			{
				DIKDoublingAtomicity(P2, P2);
			}
			for (int j = 0; j<item[number].ternaryExponent - item[number - 1].ternaryExponent; j++)
			{
				DIKTriplingAtomicity(P2, P2);
			}

			for (int i = number - 1; i>0; i--)
			{
				if (item[i].sign > 0) //Cmax>=item[i].sign
				{
					DIKAdditionAtomicity(P2, P1, P2);
				}
				else if (item[i].sign < 0) //item[i].sign>=-Cmax
				{
					DIKAdditionAtomicity(P2, negativeP, P2);
				}
				for (int j = 0; j<item[i].binaryExponent - item[i - 1].binaryExponent; j++)
				{
					DIKDoublingAtomicity(P2, P2);
				}
				for (int j = 0; j<item[i].ternaryExponent - item[i - 1].ternaryExponent; j++)
				{
					DIKTriplingAtomicity(P2, P2);
				}
			}
			if (item[0].sign >0)
			{
				DIKAdditionAtomicity(P2, P1, P2);
			}
			else if (item[0].sign <0)
			{
				DIKAdditionAtomicity(P2, negativeP, P2);
			}
			for (int j = 0; j<item[0].binaryExponent; j++)
			{
				DIKDoublingAtomicity(P2, P2);
			}
			for (int j = 0; j<item[0].ternaryExponent; j++)
			{
				DIKTriplingAtomicity(P2, P2);
			}
		}
		break;
	default: cout << "error!" << endl;
	}
	big mid = mirvar(0);
	big one = mirvar(1);

	//Big x=one; 

	/*nres_modmult(P2.Z, P2.Z, mid);
	nres_modmult(mid, P2.Z, mid);*/
	//redc(mid, mid);
	nres(one, one);
	//x=mid;
	//cout<<"Z^3:"<<x<<endl;
	if (compare(P2.Z, one) != 0)
	{
		//nres(mid,mid);
		nres_moddiv(one, P2.Z, mid);
	}
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p);
	
	nres_modmult(mid, P2.X, nP.x);
	nres_modmult(mid, P2.Y, nP.y);//nP.y=(mid*P2.Y)%p;
	nres_moddiv(nP.y, P2.Z, nP.y);
	//Big mid(0);
	//mid=inverse((P2.Z*P2.Z*P2.Z)%p,p);
	////	cout <<":" <<P2.X <<":"<<P2.Y <<":"<<P2.Z<<endl;
	//nP.x=(mid*P2.X*P2.Z)%p;
	//nP.y=(mid*P2.Y)%p;

	return 0;
}