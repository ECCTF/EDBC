#include "pointOperation.h"


int JacobiMixedAddition(projectPoint P1, affinePoint P2, projectPoint &P3)
{
	big R4=mirvar(0),R5=mirvar(0),A=mirvar(0),B=mirvar(0),M=mirvar(0),P3Z=mirvar(0),P3X=mirvar(0);//,P3Y=mirvar(0);
	nres_modmult(P1.Z,P1.Z,R5);
	//Big R5=P1.Z*P1.Z%p;//Z_1^2
	nres_modmult(P2.x, R5,B);
	nres_modsub(B,P1.X,B);//X_2Z_1^2-X_1
	// B= (P2.x*R5-P1.X)%p;//X_2Z_1^2-X_1
	nres_modmult(P1.Z,B,A);
	nres_modadd(A,A,P3Z);//P3.Z=(2*P1.Z*B)%p;

	//P3.Z=(2*P1.Z*B)%p;

	nres_modmult(P1.Z,P2.y,A);
	nres_modmult(A,R5,A);
	nres_modsub(A,P1.Y,A);
	nres_modadd(A,A,A);//2(Z_1^3Y_2-Y_1)

	//Big A=2*(P2.y*P1.Z*R5-P1.Y)%p;//2(Z_1^3Y_2-Y_1)
	nres_modadd(B,B,R5);
	nres_modmult(R5,R5,R5);

	//R5=(4*B*B)%p;//4B^2
	nres_modmult(B,R5,R4);
	//Big R4=B*R5;//4B^3

	nres_modmult(P1.X,R5,R5);
	//R5=R5*P1.X%p;//4B^2X_1
	nres_modmult(A,A,B);
	//B=A*A%p;
	nres_modsub(B,R4,M);
	nres_modsub(M,R5,M);
	nres_modsub(M,R5,P3X);

	//P3.X=(B-R4-2*R5)%p;
	nres_modsub(R5,P3X,R5);
	nres_modmult(A,R5,A);

	nres_modmult(P1.Y,R4,R4);
	nres_modadd(R4,R4,R4);
	nres_modsub(A,R4,P3.Y);//Y is not right
	//P3.Y=(A*(R5-P3.X)-2*P1.Y*R4)%p;

	copy(P3X,P3.X);
		//copy(P3Y,P3.Y);

	copy(P3Z,P3.Z);

	

	return 0;
}

int JacobiAddition(projectPoint P1, projectPoint P2, projectPoint &P3)
{
	big R4 = mirvar(0), R5 = mirvar(0), A = mirvar(0), B = mirvar(0), M = mirvar(0), P3Z = mirvar(0), P3X = mirvar(0);//,P3Y=mirvar(0);
	big Z1Z1 = mirvar(0), Z2Z2 = mirvar(0), U1 = mirvar(0), U2 = mirvar(0), S1 = mirvar(0), S2 = mirvar(0);
	big H = mirvar(0), HH = mirvar(0), HHH = mirvar(0), r = mirvar(0), V = mirvar(0);
	
	nres_modmult(P1.Z, P1.Z, Z1Z1);
	//Z1Z1=Z1^2
	nres_modmult(P2.Z, P2.Z, Z2Z2);
	//Z2Z2=Z2^2
	nres_modmult(P1.X, Z2Z2, U1);
	//U1=X1*Z2Z2
	nres_modmult(P2.X, Z1Z1, U2);
	//U2=X2*Z1Z1
	nres_modmult(P1.Y, P2.Z, S1);
	nres_modmult(S1, Z2Z2, S1);
	//S1=Y1*Z2*Z2Z2
	nres_modmult(P2.Y, P1.Z, S2);
	nres_modmult(S2, Z1Z1, S2);
	//S2=Y2*Z1*Z1Z1
	nres_modsub(U2, U1, H);
	//H=U2-U1
	nres_modmult(H, H, HH);
	nres_modmult(H, HH, HHH);
	//HH=H^2, HHH=H*HH
	nres_modsub(S2, S1, r);
	nres_modmult(U1, HH, V);
	//r=S2-S1, V=U1*HH

	nres_modmult(r, r, A);
	nres_modsub(A, HHH, A);
	nres_modadd(V, V, B);
	nres_modsub(A, B, P3.X);
	//X3=r^2-HHH-2V

	nres_modsub(V, P3.X, A);
	nres_modmult(r, A, A);
	nres_modmult(S1, HHH, B);
	nres_modsub(A, B, P3.Y);
	//Y3=r*(V-X3)-S1*HHH

	nres_modmult(P1.Z, P2.Z, A);
	nres_modmult(A, H, P3.Z);
	//Z3=Z1*Z2*H

	return 0;
}


int JacobiAdditionSubtraction(projectPoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4)
{
	//12M+6S
	//P3=P1+p2,P4=P1-P2
	big R4 = mirvar(0), R5 = mirvar(0), A = mirvar(0), B = mirvar(0), M = mirvar(0), P3Z = mirvar(0), P3X = mirvar(0);//,P3Y=mirvar(0);
	big Z1Z1 = mirvar(0), Z2Z2 = mirvar(0), U1 = mirvar(0), U2 = mirvar(0), S1 = mirvar(0), S2 = mirvar(0);
	big H = mirvar(0), HH = mirvar(0), HHH = mirvar(0), r = mirvar(0), V = mirvar(0);
	big S21 = mirvar(0), r1 = mirvar(0);

	nres_modmult(P1.Z, P1.Z, Z1Z1);
	//Z1Z1=Z1^2
	nres_modmult(P2.Z, P2.Z, Z2Z2);
	//Z2Z2=Z2^2
	nres_modmult(P1.X, Z2Z2, U1);
	//U1=X1*Z2Z2
	nres_modmult(P2.X, Z1Z1, U2);
	//U2=X2*Z1Z1
	nres_modmult(P1.Y, P2.Z, S1);
	nres_modmult(S1, Z2Z2, S1);
	//S1=Y1*Z2*Z2Z2
	nres_modmult(P2.Y, P1.Z, S2);
	nres_modmult(S2, Z1Z1, S2);
	//S2=Y2*Z1*Z1Z1
	nres_modsub(U2, U1, H);
	//H=U2-U1
	nres_modmult(H, H, HH);
	nres_modmult(H, HH, HHH);
	//HH=H^2, HHH=H*HH
	nres_modsub(S2, S1, r);
	nres_modmult(U1, HH, V);
	//r=S2-S1, V=U1*HH
	nres_modsub(R4, S2, S21);
	//S21=-S2
	nres_modsub(S21, S1, r1);
	//r1=S21-S1


	nres_modmult(r, r, A);
	nres_modsub(A, HHH, A);
	nres_modadd(V, V, B);
	nres_modsub(A, B, P3.X);
	//X3=r^2-HHH-2V,X3=P1.X+P2.X

	nres_modsub(V, P3.X, A);
	nres_modmult(r, A, A);
	nres_modmult(S1, HHH, B);
	nres_modsub(A, B, P3.Y);
	//Y3=r*(V-X3)-S1*HHH,Y3=P1.Y+P2.Y
	nres_modsub(V, P4.X, A);
	nres_modmult(r1, A, A);
	nres_modsub(A, B, P4.Y);
	//Y4=r1*(V-X4)-S1*HHH,Y4=P1.Y-P2.Y

	nres_modmult(P1.Z, P2.Z, A);
	nres_modmult(A, H, P3.Z);
	//Z3=Z1*Z2*H,Z3=P1.Z+P2.Z
	copy(P3.Z, P4.Z);

	nres_modmult(r1, r1, A);
	nres_modsub(A, HHH, A);
	nres_modadd(V, V, B);
	nres_modsub(A, B, P4.X);
	//X4=r1^2-HHH-2V,X4=P1.X-P2.X


	return 0;
}


int JacobiNegative(projectPoint P1, projectPoint &P3)
{
	big R4 = mirvar(0), R5 = mirvar(0), A = mirvar(0), B = mirvar(0), M = mirvar(0), P3Z = mirvar(0), P3X = mirvar(0);//,P3Y=mirvar(0);
	
	copy(P1.X, P3.X);
	//X3=X1
	nres_modsub(A, P1.Y, P3.Y);
	//Y3=-Y1
	copy(P1.Z, P3.Z);
	//Z3=Z1

	return 0;
}


int JacobiNegative(affinePoint P1, affinePoint &P3)
{
	big R4 = mirvar(0), R5 = mirvar(0), A = mirvar(0), B = mirvar(0), M = mirvar(0), P3Z = mirvar(0), P3X = mirvar(0);//,P3Y=mirvar(0);

	copy(P1.x, P3.x);
	//X3=X1
	nres_modsub(A, P1.y, P3.y);
	//Y3=-Y1
	//copy(P1.Z, P3.Z);
	//Z3=Z1

	return 0;
}


int JacobiDouble(projectPoint P1, projectPoint &P3)
{
	big A=mirvar(0),B=mirvar(0),C=mirvar(0),M=mirvar(0),N=mirvar(0),P3Z=mirvar(0);;
	nres_modmult(P1.Y,P1.Z, A);
	nres_modadd(A,A,P3Z);




	//P3.Z=(2*P1.Z*P1.Y)%p;//2Y_1Z_1

	nres_modmult(P1.Y,P1.Y, C);
	nres_modadd(C,C,C);
	//Big C=(2*P1.Y*P1.Y)%p;//2Y_1^2


	nres_modmult(P1.Z,P1.Z, A);
	//Big A=(P1.Z*P1.Z)%p;

	nres_modmult(P1.X,C, B);
	nres_modadd(B,B,B);
	//Big B=(2*C*P1.X)%p;

	nres_modadd(P1.X,A,M);
	nres_modsub(P1.X,A,N);
	nres_modmult(M,N,A);

	nres_modadd(A,A,M);
	nres_modadd(M,A,A);
	//A=(3*(P1.X-A)*(P1.X+A))%p;


	nres_modmult(A,A,M);
	nres_modadd(B,B,N);
	nres_modsub(M,N,P3.X);
	//P3.X=(A*A-2*B)%p;

	nres_modmult(C,C,C);
	nres_modadd(C,C,C);
	nres_modsub(B,P3.X,B);
	nres_modmult(A,B,A);

	nres_modsub(A,C,P3.Y);
	copy(P3Z,P3.Z);
	//P3.Y=(A*(B-P3.X)-2*C*C)%p;
	return 0;
}

int JacobiTriple(projectPoint P1, projectPoint &P3)
{
	big A=mirvar(0),B=mirvar(0),C=mirvar(0),D=mirvar(0),mid=mirvar(0),M=mirvar(0);
	nres_modmult(P1.Z,P1.Z, A);
	//Big A=(P1.Z*P1.Z)%p;	//Z1^2

	nres_modmult(P1.Y,P1.Y, mid);
	nres_modadd(mid,mid,mid);
	nres_modadd(mid,mid,mid);	
	//Big mid=(4*P1.Y*P1.Y)%p;//4Y_1^2


	nres_modadd(P1.X,A,B);
	nres_modsub(P1.X,A,C);
	nres_modmult(B,C,C);
	nres_modadd(C,C,B);
	nres_modadd(B,C,C);
	//	Big C=(3*(P1.X-A)*(P1.X+A))%p;//C



	nres_modmult(P1.X,mid,D);
	nres_modadd(D,D,B);
	nres_modadd(B,D,D);
	nres_modmult(C,C,B);
	nres_modsub(D,B,D);
	//Big D=(3*P1.X*mid-C*C)%p;//D

	nres_modmult(C,D,B);
	nres_modadd(B,B,A);
	//A=(2*C*D)%p;//A


	nres_modmult(P1.Z,D,B);
	nres_modadd(B,B,P3.Z);
	//P3.Z=(2*P1.Z*D)%p;//Z3

	nres_modmult(mid,mid,B);
	//Big B=(mid*mid)%p;//B


	nres_modmult(D,D,C);
	//C=D*D%p;




	nres_modsub(B,A,M);
	nres_modmult(mid, M, mid);
	nres_modmult(C, P1.X, M);	
	nres_modadd(M,mid,M);
	nres_modadd(M,M,M);
	nres_modadd(M,M,P3.X);
	//P3.X=(4*mid*(B-A)+4*P1.X*C)%p;


	nres_modmult(C,D,D);
	nres_modadd(B,B,C);
	nres_modsub(C,A,C);
	nres_modsub(A,B,B);
	nres_modmult(B,C,A);
	nres_modsub(A,D,B);
	nres_modmult(P1.Y,B,A);
	nres_modadd(A,A,A);
	nres_modadd(A,A,A);
	nres_modadd(A,A,P3.Y);
	//P3.Y=(8*P1.Y*((A-B)*(2*B-A)-C*D))%p;

	return 0;
}






int EdwardsAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d)
{
	//8M
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0), A = mirvar(0), D = mirvar(0);
	big F = mirvar(0), G= mirvar(0);

	nres_modsub(P1.Y, P1.X, H);
	//H=Y1-X1
	nres_modsub(P2.Y, P2.X, I); 
	//I=Y2-X2
	nres_modmult(H, I, A);
	//A=H*I=(Y1-X1)(Y2-X2)
	nres_modadd(P1.Y, P1.X, H);
	//H=Y1+X1
	nres_modadd(P2.Y, P2.X, I);
	//I=Y2+X2
	nres_modmult(H, I, B);
	//B=H*I=(Y1+X1)(Y2+X2)
	nres_modmult(P1.Z, P2.T, H);
	nres_modadd(H, H, C);//C=2Z1*T2
	nres_modmult(P1.T, P2.Z, H);
	nres_modadd(H, H, D);//D=2T1*Z2
	nres_modadd(D, C, E);//E=D+C
	nres_modsub(B, A, F);//F=B-A
	nres_modadd(B, A, G);//G=B+A
	nres_modsub(D, C, H);//H=D-C
	nres_modmult(E, F, P3.X);//P3.X=E*F
	nres_modmult(G, H, P3.Y);//P3.Y=G*H
	nres_modmult(E, H, P3.T);//P3.T=E*H
	nres_modmult(F, G, P3.Z);//P3.Z=F*G
	

	return 0;
}

int EdwardsMixedAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d)
{
	//P2.Z=1,7M
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0), A = mirvar(0), D = mirvar(0);
	big F = mirvar(0), G = mirvar(0);

	nres_modsub(P1.Y, P1.X, H);
	//H=Y1-X1
	nres_modsub(P2.Y, P2.X, I);
	//I=Y2-X2
	nres_modmult(H, I, A);
	//A=H*I=(Y1-X1)(Y2-X2)
	nres_modadd(P1.Y, P1.X, H);
	//H=Y1+X1
	nres_modadd(P2.Y, P2.X, I);
	//I=Y2+X2
	nres_modmult(H, I, B);
	//B=H*I=(Y1+X1)(Y2+X2)
	nres_modmult(P1.Z, P2.T, H);
	nres_modadd(H, H, C);//C=2Z1*T2
	//nres_modmult(P1.T, P2.Z, H);
	nres_modadd(P1.T, P1.T, D);//D=2T1*Z2
	nres_modadd(D, C, E);//E=D+C
	nres_modsub(B, A, F);//F=B-A
	nres_modadd(B, A, G);//G=B+A
	nres_modsub(D, C, H);//H=D-C
	nres_modmult(E, F, P3.X);//P3.X=E*F
	nres_modmult(G, H, P3.Y);//P3.Y=G*H
	nres_modmult(E, H, P3.T);//P3.T=E*H
	nres_modmult(F, G, P3.Z);//P3.Z=F*G


	return 0;
}


int EdwardsBothMixedAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d)
{
	//P2.Z=1,P1.Z=1,6M
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0), A = mirvar(0), D = mirvar(0);
	big F = mirvar(0), G = mirvar(0);

	nres_modsub(P1.Y, P1.X, H);
	//H=Y1-X1
	nres_modsub(P2.Y, P2.X, I);
	//I=Y2-X2
	nres_modmult(H, I, A);
	//A=H*I=(Y1-X1)(Y2-X2)
	nres_modadd(P1.Y, P1.X, H);
	//H=Y1+X1
	nres_modadd(P2.Y, P2.X, I);
	//I=Y2+X2
	nres_modmult(H, I, B);
	//B=H*I=(Y1+X1)(Y2+X2)
	//nres_modmult(P1.Z, P2.T, H);
	nres_modadd(P2.T, P2.T, C);//C=2Z1*T2
						 //nres_modmult(P1.T, P2.Z, H);
	nres_modadd(P1.T, P1.T, D);//D=2T1*Z2
	nres_modadd(D, C, E);//E=D+C
	nres_modsub(B, A, F);//F=B-A
	nres_modadd(B, A, G);//G=B+A
	nres_modsub(D, C, H);//H=D-C
	nres_modmult(E, F, P3.X);//P3.X=E*F
	nres_modmult(G, H, P3.Y);//P3.Y=G*H
	nres_modmult(E, H, P3.T);//P3.T=E*H
	nres_modmult(F, G, P3.Z);//P3.Z=F*G


	return 0;
}


int EdwardsDouble(projectPoint P1, projectPoint &P3, big d)
{
	//7M
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1),I= mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E); 
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G

	return 0;
}


int EdwardsMixedDouble(projectPoint P1, projectPoint &P3, big d)
{
	//7M
	big A = mirvar(0), B = mirvar(0), C = mirvar(1), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	//nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G

	return 0;
}


int EdwardsDouble(extProjectPoint P1, projectPoint &P3, big d)
{
	//7M
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G

	return 0;
}


int EdwardsMixedDouble(extProjectPoint P1, projectPoint &P3, big d)
{
	//7M
	big A = mirvar(0), B = mirvar(0), C = mirvar(1), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	//nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G

	return 0;
}



int EdwardsDouble(extProjectPoint P1, extProjectPoint &P3, big d)
{
	//8M
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G
	nres_modmult(E, H, P3.T);
	//P3.T=E*H

	return 0;
}


int EdwardsMixedDouble(extProjectPoint P1, extProjectPoint &P3, big d)
{
	//8M
	big A = mirvar(0), B = mirvar(0), C = mirvar(1), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	//nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G
	nres_modmult(E, H, P3.T);
	//P3.T=E*H

	return 0;
}



int EdwardsDouble(projectPoint P1, extProjectPoint &P3, big d)
{
	//8M
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G
	nres_modmult(E, H, P3.T);
	//P3.T=E*H

	return 0;
}


int EdwardsMixedDouble(projectPoint P1, extProjectPoint &P3, big d)
{
	//8M
	big A = mirvar(0), B = mirvar(0), C = mirvar(1), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);

	nres_modmult(P1.X, P1.X, A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y, P1.Y, B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	//nres_modmult(P1.Z, P1.Z, C);
	nres_modadd(C, C, C);
	//Big C=2Z1^2;

	nres_modmult(a, A, D);
	//D=aA;

	nres_modadd(P1.X, P1.Y, I);
	nres_modmult(I, I, G);
	nres_modsub(G, A, H);
	nres_modsub(H, B, E);
	//E=(P1.X+P1.Y)^2-A-B
	nres_modadd(D, B, G);//G=D+B
	nres_modsub(G, C, F);//F=G-C
	nres_modsub(D, B, H);//H=D-B

	nres_modmult(E, F, P3.X);
	//P3.X=E*F
	nres_modmult(G, H, P3.Y);
	//P3.Y=G*H
	nres_modmult(F, G, P3.Z);
	//P3.Z=F*G
	nres_modmult(E, H, P3.T);
	//P3.T=E*H

	return 0;
}



int EdwardsTriple(projectPoint P1, projectPoint &P3, big d)
{
	//9M+3S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0),xB= mirvar(0);
	big yB = mirvar(0), AA = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	nres_modmult(P1.Z, P1.Z, B);
	nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(H, P1.X, E);
	nres_modmult(E, F, P3.X);
	//X3=X1*(yB+AA)*F
	nres_modsub(xB,AA, H);
	nres_modmult(H, P1.Y, E);
	nres_modmult(E, G, P3.Y);
	//Y3=Y1*(xB-AA)*G
	nres_modmult(P1.Z, F, E);
	nres_modmult(E, G, P3.Z);
	//Z3=Z1*F*G

	return 0;
}


int EdwardsMixedTriple(projectPoint P1, projectPoint &P3, big d)
{
	//8M+2S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(1), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(2);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	//nres_modmult(P1.Z, P1.Z, B);
	//nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(H, P1.X, E);
	nres_modmult(E, F, P3.X);
	//X3=X1*(yB+AA)*F
	nres_modsub(xB, AA, H);
	nres_modmult(H, P1.Y, E);
	nres_modmult(E, G, P3.Y);
	//Y3=Y1*(xB-AA)*G
	//nres_modmult(P1.Z, F, E);
	nres_modmult(F, G, P3.Z);
	//Z3=Z1*F*G

	return 0;
}


int EdwardsTriple(extProjectPoint P1, projectPoint &P3, big d)
{
	//9M+3S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	nres_modmult(P1.Z, P1.Z, B);
	nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(H, P1.X, E);
	nres_modmult(E, F, P3.X);
	//X3=X1*(yB+AA)*F
	nres_modsub(xB, AA, H);
	nres_modmult(H, P1.Y, E);
	nres_modmult(E, G, P3.Y);
	//Y3=Y1*(xB-AA)*G
	nres_modmult(P1.Z, F, E);
	nres_modmult(E, G, P3.Z);
	//Z3=Z1*F*G

	return 0;
}


int EdwardsMixedTriple(extProjectPoint P1, projectPoint &P3, big d)
{
	//8M+2S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(1), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(2);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	//nres_modmult(P1.Z, P1.Z, B);
	//nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(H, P1.X, E);
	nres_modmult(E, F, P3.X);
	//X3=X1*(yB+AA)*F
	nres_modsub(xB, AA, H);
	nres_modmult(H, P1.Y, E);
	nres_modmult(E, G, P3.Y);
	//Y3=Y1*(xB-AA)*G
	//nres_modmult(P1.Z, F, E);
	nres_modmult(F, G, P3.Z);
	//Z3=Z1*F*G

	return 0;
}


int EdwardsTriple(projectPoint P1, extProjectPoint &P3, big d)
{
	//11M+3S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0), xE = mirvar(0), yH = mirvar(0), zF = mirvar(0), zG = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	nres_modmult(P1.Z, P1.Z, B);
	nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(P1.X, H, xE);
	//xE=X1*(yB+AA)
	nres_modsub(xB, AA, H);
	nres_modmult(P1.Y, H, yH);
	//yH=Y1*(xB-AA)
	nres_modmult(P1.Z, F, zF);
	//zF=Z1*F
	nres_modmult(P1.Z, G, zG);
	//zG=Z1*G
	
	nres_modmult(xE, zF, P3.X);
	//X3= xE*zF
	
	nres_modmult(yH, zG, P3.Y);
	//Y3=yH*zG
	
	nres_modmult(zF, zG, P3.Z);
	//Z3=zF*zG

	nres_modmult(xE, yH, P3.T);
	//T3=xE*yH

	return 0;
}


int EdwardsMixedTriple(projectPoint P1, extProjectPoint &P3, big d)
{
	//9M+2S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(2);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0), xE = mirvar(0), yH = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	//nres_modmult(P1.Z, P1.Z, B);
	//nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(P1.X, H, xE);
	//xE=X1*(yB+AA)
	nres_modsub(xB, AA, H);
	nres_modmult(P1.Y, H, yH);
	//yH=Y1*(xB-AA)
	//nres_modmult(P1.Z, F, zF);
	//zF=Z1*F
	//nres_modmult(P1.Z, G, zG);
	//zG=Z1*G

	nres_modmult(xE, F, P3.X);
	//X3= xE*zF

	nres_modmult(yH, G, P3.Y);
	//Y3=yH*zG

	nres_modmult(F, G, P3.Z);
	//Z3=zF*zG

	nres_modmult(xE, yH, P3.T);
	//T3=xE*yH

	return 0;
}


int EdwardsTriple(extProjectPoint P1, extProjectPoint &P3, big d)
{
	//11M+3S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(0);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0), xE = mirvar(0), yH = mirvar(0), zF = mirvar(0), zG = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	nres_modmult(P1.Z, P1.Z, B);
	nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(P1.X, H, xE);
	//xE=X1*(yB+AA)
	nres_modsub(xB, AA, H);
	nres_modmult(P1.Y, H, yH);
	//yH=Y1*(xB-AA)
	nres_modmult(P1.Z, F, zF);
	//zF=Z1*F
	nres_modmult(P1.Z, G, zG);
	//zG=Z1*G

	nres_modmult(xE, zF, P3.X);
	//X3= xE*zF

	nres_modmult(yH, zG, P3.Y);
	//Y3=yH*zG

	nres_modmult(zF, zG, P3.Z);
	//Z3=zF*zG

	nres_modmult(xE, yH, P3.T);
	//T3=xE*yH

	return 0;
}


int EdwardsMixedTriple(extProjectPoint P1, extProjectPoint &P3, big d)
{
	//9M+2S, a=-1
	big A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0), a = mirvar(-1), I = mirvar(0);
	big G = mirvar(0), F = mirvar(0), H = mirvar(2);
	big YY = mirvar(0), aXX = mirvar(0), Ap = mirvar(0), xB = mirvar(0);
	big yB = mirvar(0), AA = mirvar(0), xE = mirvar(0), yH = mirvar(0);

	nres_modmult(P1.Y, P1.Y, YY);
	//Y_1^2

	nres_modmult(P1.X, P1.X, B);
	//X_1^2
	nres_modsub(C, B, aXX);
	//aX_1^2
	nres_modadd(YY, aXX, Ap);
	//Ap=YY+aXX
	//nres_modmult(P1.Z, P1.Z, B);
	//nres_modadd(B, B, H);
	nres_modsub(H, Ap, F);
	nres_modadd(F, F, B);
	//B=2(2Z1^2-Ap)
	nres_modmult(aXX, B, xB);
	//xB=aXX*B
	nres_modmult(YY, B, yB);
	//yB=YY*B
	nres_modsub(YY, aXX, F);
	nres_modmult(Ap, F, AA);
	//AA=Ap*(YY-aXX)
	nres_modsub(AA, yB, F);
	//F=AA-yB
	nres_modadd(AA, xB, G);
	//G=AA+xB
	nres_modadd(yB, AA, H);
	nres_modmult(P1.X, H, xE);
	//xE=X1*(yB+AA)
	nres_modsub(xB, AA, H);
	nres_modmult(P1.Y, H, yH);
	//yH=Y1*(xB-AA)
	//nres_modmult(P1.Z, F, zF);
	//zF=Z1*F
	//nres_modmult(P1.Z, G, zG);
	//zG=Z1*G

	nres_modmult(xE, F, P3.X);
	//X3= xE*zF

	nres_modmult(yH, G, P3.Y);
	//Y3=yH*zG

	nres_modmult(F, G, P3.Z);
	//Z3=zF*zG

	nres_modmult(xE, yH, P3.T);
	//T3=xE*yH

	return 0;
}



int EdwardsNegative(affinePoint P1, affinePoint &P3, big d)
{
	// a=-1
	big A = mirvar(0);
	nres_modsub(A, P1.x, P3.x);
	//X3=-X1
	nres_modadd(A, P1.y, P3.y);
	//Y3=Y1

	return 0;
}


int EdwardsNegative(projectPoint P1, projectPoint &P3, big d)
{
	// a=-1
	big A = mirvar(0);
	nres_modsub(A, P1.X, P3.X);
	//X3=-X1
	nres_modadd(A, P1.Y, P3.Y);
	//Y3=Y1
	nres_modadd(A, P1.Z, P3.Z);
	//Z3=Z1
	return 0;
}


int EdwardsNegative(extProjectPoint P1, extProjectPoint &P3, big d)
{
	// a=-1
	big A = mirvar(0);
	nres_modsub(A, P1.X, P3.X);
	//X3=-X1
	nres_modadd(A, P1.Y, P3.Y);
	//Y3=Y1
	nres_modadd(A, P1.Z, P3.Z);
	//Z3=Z1
	nres_modsub(A, P1.T, P3.T);
	//T3=-T1
	return 0;
}


int EdwardsAdditionSubtraction(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4, big d)
{
	//13M
	//P4=P1-P2,P3=P1+P2
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0), A = mirvar(0), D = mirvar(0);
	big F = mirvar(0), G = mirvar(0), A1 = mirvar(0), B1 = mirvar(0), C1 = mirvar(0), D1 = mirvar(0), E1 = mirvar(0);
	big F1 = mirvar(0), G1 = mirvar(0), H1 = mirvar(0), zero = mirvar(0);

	nres_modsub(P1.Y, P1.X, H);
	//H=Y1-X1
	nres_modsub(P2.Y, P2.X, I);
	//I=Y2-X2
	nres_modadd(P2.Y, P2.X, A1);
	//A1=Y2+X2
	nres_modmult(H, A1, A1);
	//A1=H*A1=(Y1-X1)(Y2+X2)
	nres_modmult(H, I, A);
	//A=H*I=(Y1-X1)(Y2-X2)
	nres_modadd(P1.Y, P1.X, H);
	//H=Y1+X1
	nres_modadd(P2.Y, P2.X, I);
	//I=Y2+X2
	nres_modmult(H, I, B);
	//B=H*I=(Y1+X1)(Y2+X2)
	nres_modsub(P2.Y, P2.X, I);
	//I=Y2-X2
	nres_modmult(H, I, B1);
	//B1=H*I=(Y1+X1)(Y2-X2)
	nres_modmult(P1.Z, P2.T, H);
	nres_modadd(H, H, C);//C=2Z1*T2
	nres_modsub(zero, C, C1);
	//C1=-C
	nres_modmult(P1.T, P2.Z, H);
	nres_modadd(H, H, D);//D=2T1*Z2
	nres_modadd(D, C, E);//E=D+C
	nres_modadd(D, C1, E1);//E1=D+C1
	nres_modsub(B, A, F);//F=B-A
	nres_modsub(B1, A1, F1);//F1=B1-A1
	nres_modadd(B, A, G);//G=B+A
	nres_modadd(B1, A1, G1);//G1=B1+A1
	nres_modsub(D, C, H);//H=D-C
	nres_modsub(D, C1, H1);//H1=D-C1

	nres_modmult(E, F, P3.X);//P3.X=E*F
	nres_modmult(G, H, P3.Y);//P3.Y=G*H
	nres_modmult(E, H, P3.T);//P3.T=E*H
	nres_modmult(F, G, P3.Z);//P3.Z=F*G

	nres_modmult(E1, F1, P4.X);//P4.X=E1*F1
	nres_modmult(G1, H1, P4.Y);//P4.Y=G1*H1
	nres_modmult(F1, G1, P4.Z);//P4.Z=F1*G1
	copy(P3.T, P4.T);
	//P4.T=P3.T

	return 0;
}




int InvEdwardsMixedAddition(projectPoint P1, affinePoint P2, projectPoint &P3, big d)
{

	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);

	//Big A=P1.Z;//Z_1Z_2,Z_2=1
	nres_modmult(P1.Z, P1.Z, B);
	nres_modmult(B, d, B);
	//Big B=(d*A*A)%p;//dA^2,d=1


	nres_modmult(P1.X, P2.x, C);
	//Big C=(P1.X*P2.x)%p;//X_1X_2

	nres_modmult(P1.Y, P2.y, I);
	//Big D=(P1.Y*P2.y)%p;//Y_1Y_2


	nres_modmult(C, I, E);
	//Big E=(C*D)%p;//C*D

	nres_modsub(C, I, H);
	//Big H=C-D;//C-D

	nres_modadd(C, I, P3X);//C+D
	nres_modadd(P1.X, P1.Y, C);
	nres_modadd(P2.x, P2.y, I);
	nres_modmult(C, I, I);
	nres_modsub(I, P3X, I);
	//Big I=((P1.X+P1.Y)*(P2.x+P2.y)-C-D)%p;//(X_1+Y_1)(X_2+Y_2)-C-D


	nres_modadd(E, B, C);
	nres_modmult(C, H, P3.X);
	//P3.X=((E+B)*H)%p;

	nres_modsub(E, B, C);
	nres_modmult(C, I, P3.Y);
	//P3.Y=((E-B)*I)%p;

	nres_modmult(P1.Z, H, C);
	nres_modmult(C, I, P3.Z);
	//P3.Z=(A*H*I)%p;
	return 0;
}


int InvEdwardsAddition(projectPoint P1, projectPoint P2, projectPoint &P3, big d,big a)
{
	//9M+1S
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big A = mirvar(0), D = mirvar(0), F = mirvar(0);

	//Big A=P1.Z;//Z_1Z_2,Z_2=1
	nres_modmult(P1.Z, P2.Z, A);
	nres_modmult(A, A, B);
	nres_modmult(B, d, B);
	//Big B=(d*A*A)%p;//dA^2,d=1


	nres_modmult(P1.X, P2.X, C);
	//Big C=(P1.X*P2.X)%p;//X_1X_2

	nres_modmult(P1.Y, P2.Y, D);
	//Big D=(P1.Y*P2.Y)%p;//Y_1Y_2


	nres_modmult(C, D, E);
	//Big E=(C*D)%p;//C*D

	nres_modmult(a, D, H);
	nres_modsub(C, H, H);
	//Big H=C-aD;//C-D

	nres_modadd(P1.X, P1.Y, F);
	nres_modadd(P2.X, P2.Y, I);
	nres_modmult(F, I, I);
	nres_modsub(I, C, I);
	nres_modsub(I, D, I);
	//Big I=((P1.X+P1.Y)*(P2.X+P2.Y)-C-D)%p;//(X_1+Y_1)(X_2+Y_2)-C-D


	nres_modadd(E, B, C);
	nres_modmult(C, H, P3.X);
	//P3.X=((E+B)*H)%p;

	nres_modsub(E, B, C);
	nres_modmult(C, I, P3.Y);
	//P3.Y=((E-B)*I)%p;

	nres_modmult(A, H, C);
	nres_modmult(C, I, P3.Z);
	//P3.Z=(A*H*I)%p;
	return 0;
}


int InvEdwardsAdditionSubtraction(projectPoint P1, projectPoint P2, projectPoint &P3,projectPoint &P4, big d, big a)
{
	//9M+1S(+5M)
	//P3=P1+p2,P4=P1-P2
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big A = mirvar(0), D = mirvar(0), F = mirvar(0);
	big C1 = mirvar(0), E1 = mirvar(0), H1 = mirvar(0), I1 = mirvar(0);
	big Zero = mirvar(0);

	//Big A=P1.Z;//Z_1Z_2,Z_2=1
	nres_modmult(P1.Z, P2.Z, A);
	nres_modmult(A, A, B);
	nres_modmult(B, d, B);
	//Big B=(d*A*A)%p;//dA^2,d=1


	nres_modmult(P1.X, P2.X, C);
	//Big C=(P1.X*P2.X)%p;//X_1X_2
	nres_modsub(Zero, C, C1);
	//C1=-C

	nres_modmult(P1.Y, P2.Y, D);
	//Big D=(P1.Y*P2.Y)%p;//Y_1Y_2

	nres_modmult(C, D, E);
	//Big E=(C*D)%p;//C*D
	nres_modsub(Zero, E, E1);
	//E1=-E

	nres_modmult(a, D, H);
	nres_modsub(C1, H, H1);
	//Big H1=C1-aD;//C-D
	nres_modsub(C, H, H);
	//Big H=C-aD;//C-D

	nres_modadd(P1.X, P1.Y, F);
	nres_modadd(P2.X, P2.Y, I);
	nres_modmult(F, I, I);
	nres_modsub(I, C, I);
	nres_modsub(I, D, I);
	//Big I=((P1.X+P1.Y)*(P2.X+P2.Y)-C-D)%p;//(X_1+Y_1)(X_2+Y_2)-C-D
	nres_modsub(P2.Y, P2.X, I1);
	nres_modmult(F, I1, I1);
	nres_modsub(I1, C1, I1);
	nres_modsub(I1, D, I1);
	//Big I1=((P1.X+P1.Y)*(P2.Y-P2.X)-C1-D)%p;//(X_1+Y_1)(X_2+Y_2)-C-D

	nres_modadd(E, B, C);
	nres_modmult(C, H, P3.X);
	//P3.X=((E+B)*H)%p;
	nres_modadd(E1, B, C);
	nres_modmult(C, H1, P4.X);
	//P4.X=((E1+B)*H1)%p;

	nres_modsub(E, B, C);
	nres_modmult(C, I, P3.Y);
	//P3.Y=((E-B)*I)%p;
	nres_modsub(E1, B, C);
	nres_modmult(C, I1, P4.Y);
	//P4.Y=((E1-B)*I1)%p;

	nres_modmult(A, H, C);
	nres_modmult(C, I, P3.Z);
	//P3.Z=(A*H*I)%p;
	nres_modmult(A, H1, C);
	nres_modmult(C, I1, P4.Z);
	//P4.Z=(A*H1*I1)%p;

	return 0;
}


int InvEdwardsNegative(projectPoint P1, projectPoint &P3)
{
	
	nres_negate(P1.X, P3.X);
	copy(P1.Y, P3.Y);
	copy(P1.Z, P3.Z);

	return 0;
}


int InvEdwardsNegative(affinePoint P1, affinePoint &P3)
{

	nres_negate(P1.x, P3.x);
	copy(P1.y, P3.y);
	//copy(P1.Z, P3.Z);

	return 0;
}


int InvEdwardsDouble(projectPoint P1, projectPoint &P3, big d)
{
	big A=mirvar(0),B=mirvar(0),C=mirvar(0),D=mirvar(0),E=mirvar(0);
	nres_modmult(P1.X,P1.X,A);
	//Big A=(P1.X*P1.X)%p;//X_1^2

	nres_modmult(P1.Y,P1.Y,E);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modadd(A,E,C);
	//Big C=A+B;

	nres_modsub(A,E,D);
	//Big D=A-B;

	nres_modmult(P1.X,P1.Y,E);
	nres_modadd(E,E,E);
	//Big E=(2*P1.X*P1.Y)%p;//2X_1Y_1

	nres_modmult(P1.Z,P1.Z,A);
	//A=P1.Z*P1.Z%p;

	nres_modmult(C,D,P3.X);
	//P3.X=(C*D)%p;

	nres_modmult(D,E,P3.Z);
	//P3.Z=(D*E)%p;

	nres_modmult(A,d,D);

	nres_modadd(D,D,D);
	nres_modsub(C,D,D);
	nres_modmult(D,E,P3.Y);
	//P3.Y=((C-2*d*(A))*E)%p;


	return 0;
}

int InvEdwardsTriple(projectPoint P1, projectPoint &P3, big d)
{

	big A=mirvar(0),B=mirvar(0),C=mirvar(0),D=mirvar(0),H=mirvar(0),P=mirvar(0),Q=mirvar(0);

	nres_modmult(P1.X,P1.X,A);
	//Big A=(P1.X*P1.X)%p;//X_1^2


	nres_modmult(P1.Y,P1.Y,B);
	//Big B=(P1.Y*P1.Y)%p;//Y_1^2

	nres_modmult(P1.Z,P1.Z,C);
	//Big C=(P1.Z*P1.Z)%p;//Z_1^2

	nres_modadd(A,B,D);
	//Big D=A+B;


	nres_modmult(d,C,P);
	nres_modsub(D,P,P);
	nres_modadd(P,P,P);
	nres_modadd(P,P,Q);
	//Big E=4*(D-d*C);

	nres_modsub(B,A,H);
	nres_modmult(H,D,H);
	nres_modadd(H,H,H);
	//Big H=(2*D*(B-A))%p;


	nres_modmult(D,D,C);
	//C=(D*D)%p;

	nres_modmult(A,Q,P);
	nres_modsub(C,P,P);
	//Big P=(C-A*E)%p;

	nres_modmult(B,Q,Q);
	nres_modsub(C,Q,Q);
	//Big Q=(C-B*E)%p;



	nres_modadd(H,Q,A);
	nres_modmult(A,Q,B);
	nres_modmult(B,P1.X,P3.X);
	//P3.X=((H+Q)*Q*P1.X)%p;


	nres_modsub(H,P,A);
	nres_modmult(A,P,B);
	nres_modmult(B,P1.Y,P3.Y);
	//P3.Y=((H-P)*P*P1.Y)%p;

	nres_modmult(P,Q,B);
	nres_modmult(B,P1.Z,P3.Z);

	//P3.Z=(P*Q*P1.Z)%p;
	return 0;
}



int DIKMixedAddition(extProjectPoint P1, affinePoint P2, extProjectPoint &P3)
{
	//T=ZZ,.T is .ZZ
	//7M+4S
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big a = mirvar(6), four = mirvar(4);//a=3a=3*2,is curve parameter
	big A = mirvar(0), D = mirvar(0), F = mirvar(0), F4 = mirvar(0), G = mirvar(0), H2 = mirvar(0);
	big E3 = mirvar(0);

	nres_modmult(P2.x, P1.T, A);
	//A=X2*ZZ1
	nres_modmult(P2.y, P1.T, B);
	nres_modmult(B, P1.Z, B);
	//Big B=Y2*ZZ1*Z1
	nres_modsub(P1.X, A, C);
	//C=X1-A
	nres_modsub(P1.Y, B, D);
	nres_modadd(D, D, D);
	//C=2*(Y1-B)
	nres_modmult(C, C, F);
	//F=C^2
	nres_modmult(four, F, F4);
	//F4=4*F
	nres_modadd(P1.Z, C, I);
	nres_modmult(I, I, I);
	nres_modsub(I, P1.T, I);
	nres_modsub(I, F, P3.Z);
	//Z3=(Z1+C)^2-ZZ1-F
	nres_modmult(P3.Z, P3.Z, E);
	//E=Z3^2
	nres_modmult(C, F4, G);
	//G=C*F4
	nres_modmult(A, F4, H);
	//H=A*F4

	nres_modmult(D, D, I);
	nres_modsub(I, G, I);
	nres_modadd(H, H, H2);
	nres_modmult(a,E, E3);
	nres_modsub(I, H2, I);
	nres_modsub(I, E3, P3.X);
	//X3=D^2-G-2H-6E
	
	nres_modsub(H, P3.X, I);
	nres_modmult(D, I, I);
	nres_modmult(B, G, H2);
	nres_modadd(H2, H2, H2);
	nres_modsub(I, H2, P3.Y);
	//Y3=D*(H-X3)-2*B*G

	copy(E, P3.T);
	//ZZ3=E

	return 0;
}


int DIKAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3)
{
	//T=ZZ,.T is .ZZ
	//11M+6S
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big a = mirvar(6), four = mirvar(4), eight = mirvar(8), sixteen = mirvar(16);//a=3a=3*2,is curve parameter
	big A = mirvar(0), D = mirvar(0), F = mirvar(0), F4 = mirvar(0), G = mirvar(0), H2 = mirvar(0);
	big E3 = mirvar(0), ZZZ2 = mirvar(0), ZZZ1 = mirvar(0), CC = mirvar(0), CC8 = mirvar(0), AA = mirvar(0);


	nres_modmult(P2.Z, P2.T, ZZZ2);
	//ZZZ2=Z2*ZZ2
	nres_modmult(P1.Z, P1.T, ZZZ1);
	//ZZZ1=ZZ1*Z1
	nres_modmult(P2.X, P1.T, A);
	nres_modmult(P1.X, P2.T, B);
	nres_modsub(A, B, A);
	//A=X2*ZZ1-X1*ZZ2
	nres_modmult(P2.Y, ZZZ1, B);
	nres_modmult(P1.Y, ZZZ2, D);
	nres_modsub(B, D, B);
	//B=Y2*ZZZ1-Y1*ZZZ2
	nres_modmult(A, A, AA);
	//AA=A^2
	nres_modadd(P2.Z, A, C);
	nres_modmult(C, C, C);
	nres_modsub(C, P2.T, C);
	nres_modsub(C, AA, C);
	//C=(Z2+A)^2-ZZ2-AA
	nres_modmult(C, C, CC);
	//CC=C^2
	nres_modmult(eight, CC, CC8);
	//CC8=8*CC
	nres_modmult(P1.X, CC8, D);
	//D=X1*CC8

	nres_modadd(P1.Z, C, H);
	nres_modmult(H, H, H);
	nres_modsub(H, P1.T, H);
	nres_modsub(H, CC, P3.Z);
	//Z3=(Z1+C)-ZZ1-CC

	nres_modmult(P3.Z, P3.Z, E);
	//E=Z3^2
	nres_modmult(B, B, I);
	nres_modmult(A, AA, G);
	nres_modsub(I, G, I);//I=B^2-A*AA
	nres_modmult(sixteen, I, I);
	nres_modmult(a, E, E3);
	nres_modsub(I, D, I);
	nres_modsub(I, E3, P3.X);
	//X3=16*(B^2-A*AA)-D-6E

	nres_modadd(P3.X, P3.X, H);
	nres_modsub(D, H, I);
	nres_modmult(B, I, I);//I=B*(D-2X3)
	nres_modadd(I, I, I);
	nres_modmult(P1.Y, C, G);
	nres_modmult(G, CC8, G);//G=Y1*C*CC8
	nres_modsub(I, G, P3.Y);
	//Y3=2*B*(D-2X3)-Y1*C*CC8

	copy(E, P3.T);
	//ZZ3=E

	return 0;
}


int DIKAdditionSubtraction(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4)
{
	//T=ZZ,.T is .ZZ
	//11M+6S
	//P3=P1+P2, P4=P1-P2
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big a = mirvar(6), four = mirvar(4), eight = mirvar(8), sixteen = mirvar(16);//a=3a=3*2,is curve parameter
	big A = mirvar(0), D = mirvar(0), F = mirvar(0), F4 = mirvar(0), G = mirvar(0), H2 = mirvar(0);
	big E3 = mirvar(0), ZZZ2 = mirvar(0), ZZZ1 = mirvar(0), CC = mirvar(0), CC8 = mirvar(0), AA = mirvar(0);
	big B1 = mirvar(0),Zero= mirvar(0);


	nres_modmult(P2.Z, P2.T, ZZZ2);
	//ZZZ2=Z2*ZZ2
	nres_modmult(P1.Z, P1.T, ZZZ1);
	//ZZZ1=ZZ1*Z1
	nres_modmult(P2.X, P1.T, A);
	nres_modmult(P1.X, P2.T, B);
	nres_modsub(A, B, A);
	//A=X2*ZZ1-X1*ZZ2
	nres_modmult(P2.Y, ZZZ1, B);
	nres_modsub(Zero, B, B1);
	nres_modmult(P1.Y, ZZZ2, D);
	nres_modsub(B, D, B);
	//B=Y2*ZZZ1-Y1*ZZZ2
	nres_modsub(B1, D, B1);
	//B1=-Y2*ZZZ1-Y1*ZZZ2
	nres_modmult(A, A, AA);
	//AA=A^2
	nres_modadd(P2.Z, A, C);
	nres_modmult(C, C, C);
	nres_modsub(C, P2.T, C);
	nres_modsub(C, AA, C);
	//C=(Z2+A)^2-ZZ2-AA
	nres_modmult(C, C, CC);
	//CC=C^2
	nres_modmult(eight, CC, CC8);
	//CC8=8*CC
	nres_modmult(P1.X, CC8, D);
	//D=X1*CC8

	nres_modadd(P1.Z, C, H);
	nres_modmult(H, H, H);
	nres_modsub(H, P1.T, H);
	nres_modsub(H, CC, P3.Z);
	//Z3=(Z1+C)-ZZ1-CC
	copy(P3.Z,P4.Z);
	//Z4=Z3

	nres_modmult(P3.Z, P3.Z, E);
	//E=Z3^2
	nres_modmult(B, B, I);
	nres_modmult(A, AA, G);//G=A*AA
	nres_modsub(I, G, I);//I=B^2-A*AA
	nres_modmult(sixteen, I, I);
	nres_modmult(a, E, E3);
	nres_modsub(I, D, I);
	nres_modsub(I, E3, P3.X);
	//X3=16*(B^2-A*AA)-D-6E
	nres_modmult(B1, B1, I);
	nres_modsub(I, G, I);//I=B1^2-A*AA
	nres_modmult(sixteen, I, I);
	nres_modsub(I, D, I);
	nres_modsub(I, E3, P4.X);
	////X4=16*(B1^2-A*AA)-D-6E

	nres_modadd(P3.X, P3.X, H);
	nres_modsub(D, H, I);
	nres_modmult(B, I, I);//I=B*(D-2X3)
	nres_modadd(I, I, I);
	nres_modmult(P1.Y, C, G);
	nres_modmult(G, CC8, G);//G=Y1*C*CC8
	nres_modsub(I, G, P3.Y);
	//Y3=2*B*(D-2X3)-Y1*C*CC8
	nres_modadd(P4.X, P4.X, H);
	nres_modsub(D, H, I);
	nres_modmult(B1, I, F4);//F4=B1*(D-2X3)
	nres_modadd(F4, F4, F4);
	nres_modsub(F4, G, P4.Y);
	//Y4=2*B1*(D-2X4)-Y1*C*CC8

	copy(E, P3.T);
	//ZZ3=E
	copy(P3.T, P4.T);
	//T4=T3

	return 0;
}


int DIKDoubling(extProjectPoint P1,extProjectPoint &P3)
{
	//T=ZZ,.T is .ZZ
	//2M+7S
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big six = mirvar(6), three= mirvar(3),four = mirvar(4), eight = mirvar(8), sixteen = mirvar(16);//a=3a=3*2,is curve parameter
	big A = mirvar(0), D = mirvar(0), F = mirvar(0), F2 = mirvar(0), G = mirvar(0), H2 = mirvar(0);
	big E3 = mirvar(0), ZZZ2 = mirvar(0), ZZZ1 = mirvar(0), CC = mirvar(0), CC8 = mirvar(0), AA = mirvar(0);
	big B1 = mirvar(0), Zero = mirvar(0);


	nres_modmult(P1.X, P1.X, A);
	//A=X1^2
	nres_modadd(P1.X, P1.T, B);
	nres_modmult(four, P1.T, C);
	nres_modmult(C, B, B);
	//B=a2*ZZ1*(X1+ZZ1)
	nres_modadd(A, B,C);
	nres_modmult(three, C, C);
	//C=3*(A+B)
	nres_modmult(P1.Y, P1.Y, D);
	//D=Y1^2
	nres_modmult(D, D, E);
	//E=D^2

	nres_modadd(P1.Y, P1.Z, H);
	nres_modmult(H, H, H);
	nres_modsub(H, D, H);
	nres_modsub(H, P1.T, P3.Z);
	//Z3=(Y1+Z1)^2-D-ZZ1

	nres_modmult(P3.Z, P3.Z,P3.T);
	//T3=Z3^2

	nres_modadd(P1.X, D, H);
	nres_modmult(H, H, H);//H=(X1+D)^2
	nres_modsub(H, A, I);
	nres_modsub(I, E, I);//I=(X1+D)^2-A-E
	nres_modadd(I, I, F);
	//F=2*((X1+D)^2-A-E)

	nres_modmult(C, C, H);
	nres_modmult(six, P3.T, I);//I=a3*ZZ3
	nres_modadd(F, F, F2);//F2=2F
	nres_modsub(H, I, H);
	nres_modsub(H, F2, P3.X);
	//X3=C^2-a3*ZZ3-2*F

	nres_modsub(F, P3.X, H);
	nres_modmult(C, H, H);//H=C*(F-X3)
	nres_modmult(eight, E, I);//I=8*E
	nres_modsub(H, I, P3.Y);
	//Y3=C*(F-X3)-8E

	return 0;
}


int DIKTripling(extProjectPoint P1, extProjectPoint &P3)
{
	//T=ZZ,.T is .ZZ
	//6M+6S
	big C = mirvar(0), B = mirvar(0), E = mirvar(0), H = mirvar(0), I = mirvar(0), P3X = mirvar(0);
	big six = mirvar(6), three = mirvar(3), four = mirvar(4), eight = mirvar(8), sixteen = mirvar(16);//a=3a=3*2,is curve parameter
	big A = mirvar(0), D = mirvar(0), F = mirvar(0), F2 = mirvar(0), G = mirvar(0), H2 = mirvar(0);
	big E3 = mirvar(0), YY1 = mirvar(0), Xt = mirvar(0), Yt = mirvar(0), Zt = mirvar(0), Ca = mirvar(0);
	big B1 = mirvar(0), Zero = mirvar(0), nine = mirvar(9);


	nres_modmult(P1.Y, P1.Y, YY1);
	//YY1=Y1^2
	nres_modmult(three, P1.T, B);
	nres_modadd(P1.X, B, A);//A=X1+3*ZZ1
	nres_modmult(A, A, A);
	//A=(X1+3*ZZ1)^2
	nres_modadd(P1.T,P1.T, B);
	nres_modmult(B, A, B);
	//B=a*ZZ1*A
	nres_modadd(YY1, B, Xt);
	//Xt=YY1+B
	nres_modmult(three, B, H);
	nres_modsub(YY1, H, I);//I=YY1-3*B
	nres_modmult(P1.Y, I, Yt);
	//Yt=Y1*(YY1-3*B)
	nres_modmult(P1.X, P1.Z, Zt);
	//Zt=X1*Z1
	nres_modmult(Zt, Zt, C);
	//C=Zt^2
	nres_modadd(C, C, Ca);
	//Ca=a*C

	nres_modmult(four, Ca, G);
	nres_modmult(nine, C, H);
	nres_modsub(G, H, I);//I=4*Ca-9*C
	nres_modsub(I, Xt, I);
	nres_modmult(I, I, I);//I=(4*Ca-9*C-Xt)^2
	nres_modmult(three, Ca, F);
	nres_modmult(F, I, E);
	//E=3*Ca*(4*Ca-9*C-Xt)^2

	nres_modmult(Yt, Yt, I);
	nres_modsub(I, E, P3.X);
	//X3=Yt^2-E

	nres_modmult(four, E, I);
	nres_modadd(P3.X, I, I);
	nres_modmult(Yt, I, P3.Y);
	//Y3=Yt*(X3+4*E)

	nres_modmult(three, Xt, I);
	nres_modmult(I, Zt, P3.Z);
	//Z3=3*Xt*Zt

	nres_modmult(P3.Z,P3.Z,P3.T);
	//ZZ#=Z3^2

	return 0;
}


int DIKNegative(extProjectPoint P1, extProjectPoint &P3)
{
	//T=ZZ,.T is .ZZ


	copy(P1.X, P3.X);
	copy(P1.T, P3.T);
	copy(P1.Z, P3.Z);
	nres_negate(P1.Y, P3.Y);

	return 0;
}