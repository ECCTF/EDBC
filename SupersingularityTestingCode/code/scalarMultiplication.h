#include "doublebasechain.h"
#include "pointOperation.h"
#include "pointOperationAtomicity.h"
#pragma once


int scalarMultiplicationJ(affinePoint P, affinePoint &nP, int *item, int number, int Atomicity = 0);


int scalarMultiplicationJ(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, int Atomicity = 0);

int scalarMultiplicationAd(affinePoint P, affinePoint &nP, int *item, int number, int Atomicity = 0);


int scalarMultiplicationAd(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, int Atomicity = 0);


int windowScalarMultiplicationJ(affinePoint P, affinePoint &nP, int *item, int number, projectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationJ(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, projectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationEx(affinePoint P, affinePoint &nP, int *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationEx(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationInvE(affinePoint P, affinePoint &nP, int *item, int number, projectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationInvE(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, projectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationD(affinePoint P, affinePoint &nP, int *item, int number, extProjectPoint *Pre, int Cmax, int Atomicity = 0);


int windowScalarMultiplicationD(affinePoint P, affinePoint &nP, doubleBaseType *item, int number, 
	extProjectPoint *Pre, int Cmax, int Atomicity = 0);


// Compute X:Z only. The caller owns all inputs and output storage.
// Return -1 for a negative scalar; otherwise return 0. Zero maps to infinity.
int MontgomeryLadder(projectPoint P1, big number, projectPoint &Q, big ConstantA);
