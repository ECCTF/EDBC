#include "pointOperation.h"
#include "scalarMultiplication.h"
#pragma once
 

int productTree(big A, int &IsSupersingularity);
int newProductTree(big A, doubleBaseType DBCTerm[], int DBCLength, int &IsSupersingularity);

int pairingMontgomeryLadder(big A, int &IsSupersingularity);