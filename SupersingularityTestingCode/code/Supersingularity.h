#include "pointOperation.h"
#include "scalarMultiplication.h"
#pragma once
 
int ConstructEdwardsPoint(big A, affinePoint &P);
int productTree(big A, int &IsSupersingularity);
int productTree(affinePoint P, big A, int &IsSupersingularity);
int newProductTree(big A, int &IsSupersingularity);
int newProductTree(affinePoint P, big A, int &IsSupersingularity, int &costInM);

int pairingMontgomeryLadder(big A, int &IsSupersingularity);