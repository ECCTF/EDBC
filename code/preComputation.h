#include "pointOperation.h"

#pragma once
int preComputationW(affinePoint P1, projectPoint *P, long long &CPUcycles, int Cmax);

int preComputationE(affinePoint P1, projectPoint *P, long long &CPUcycles, int Cmax, big d, big a);

int preComputationD(affinePoint P1, extProjectPoint *P, long long &CPUcycles, int Cmax);

int preComputationEx(affinePoint P1, extProjectPoint *P, long long &CPUcycles, int Cmax, big d);