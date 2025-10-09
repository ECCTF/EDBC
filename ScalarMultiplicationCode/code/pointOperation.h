#include "big.h"
#pragma once

struct affinePoint
{
	big x;
	big y;
};


struct projectPoint
{
	big X;
	big Y;
	big Z;
};

struct extProjectPoint
{//T=XY/Z
	big X;
	big Y;
	big Z;
	big T;
};

#include "pointOperation.h"

int JacobiMixedAddition(projectPoint P1, affinePoint P2, projectPoint &P3);

int JacobiAddition(projectPoint P1, projectPoint P2, projectPoint &P3);

int JacobiAdditionSubtraction(projectPoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4);

int JacobiNegative(projectPoint P1, projectPoint &P3);

int JacobiNegative(affinePoint P1, affinePoint &P3);

int JacobiDouble(projectPoint P1, projectPoint &P3);

int JacobiTriple(projectPoint P1, projectPoint &P3);


int EdwardsAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsMixedAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsBothMixedAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsDouble(projectPoint P1, projectPoint &P3, big d);

int EdwardsMixedDouble(projectPoint P1, projectPoint &P3, big d);

int EdwardsDouble(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsMixedDouble(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsDouble(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedDouble(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsDouble(projectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedDouble(projectPoint P1, extProjectPoint &P3, big d);

int EdwardsTriple(projectPoint P1, projectPoint &P3, big d);

int EdwardsMixedTriple(projectPoint P1, projectPoint &P3, big d);

int EdwardsTriple(extProjectPoint P1, projectPoint &P3, big d);


int EdwardsMixedTriple(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsTriple(projectPoint P1, extProjectPoint &P3, big d);

int EdwardsMixedTriple(projectPoint P1, extProjectPoint &P3, big d);


int EdwardsTriple(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedTriple(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsNegative(affinePoint P1, affinePoint &P3, big d);


int EdwardsNegative(projectPoint P1, projectPoint &P3, big d);

int EdwardsNegative(extProjectPoint P1, extProjectPoint &P3, big d);

int EdwardsAdditionSubtraction(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4, big d);


int InvEdwardsMixedAddition(projectPoint P1, affinePoint P2, projectPoint &P3, big d);

int InvEdwardsAddition(projectPoint P1, projectPoint P2, projectPoint &P3, big d, big a);

int InvEdwardsAdditionSubtraction(projectPoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4, big d, big a);

int InvEdwardsNegative(projectPoint P1, projectPoint &P3);

int InvEdwardsNegative(affinePoint P1, affinePoint &P3);

int InvEdwardsDouble(projectPoint P1, projectPoint &P3, big d);


int InvEdwardsTriple(projectPoint P1, projectPoint &P3, big d);


int DIKMixedAddition(extProjectPoint P1, affinePoint P2, extProjectPoint &P3);


int DIKAddition(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3);


int DIKAdditionSubtraction(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4);


int DIKDoubling(extProjectPoint P1, extProjectPoint &P3);


int DIKTripling(extProjectPoint P1, extProjectPoint &P3);

int DIKNegative(extProjectPoint P1, extProjectPoint &P3);