#include "big.h"
#pragma once


#include "pointOperation.h"

int JacobiMixedAdditionAtomicity(projectPoint P1, affinePoint P2, projectPoint &P3);

int JacobiAdditionAtomicity(projectPoint P1, projectPoint P2, projectPoint &P3);

int JacobiAdditionSubtractionAtomicity(projectPoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4);

int JacobiNegativeAtomicity(projectPoint P1, projectPoint &P3);

int JacobiNegativeAtomicity(affinePoint P1, affinePoint &P3);

int JacobiDoubleAtomicity(projectPoint P1, projectPoint &P3);

int JacobiTripleAtomicity(projectPoint P1, projectPoint &P3);


int EdwardsAdditionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsMixedAdditionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsBothMixedAdditionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, big d);

int EdwardsDoubleAtomicity(projectPoint P1, projectPoint &P3, big d);

int EdwardsMixedDoubleAtomicity(projectPoint P1, projectPoint &P3, big d);

int EdwardsDoubleAtomicity(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsMixedDoubleAtomicity(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsDoubleAtomicity(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedDoubleAtomicity(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsDoubleAtomicity(projectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedDoubleAtomicity(projectPoint P1, extProjectPoint &P3, big d);

int EdwardsTripleAtomicity(projectPoint P1, projectPoint &P3, big d);

int EdwardsMixedTripleAtomicity(projectPoint P1, projectPoint &P3, big d);

int EdwardsTripleAtomicity(extProjectPoint P1, projectPoint &P3, big d);


int EdwardsMixedTripleAtomicity(extProjectPoint P1, projectPoint &P3, big d);

int EdwardsTripleAtomicity(projectPoint P1, extProjectPoint &P3, big d);

int EdwardsMixedTripleAtomicity(projectPoint P1, extProjectPoint &P3, big d);


int EdwardsTripleAtomicity(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsMixedTripleAtomicity(extProjectPoint P1, extProjectPoint &P3, big d);


int EdwardsNegativeAtomicity(affinePoint P1, affinePoint &P3, big d);


int EdwardsNegativeAtomicity(projectPoint P1, projectPoint &P3, big d);

int EdwardsNegativeAtomicity(extProjectPoint P1, extProjectPoint &P3, big d);

int EdwardsAdditionSubtractionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4, big d);


int InvEdwardsMixedAdditionAtomicity(projectPoint P1, affinePoint P2, projectPoint &P3, big d);

int InvEdwardsAdditionAtomicity(projectPoint P1, projectPoint P2, projectPoint &P3, big d, big a);

int InvEdwardsAdditionSubtractionAtomicity(projectPoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4, big d, big a);

int InvEdwardsNegativeAtomicity(projectPoint P1, projectPoint &P3);

int InvEdwardsNegativeAtomicity(affinePoint P1, affinePoint &P3);

int InvEdwardsDoubleAtomicity(projectPoint P1, projectPoint &P3, big d);


int InvEdwardsTripleAtomicity(projectPoint P1, projectPoint &P3, big d);


int DIKMixedAdditionAtomicity(extProjectPoint P1, affinePoint P2, extProjectPoint &P3);


int DIKAdditionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3);


int DIKAdditionSubtractionAtomicity(extProjectPoint P1, extProjectPoint P2, extProjectPoint &P3, extProjectPoint &P4);


int DIKDoublingAtomicity(extProjectPoint P1, extProjectPoint &P3);


int DIKTriplingAtomicity(extProjectPoint P1, extProjectPoint &P3);

int DIKNegativeAtomicity(extProjectPoint P1, extProjectPoint &P3);