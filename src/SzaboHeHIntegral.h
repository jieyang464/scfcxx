#pragma once

#include "IIntegralProvider.h"

class SzaboHeHIntegralProvider final : public IIntegralProvider {
public:
    Integrals ComputeIntegrals() const override;
    IntegralDerivatives ComputeFirstDerivatives() const override;
    double ComputeNuclearRepulsionEnergy() const override {
        return 1.3669;  // from Szabo's book
    }
};
