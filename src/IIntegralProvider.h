#pragma once

#include <functional>
#include "types.h"

struct Integrals {
    T4 eri;
    T2 hcore;
    T2 overlap;
};

struct IntegralDerivatives {
    std::vector<T4> d_eri;
    std::vector<T2> d_hcore;
    std::vector<T2> d_overlap;
};

struct IIntegralProvider {
    virtual ~IIntegralProvider() = default;

    virtual Integrals ComputeIntegrals() const = 0;

    // Default: not implemented unless overridden.
    virtual IntegralDerivatives ComputeFirstDerivatives() const {
        throw std::runtime_error("Integral derivatives not implemented by this provider");
    }

    virtual double ComputeNuclearRepulsionEnergy() const = 0;
};