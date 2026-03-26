#include "SzaboHeHIntegral.h"
#include <iostream>

Integrals SzaboHeHIntegrals() {
    Integrals out;

    // Szabo/Ostlund HeH+ STO-3G example (R = 1.4632 bohr).
    out.hcore.resize(2, 2);
    out.hcore.setZero();
    out.hcore(0, 0) = -2.652744703;
    out.hcore(0, 1) = -1.347205024;
    out.hcore(1, 0) = -1.347205024;
    out.hcore(1, 1) = -1.731828436;

    // Chemist notation: (mu nu | lambda sigma), 0-based indices.
    out.eri.resize(2, 2, 2, 2);
    out.eri.setZero();
    out.eri(0, 0, 0, 0) = 1.307152;
    out.eri(0, 0, 0, 1) = 0.437279;
    out.eri(0, 0, 1, 0) = 0.437279;
    out.eri(0, 0, 1, 1) = 0.605703;
    out.eri(0, 1, 0, 0) = 0.437279;
    out.eri(0, 1, 0, 1) = 0.177267;
    out.eri(0, 1, 1, 0) = 0.177267;
    out.eri(0, 1, 1, 1) = 0.311795;
    out.eri(1, 0, 0, 0) = 0.437279;
    out.eri(1, 0, 0, 1) = 0.177267;
    out.eri(1, 0, 1, 0) = 0.177267;
    out.eri(1, 0, 1, 1) = 0.311795;
    out.eri(1, 1, 0, 0) = 0.605703;
    out.eri(1, 1, 0, 1) = 0.311795;
    out.eri(1, 1, 1, 0) = 0.311795;
    out.eri(1, 1, 1, 1) = 0.774608;

    out.overlap.resize(2, 2);
    out.overlap.setZero();
    out.overlap(0, 0) = 1.0;
    out.overlap(0, 1) = 0.4508;
    out.overlap(1, 0) = 0.4508;
    out.overlap(1, 1) = 1.0;

    return out;
}

IntegralDerivatives SzaboHeHIntegralDerivatives() {
    IntegralDerivatives out;
    std::cout << "SzaboHeHIntegralDerivatives: returning zero derivatives for testing\n";   
    out.d_hcore.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& m : out.d_hcore) {
        m.resize(2, 2);
        m.setZero();
    }

    out.d_eri.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& t : out.d_eri) {
        t.resize(2, 2, 2, 2);
        t.setZero();
    }

    out.d_overlap.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& m : out.d_overlap) {
        m.resize(2, 2);
        m.setZero();
    }

    return out;
}

Integrals SzaboHeHIntegralProvider::ComputeIntegrals() const {
    return SzaboHeHIntegrals();
}

IntegralDerivatives SzaboHeHIntegralProvider::ComputeFirstDerivatives() const {
    return SzaboHeHIntegralDerivatives();
}
