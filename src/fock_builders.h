#pragma once

#include <memory>
#include "types.h"
#include "xc/vxc_evaluator.h"

struct SCFResults;

// Abstract interface used by the SCF driver to build Fock matrices.
class IFockBuilder {
public:
    virtual ~IFockBuilder() = default;
    // Build Fock matrices inside `scfResults.fockMatrices` and return a
    // scalar energy correction that the driver should add to the SCF energy
    // after the standard 0.5*Tr[D(H+F)] evaluation (used for XC corrections).
    virtual double build_fock_and_energy(SCFResults& scfResults) = 0;
};

// Unrestricted HF builder: Fa = H + J - K_a, Fb = H + J - K_b
class UHFBuilder : public IFockBuilder {
public:
    UHFBuilder() = default;
    double build_fock_and_energy(SCFResults& scfResults) override;
};

// Unrestricted KS builder: builds HF-like Coulomb/Exchange terms and adds
// Vxc from the configured vxc functor. Returns the XC energy correction
// (Exc - 0.5 * Tr[P Vxc]) so the driver can form the correct total energy.
class UKSBuilder : public IFockBuilder {
public:
    UKSBuilder(VxcFunctor vxc_functor = make_stub_vxc_functor(),
               double exact_exchange_fraction = 0.0);
    double build_fock_and_energy(SCFResults& scfResults) override;

private:
    VxcFunctor vxc_functor_;
    double exact_exchange_fraction_;
};