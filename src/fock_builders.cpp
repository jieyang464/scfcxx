#include "fock_builders.h"
#include "scf.h"
#include "xc/vxc_evaluator.h"

#include <algorithm>

// UHF builder implementation.
// Fock matrices are constructed as:
//   F_a = H + J - K_a
//   F_b = H + J - K_b
// where J_{pq} = \nabla_{\mu,\nu,\lambda,\sigma} (D_a + D_b)_{\nu\mu} (eri_{\mu\nu\lambda\sigma})
// and K_{pq} acts only on like-spin density (here provided already as Ka/Kb).
// In our code, `scfResults.jkResults` is populated by `JKEngine` and contains
// J, Ka, Kb. We therefore assemble F by simple matrix arithmetic.
double UHFBuilder::build_fock_and_energy(SCFResults& scfResults) {
    const T2& H = scfResults.integrals.hcore;
    const T2& J = scfResults.jkResults.J;
    const T2& Ka = scfResults.jkResults.Ka;
    const T2& Kb = scfResults.jkResults.Kb;

    // Fa = H + J - Ka
    scfResults.fockMatrices.Fa = H + J - Ka;
    // Fb = H + J - Kb
    scfResults.fockMatrices.Fb = H + J - Kb;

    // HF-only builder returns zero extra energy; the driver computes the
    // canonical HF energy via 0.5 * Tr[D (H + F)].
    return 0.0;
}

UKSBuilder::UKSBuilder(VxcFunctor vxc_functor, double exact_exchange_fraction)
    : vxc_functor_(std::move(vxc_functor)),
      exact_exchange_fraction_(exact_exchange_fraction) {}


// UKS builder implementation.
// Build Fock matrices similarly to UHF but also add the XC potential Vxc
// evaluated by the configured vxc functor. The driver will compute
// E_scf = 0.5 * Tr[D (H + F)] and we return the correction
//   E_xc_corr = Exc - 0.5 * Tr[P Vxc]
// so that the final total energy becomes the usual DFT energy
//   E = 0.5 Tr[D(H+F)] + (Exc - 0.5 Tr[P Vxc]) = Tr[DH] + 0.5 Tr[D J] + Exc
// (the trace `Tr[P Vxc]` is approximated by the grid routine and returned
// in `tr_PVxc`).
double UKSBuilder::build_fock_and_energy(SCFResults& scfResults) {
    const T2& H = scfResults.integrals.hcore;
    const T2& J = scfResults.jkResults.J;
    const T2& Ka = scfResults.jkResults.Ka;
    const T2& Kb = scfResults.jkResults.Kb;

    // Start with HF-like terms, allowing for fractional exact exchange
    // if requested in settings: F = H + J - ax * K.
    const double ax = exact_exchange_fraction_;

    scfResults.fockMatrices.Fa = H + J - ax * Ka;
    scfResults.fockMatrices.Fb = H + J - ax * Kb;

    // Prepare Vxc buffers (initialized to zero by contract).
    const Eigen::Index nbf = H.dimension(0);
    scfResults.vxcResults.Vxc_a.resize(nbf, nbf);
    scfResults.vxcResults.Vxc_b.resize(nbf, nbf);
    scfResults.vxcResults.Vxc_a.setZero();
    scfResults.vxcResults.Vxc_b.setZero();
    scfResults.vxcResults.Exc = 0.0;
    scfResults.vxcResults.tr_PVxc = 0.0;

    UksDensityInput din{ scfResults.densityMatrices.Da,
                          scfResults.densityMatrices.Db,
                          scfResults };
    UksVxcOutput dout(scfResults.vxcResults.Vxc_a,
                      scfResults.vxcResults.Vxc_b,
                      scfResults.vxcResults.Exc,
                      scfResults.vxcResults.tr_PVxc);

    // Call the configured Vxc functor (may be a libxc-backed evaluator or a stub).
    vxc_functor_(din, dout);

    // Add Vxc to Fock matrices: F <- F + Vxc
    scfResults.fockMatrices.Fa += scfResults.vxcResults.Vxc_a;
    scfResults.fockMatrices.Fb += scfResults.vxcResults.Vxc_b;

    // Return the XC energy correction: Exc - 0.5 * Tr[P Vxc]. The driver will
    // add this to the computed 0.5*Tr[D(H+F)].
    return scfResults.vxcResults.Exc - 0.5 * scfResults.vxcResults.tr_PVxc;
}
