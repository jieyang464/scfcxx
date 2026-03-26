#include "scf.h"
#include "linalg.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

// Build density matrix D = C_occ · C_occᵀ from coefficient matrix C.
// n_occ is the number of occupied spin orbitals (alpha or beta).
inline T2 BuildDensity(const T2& C, int n_occ) {
    return BuildDensityMatrix(C, n_occ);
}

inline void ResizeSquareAndZero(T2& matrix, Eigen::Index nbf) {
    matrix.resize(nbf, nbf);
    matrix.setZero();
}

inline double TraceProduct(const T2& A, const T2& B) {
    const Eigen::Index nbf = A.dimension(0);
    double trace = 0.0;
    for (Eigen::Index i = 0; i < nbf; ++i) {
        for (Eigen::Index j = 0; j < nbf; ++j) {
            trace += A(i, j) * B(j, i);
        }
    }
    return trace;
}

// Diagonalize Fa/Fb in the metric of S using the helper in linalg.h.
inline void DiagonalizeFockMatrices(SCFResults& scfResults) {
    const T2& S = scfResults.integrals.overlap;

    auto diag_a = DiagonalizeInSMetric(scfResults.fockMatrices.Fa, S);
    scfResults.moCoefficients.Ca = std::move(diag_a.C);
    scfResults.fockEigenvalues.eps_a = std::move(diag_a.eps);

    auto diag_b = DiagonalizeInSMetric(scfResults.fockMatrices.Fb, S);
    scfResults.moCoefficients.Cb = std::move(diag_b.C);
    scfResults.fockEigenvalues.eps_b = std::move(diag_b.eps);
}

// Build (alpha/beta) density matrices from the current MO coefficients.
inline void BuildDensityMatrices(SCFResults& scfResults) {
    const int Na = scfResults.settings.Na;
    const int Nb = scfResults.settings.Nb;

    int n_occ_a = Na;
    int n_occ_b = Nb;

    // Break spin symmetry on the first SCF cycle if requested.
    // (Add one extra alpha electron and remove one beta electron.)
    if (scfResults.settings.break_symmetry && scfResults.iteration == 0 && Nb > 0) {
        n_occ_a = Na + 1;
        n_occ_b = Nb - 1;
    }

    const Eigen::Index nbf = scfResults.integrals.hcore.dimension(0);
    n_occ_a = std::max<Eigen::Index>(0, std::min<Eigen::Index>(n_occ_a, nbf));
    n_occ_b = std::max<Eigen::Index>(0, std::min<Eigen::Index>(n_occ_b, nbf));

    scfResults.densityMatrices.Da = BuildDensity(scfResults.moCoefficients.Ca, n_occ_a);
    scfResults.densityMatrices.Db = BuildDensity(scfResults.moCoefficients.Cb, n_occ_b);
}

} // namespace

void JKEngine(SCFResults& scfResults) {
    const T2& Da = scfResults.densityMatrices.Da;
    const T2& Db = scfResults.densityMatrices.Db;
    const T4& eri = scfResults.integrals.eri;

    const Eigen::array<Eigen::IndexPair<int>, 2> j_pairs = {
        Eigen::IndexPair<int>(2, 0), Eigen::IndexPair<int>(3, 1)};
    const Eigen::array<Eigen::IndexPair<int>, 2> k_pairs = {
        Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(3, 1)};

    scfResults.jkResults.J  = eri.contract(Da + Db, j_pairs);
    scfResults.jkResults.Ka = eri.contract(Da, k_pairs);
    scfResults.jkResults.Kb = eri.contract(Db, k_pairs);
}

void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao, T2& h_mo, T4& eri_mo) {
    const Eigen::Index nbf = C.dimension(0);
    h_mo.resize(nbf, nbf);
    eri_mo.resize(nbf, nbf, nbf, nbf);

    // h_mo[p,q] = sum_{mu,nu} C(mu,p) h_ao(mu,nu) C(nu,q)
    // t(mu,q) = sum_nu h_ao(mu,nu) C(nu,q)
    const Eigen::array<Eigen::IndexPair<int>, 1> h1 = {
        Eigen::IndexPair<int>(1, 0)};
    T2 t = h_ao.contract(C, h1); // (mu,q)

    const Eigen::array<Eigen::IndexPair<int>, 1> h2 = {
        Eigen::IndexPair<int>(0, 0)};
    h_mo = C.contract(t, h2); // (p,q)

    // eri_mo[p,q,r,s] = sum_{mu,nu,lam,sig} C(mu,p) C(nu,q) C(lam,r) C(sig,s) eri_ao[mu,nu,lam,sig]
    // Do the 4-index transform explicitly one AO index at a time.
    // This avoids relying on implicit tensor-axis ordering in chained contractions.
    T4 t1(nbf, nbf, nbf, nbf);  // (p, nu, lam, sig)
    T4 t2(nbf, nbf, nbf, nbf);  // (p, q, lam, sig)
    T4 t3(nbf, nbf, nbf, nbf);  // (p, q, r, sig)
    t1.setZero();
    t2.setZero();
    t3.setZero();
    eri_mo.setZero();

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index nu = 0; nu < nbf; ++nu) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index mu = 0; mu < nbf; ++mu) {
                        value += C(mu, p) * eri_ao(mu, nu, lam, sig);
                    }
                    t1(p, nu, lam, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index nu = 0; nu < nbf; ++nu) {
                        value += C(nu, q) * t1(p, nu, lam, sig);
                    }
                    t2(p, q, lam, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index r = 0; r < nbf; ++r) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                        value += C(lam, r) * t2(p, q, lam, sig);
                    }
                    t3(p, q, r, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index r = 0; r < nbf; ++r) {
                for (Eigen::Index s = 0; s < nbf; ++s) {
                    double value = 0.0;
                    for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                        value += C(sig, s) * t3(p, q, r, sig);
                    }
                    eri_mo(p, q, r, s) = value;
                }
            }
        }
    }
}

void ComputeEnergy(SCFResults& scfResults) {
    const T2& Da = scfResults.densityMatrices.Da;
    const T2& Db = scfResults.densityMatrices.Db;
    const T2& H  = scfResults.integrals.hcore;
    const T2& Fa = scfResults.fockMatrices.Fa;
    const T2& Fb = scfResults.fockMatrices.Fb;

    const double electronic_energy = 0.5 * (TraceProduct(Da, H + Fa) + TraceProduct(Db, H + Fb));
    scfResults.energy = electronic_energy + scfResults.nuclear_repulsion;
}


void SCFLoop(SCFResults& scfResults, IFockBuilder& fock_builder) {
    // Initial guess: diagonalize Hcore in S metric and build densities.
    const T2& Hcore = scfResults.integrals.hcore;
    const T2& S     = scfResults.integrals.overlap;

    auto diag0 = DiagonalizeInSMetric(Hcore, S);
    scfResults.moCoefficients.Ca = diag0.C;
    scfResults.moCoefficients.Cb = diag0.C;
    scfResults.fockEigenvalues.eps_a = diag0.eps;
    scfResults.fockEigenvalues.eps_b = diag0.eps;

    scfResults.iteration = 0;
    BuildDensityMatrices(scfResults);

    scfResults.energy = 0.0;
    double prev_energy = 0.0;

    const int max_iter = 100;
    const double energy_tol = 1e-8;

    for (int iter = 0; iter < max_iter; ++iter) {
        scfResults.iteration = iter;

        JKEngine(scfResults);

        double xc_correction = fock_builder.build_fock_and_energy(scfResults);

        DiagonalizeFockMatrices(scfResults);
        BuildDensityMatrices(scfResults);

        ComputeEnergy(scfResults);
        scfResults.energy += xc_correction;

        // Print densities + energy each cycle.
        const auto& Da = scfResults.densityMatrices.Da;
        const auto& Db = scfResults.densityMatrices.Db;
        const Eigen::Index nbf = Da.dimension(0);
        std::cout << "SCF cycle " << iter << " energy=" << scfResults.energy << "\n";
        std::cout << "Da:\n";
        for (Eigen::Index i = 0; i < nbf; ++i) {
            for (Eigen::Index j = 0; j < nbf; ++j) {
                std::cout << Da(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Db:\n";
        for (Eigen::Index i = 0; i < nbf; ++i) {
            for (Eigen::Index j = 0; j < nbf; ++j) {
                std::cout << Db(i, j) << " ";
            }
            std::cout << "\n";
        }

        if (iter > 0) {
            if (std::abs(scfResults.energy - prev_energy) < energy_tol) {
                break;
            }
        }

        prev_energy = scfResults.energy;
    }
}
