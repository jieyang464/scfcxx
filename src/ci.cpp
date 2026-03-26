#include "ci.h"
#include "linalg.h"

#include <algorithm>
#include <unordered_map>

namespace {

// Generate all combinations of `k` indices out of `n` in ascending order.
// This is a simple recursive generator.
void GenerateCombinations(int n, int k, int start, std::vector<int>& current,
                          std::vector<std::vector<int>>& out) {
    if (k == 0) {
        out.push_back(current);
        return;
    }
    for (int i = start; i <= n - k; ++i) {
        current.push_back(i);
        GenerateCombinations(n, k - 1, i + 1, current, out);
        current.pop_back();
    }
}

// Convert a determinant (alpha/beta occupancy) to a sorted list of spin-orbital indices.
// Spin orbitals are numbered: (0,alpha), (0,beta), (1,alpha), (1,beta), ...
std::vector<int> BuildSpinOrbitalList(const CIDeterminant& det) {
    const int n = static_cast<int>(det.alpha.size());
    std::vector<int> out;
    out.reserve(2 * n);
    for (int p = 0; p < n; ++p) {
        if (det.alpha[p]) out.push_back(2 * p + 0);
        if (det.beta[p])  out.push_back(2 * p + 1);
    }
    return out;
}

// Spin-orbital integral <pq|rs> in the MO basis (spatial integrals).
// Uses the convention that spin-orbital index is (2*p + spin), spin=0 (alpha), 1 (beta).
inline double SpinOrbitalERI(int p, int q, int r, int s, const T4& eri_mo) {
    const int spin_p = p & 1;
    const int spin_q = q & 1;
    const int spin_r = r & 1;
    const int spin_s = s & 1;

    // Nonzero only when spins match for (p,r) and (q,s)
    if (spin_p != spin_r) return 0.0;
    if (spin_q != spin_s) return 0.0;

    const int p_sp = p >> 1;
    const int q_sp = q >> 1;
    const int r_sp = r >> 1;
    const int s_sp = s >> 1;

    return eri_mo(p_sp, r_sp, q_sp, s_sp);
}

// One-electron integral in spin-orbital basis: h_pq = h_spatial(p/2,q/2) if spins match.
inline double SpinOrbitalH(int p, int q, const T2& h_mo) {
    if ((p & 1) != (q & 1)) return 0.0;
    return h_mo(p >> 1, q >> 1);
}

// Single excitation and double-excitation contribution helpers are declared in ci.h
// and defined in global scope below (outside anonymous namespace).

// Helper: annihilate occupied spin orbital `p` from `occ`, and compute phase.
// Returns phase (+1 / -1) and writes result state to `out`.
// If p is not occupied, returns 0 and leaves out unchanged.
int ApplyAnnihilation(const std::vector<int>& occ, int p, std::vector<int>& out) {
    auto it = std::find(occ.begin(), occ.end(), p);
    if (it == occ.end()) return 0;
    int idx = static_cast<int>(it - occ.begin());
    out = occ;
    out.erase(out.begin() + idx);
    return (idx % 2 == 0) ? +1 : -1;
}

// Helper: create spin orbital `q` into `occ`, and compute phase.
// Returns phase (+1 / -1) and writes result state to `out`.
// If q is already occupied, returns 0.
int ApplyCreation(const std::vector<int>& occ, int q, std::vector<int>& out) {
    if (std::find(occ.begin(), occ.end(), q) != occ.end()) return 0;
    auto it = std::lower_bound(occ.begin(), occ.end(), q);
    int idx = static_cast<int>(it - occ.begin());
    out = occ;
    out.insert(out.begin() + idx, q);
    return (idx % 2 == 0) ? +1 : -1;
}

// Phase for single excitation p -> q from occ_i to occ_j using fermionic creation/annihilation rules.
int ExcitationPhase(const std::vector<int>& occ_i, int p, int q, const std::vector<int>& occ_j) {
    std::vector<int> state1;
    int phase = ApplyAnnihilation(occ_i, p, state1);
    if (phase == 0) return 0;

    std::vector<int> state2;
    int phase2 = ApplyCreation(state1, q, state2);
    if (phase2 == 0) return 0;

    if (state2 != occ_j) return 0;
    return phase * phase2;
}

// Phase for double excitation p,q -> a,b from occ_i to occ_j.
int DoubleExcitationPhase(const std::vector<int>& occ_i,
                          int p, int q,
                          int a, int b,
                          const std::vector<int>& occ_j) {
    std::vector<int> state1;
    int phase = ApplyAnnihilation(occ_i, p, state1);
    if (phase == 0) return 0;
    std::vector<int> state2;
    int phase2 = ApplyAnnihilation(state1, q, state2);
    if (phase2 == 0) return 0;
    std::vector<int> state3;
    int phase3 = ApplyCreation(state2, a, state3);
    if (phase3 == 0) return 0;
    std::vector<int> state4;
    int phase4 = ApplyCreation(state3, b, state4);
    if (phase4 == 0) return 0;

    // Order of applications in the formula matters. We'll also try swapped creation order
    // if needed, to match the sign convention.
    if (state4 == occ_j) {
        return phase * phase2 * phase3 * phase4;
    }

    // Try alternate order a <-> b creation (symmetric two-electron matrix element handles exchange later).
    std::vector<int> state3b;
    int phase3b = ApplyCreation(state2, b, state3b);
    if (phase3b == 0) return 0;
    std::vector<int> state4b;
    int phase4b = ApplyCreation(state3b, a, state4b);
    if (state4b == occ_j) {
        return phase * phase2 * phase3b * phase4b;
    }

    return 0;
}

// Close anonymous namespace before defining public global functions
} // namespace

// Global definitions corresponding to declarations in ci.h
// They are implemented with internal use of SpinOrbital* helpers above.

double diagonal_contribution(const std::vector<int>& occ_i,
                             const T2& h_mo, const T4& eri_mo) {
    double diagonal = 0.0;
    for (int p : occ_i) {
        diagonal += SpinOrbitalH(p, p, h_mo);
    }

    for (int p : occ_i) {
        for (int q : occ_i) {
            diagonal += 0.5 * (SpinOrbitalERI(p, q, p, q, eri_mo)
                               - SpinOrbitalERI(p, q, q, p, eri_mo));
        }
    }

    return diagonal;
}

double single_contribution(int p, int q, const std::vector<int>& occ_i,
                           const T2& h_mo, const T4& eri_mo) {
    double hij = SpinOrbitalH(p, q, h_mo);
    for (int r : occ_i) {
        if (r == p) continue;
        hij += SpinOrbitalERI(p, r, q, r, eri_mo) - SpinOrbitalERI(p, r, r, q, eri_mo);
    }
    return hij;
}

// Double excitation contribution (same spin p,q -> a,b): (pq|ab) - (pq|ba)
double double_same_spin_contribution(int p, int q, int a, int b,
                                     const T4& eri_mo) {
    return SpinOrbitalERI(p, q, a, b, eri_mo) - SpinOrbitalERI(p, q, b, a, eri_mo);
}

// Double excitation contribution (different spin p,q -> a,b): (pq|ab)
double double_diff_spin_contribution(int p, int q, int a, int b,
                                     const T4& eri_mo) {
    return SpinOrbitalERI(p, q, a, b, eri_mo);
}

void CISpace::Build(int n_orbitals_, int n_alpha_, int n_beta_) {
    n_orbitals = n_orbitals_;
    n_alpha = n_alpha_;
    n_beta = n_beta_;

    dets.clear();

    // Generate alpha and beta occupation combinations.
    std::vector<std::vector<int>> alpha_combs;
    std::vector<std::vector<int>> beta_combs;

    std::vector<int> tmp;
    GenerateCombinations(n_orbitals, n_alpha, 0, tmp, alpha_combs);
    tmp.clear();
    GenerateCombinations(n_orbitals, n_beta, 0, tmp, beta_combs);

    // Build determinants from all alpha/beta pairs.
    for (auto& alpha_occ : alpha_combs) {
        for (auto& beta_occ : beta_combs) {
            CIDeterminant det;
            det.alpha.assign(n_orbitals, 0);
            det.beta.assign(n_orbitals, 0);
            for (int i : alpha_occ) det.alpha[i] = 1;
            for (int i : beta_occ) det.beta[i] = 1;
            dets.push_back(std::move(det));
        }
    }
}

void BuildCIHamiltonian(CIResult& ciResult) {
    const int nb = static_cast<int>(ciResult.space.dets.size());
    ciResult.H.resize(nb, nb);
    ciResult.H.setZero();

    // Precompute spin-orbital occupation lists for each determinant.
    std::vector<std::vector<int>> occ_lists;
    occ_lists.reserve(nb);
    for (const auto& det : ciResult.space.dets) {
        occ_lists.push_back(BuildSpinOrbitalList(det));
    }

    const T2& h_mo = ciResult.hcore_mo;
    const T4& eri_mo = ciResult.eri_mo;

    // Build full CI Hamiltonian via Slater-Condon rules.
    for (int i = 0; i < nb; ++i) {
        const auto& occ_i = occ_lists[i];

        ciResult.H(i, i) = ::diagonal_contribution(occ_i, h_mo, eri_mo);

        // Off-diagonal elements
        for (int j = i + 1; j < nb; ++j) {
            const auto& occ_j = occ_lists[j];

            // Determine how many spin orbitals differ.
            std::vector<int> diff_i;
            std::vector<int> diff_j;

            auto it_i = occ_i.begin();
            auto it_j = occ_j.begin();

            while (it_i != occ_i.end() || it_j != occ_j.end()) {
                if (it_i != occ_i.end() && (it_j == occ_j.end() || *it_i < *it_j)) {
                    diff_i.push_back(*it_i);
                    ++it_i;
                } else if (it_j != occ_j.end() && (it_i == occ_i.end() || *it_j < *it_i)) {
                    diff_j.push_back(*it_j);
                    ++it_j;
                } else {
                    ++it_i;
                    ++it_j;
                }
            }

            if (diff_i.size() != diff_j.size()) continue;
            if (diff_i.size() > 2) continue; // more than double excitation

            double hij = 0.0;
            int phase = 0;

            if (diff_i.empty()) {
                // should already be diagonal
                continue;
            } else if (diff_i.size() == 1) {
                // Single excitation
                int p = diff_i[0];
                int q = diff_j[0];

                phase = ExcitationPhase(occ_i, p, q, occ_j);
                if (phase == 0) continue;

                hij = ::single_contribution(p, q, occ_i, h_mo, eri_mo);
            } else if (diff_i.size() == 2) {
                // Double excitation: i,j -> a,b
                int p = diff_i[0];
                int q = diff_i[1];
                int a = diff_j[0];
                int b = diff_j[1];

                const bool same_spin_double =
                    (p & 1) == (q & 1) && (a & 1) == (b & 1) && ((p & 1) == (a & 1));

                // For mixed-spin doubles, pair removed and added spin orbitals by spin.
                // The determinant differences are sorted by index, which is not necessarily
                // the same as the physically meaningful alpha/beta pairing.
                if (!same_spin_double) {
                    if ((p & 1) != (a & 1) || (q & 1) != (b & 1)) {
                        std::swap(a, b);
                    }
                }

                phase = DoubleExcitationPhase(occ_i, p, q, a, b, occ_j);
                if (phase == 0) continue;

                if (same_spin_double) {
                    hij = ::double_same_spin_contribution(p, q, a, b, eri_mo);
                } else {
                    hij = ::double_diff_spin_contribution(p, q, a, b, eri_mo);
                }
            }

            ciResult.H(i, j) = phase * hij;
            ciResult.H(j, i) = ciResult.H(i, j);
        }
    }
}

void DiagonalizeCI(CIResult& ciResult) {
    const int n = static_cast<int>(ciResult.H.dimension(0));
    T2 S = Eye(n); // identity overlap
    auto diag = DiagonalizeInSMetric(ciResult.H, S);

    ciResult.eigenvalues = std::move(diag.eps);
    ciResult.cicoeff   = std::move(diag.C);
    ciResult.settings.n_orbitals = ciResult.space.n_orbitals;
    ciResult.settings.n_alpha = ciResult.space.n_alpha;
    ciResult.settings.n_beta = ciResult.space.n_beta;
}
