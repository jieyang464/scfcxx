//naive CI implementation, precompute the full H and diagonalize

#include "linalg.h"
#include "types.h"
 
#include <vector>

// A determinant is an (alpha, beta) occupation pattern over n_orbitals.
// Each string is a length-n_orbitals vector of 0/1 occupancy values.
// This makes the code easy to read and avoids bit-twiddling.
using DetString = std::vector<int>;

struct CIDeterminant {
    DetString alpha;  // alpha[i] == 1 means MO i occupied by alpha
    DetString beta;   // beta[i]  == 1 means MO i occupied by beta
};

// Full-CI configuration space
struct CISpace {
    int n_orbitals;
    int n_alpha;
    int n_beta;

    std::vector<CIDeterminant> dets; // all determinants in the CI space

    // Build the determinant list (all alpha/beta combos with correct occupations)
    void Build(int n_orbitals_, int n_alpha_, int n_beta_);
};

// Full-CI results (ground state energy + wavefunction)

struct CISettings {
    int n_orbitals;
    int n_alpha;
    int n_beta;
};

struct CIResult {
    CISettings settings;
    CISpace space;          // determinant basis
    T2 H;                   // CI Hamiltonian matrix in the determinant basis
    T2 hcore_mo;
    T4 eri_mo;
    T2 cicoeff;             // CI eigenvectors; columns are CI states in determinant basis
    T2 eigenvalues;         // CI eigenvalues stored as a diagonal matrix
};

// Diagonal CI matrix element for determinant I:
//   <I|H|I> = sum_{p in occ_i} h_p,p
//          + 1/2 * sum_{p in occ_i} sum_{q in occ_i}
//              [ (p q | p q) - (p q | q p) ]
// occ_i is a spin-orbital list (alpha/beta encoded as even/odd indices).
double diagonal_contribution(const std::vector<int>& occ_i,
                             const T2& h_mo, const T4& eri_mo);


// Matrix element between determinants differing by a single spin-orbital excitation p->q:
//   <I|H|J> = h_pq + sum_{r in occ_i, r!=p} [ (p r | q r) - (p r | r q) ]
// where integrals are spin-orbital values (h_ij and (ij|kl) in MO basis).
// Phase due to fermionic excitation is applied externally by caller.
double single_contribution(int p, int q, const std::vector<int>& occ_i,
                           const T2& h_mo, const T4& eri_mo);

// Matrix element for same-spin double excitation p,q -> a,b (all alpha or all beta):
//   <I|H|J> = (p q | a b) - (p q | b a)
// (Spin orbital exchange term exists only when same spin on both excitations.)
double double_same_spin_contribution(int p, int q, int a, int b,
                                     const T4& eri_mo);  

// Matrix element for mixed-spin double excitation (one alpha, one beta):
//   <I|H|J> = (p q | a b)
// Exchange term is zero because cross-spin exchange is forbidden in spin-orbital integrals.
double double_diff_spin_contribution(int p, int q, int a, int b,
                                     const T4& eri_mo);
                                     
 void BuildCIHamiltonian(CIResult& ciResult);
//form the alpha-beta string tensor and flattern to define the determinant basis
//[[a b] [a b]] -> [aa, ab, ba, bb]
//build Hamiltonian matrix in the determinant basis using Slater-Condon rules and the MO integrals

 void DiagonalizeCI(CIResult& ciResult);


