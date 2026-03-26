#include "Libint2IntegralProvider.h"

#if SCFCXX_ENABLE_LIBINT2_PROVIDER

#include <libint2.hpp>

#include <algorithm>
#include <stdexcept>

struct Libint2IntegralProvider::Impl {
    Molecule molecule;
    std::vector<std::string> basis_by_atom;
    IntegralBuildOptions options;

    libint2::BasisSet basis;
    std::vector<std::pair<double, std::array<double, 3>>> nuclear_charges;
    std::vector<long> shell2atom;  // per-shell mapping to molecule atom index

    Impl(const Molecule& molecule_, 
        std::vector<std::string> basis_by_atom_,
         IntegralBuildOptions options_)
        : molecule(molecule_),
          basis_by_atom(std::move(basis_by_atom_)),
          options(options_) {
        libint2::initialize();

        // Build the libint2 atom list.
        std::vector<libint2::Atom> atoms;
        atoms.reserve(molecule.atoms.size());
        nuclear_charges.reserve(molecule.atoms.size());
        for (auto const& a : molecule.atoms) {
            const auto& libint_atom = a.atom;
            atoms.push_back(libint_atom);

            nuclear_charges.push_back({
                static_cast<double>(libint_atom.atomic_number),
                {libint_atom.x, libint_atom.y, libint_atom.z}});
        }

        // Determine a basis-set name. If per-atom basis is provided,
        // we require all entries are identical (for now).
        std::string basis_name = "";
        if (!this->basis_by_atom.empty()) {
            basis_name = this->basis_by_atom.front();
            if (!std::all_of(this->basis_by_atom.begin(), this->basis_by_atom.end(),
                             [&](const std::string& b) { return b == basis_name; })) {
                throw std::invalid_argument(
                    "Libint2IntegralProvider currently only supports a single basis set name for all atoms.");
            }
        }
        if (basis_name.empty()) {
            basis_name = "sto-3g";
        }

        basis = libint2::BasisSet(basis_name, atoms, true);
        shell2atom = basis.shell2atom(atoms);
    }
};

Libint2IntegralProvider::Libint2IntegralProvider(
    const Molecule& molecule, std::vector<std::string> basis_by_atom,
    IntegralBuildOptions options)
    : impl_(std::make_unique<Impl>(molecule, std::move(basis_by_atom),
                                   options)) {}

Libint2IntegralProvider::~Libint2IntegralProvider() = default;

Libint2IntegralProvider::Libint2IntegralProvider(const Libint2IntegralProvider& other)
    : impl_(std::make_unique<Impl>(*other.impl_)) {}

Libint2IntegralProvider& Libint2IntegralProvider::operator=(
    const Libint2IntegralProvider& other) {
    impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}

Libint2IntegralProvider::Libint2IntegralProvider(Libint2IntegralProvider&&) noexcept = default;
Libint2IntegralProvider& Libint2IntegralProvider::operator=(
    Libint2IntegralProvider&&) noexcept = default;

void Libint2IntegralProvider::SetMolecule(const Molecule& molecule) {
    impl_->molecule = molecule;
    // Rebuild basis and nuclear charges to keep in sync.
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

void Libint2IntegralProvider::SetBasisByAtom(std::vector<std::string> basis_by_atom) {
    impl_->basis_by_atom = std::move(basis_by_atom);
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

void Libint2IntegralProvider::SetOptions(IntegralBuildOptions options) {
    impl_->options = options;
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

static void FillSymmetricMatrix(Eigen::Tensor<double, 2>& mat,
                                size_t i0,
                                size_t j0,
                                int n1,
                                int n2,
                                const double* buf) {
    for (int p = 0; p < n1; ++p) {
        for (int q = 0; q < n2; ++q) {
            const double val = buf[p * n2 + q];
            mat(i0 + p, j0 + q) = val;
            if (i0 + p != j0 + q) {
                mat(j0 + q, i0 + p) = val;
            }
        }
    }
}

Integrals Libint2IntegralProvider::ComputeIntegrals() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());

    // One-body integrals
    const auto max_nprim = basis.max_nprim();
    const auto max_l = basis.max_l();

    const int max_l_i = static_cast<int>(max_l);
    // Compute only the 0th-order (non-derivative) integrals here.
    libint2::Engine overlap_engine(libint2::Operator::overlap, max_nprim,
                                   max_l_i, 0);
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, max_nprim,
                                   max_l_i, 0);
    libint2::Engine nuclear_engine(
        libint2::Operator::nuclear, max_nprim, max_l_i, 0,
        std::numeric_limits<double>::epsilon(), impl_->nuclear_charges);

    Integrals out;
    out.overlap = T2(nbf, nbf);
    out.overlap.setZero();
    out.hcore = T2(nbf, nbf);
    out.hcore.setZero();
    out.eri = T4(nbf, nbf, nbf, nbf);
    out.eri.setZero();

    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    // Build all one-body integrals.
    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];

            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            const int n1 = static_cast<int>(shell1.size());
            const int n2 = static_cast<int>(shell2.size());

            // Overlap
            auto ov_results = overlap_engine.compute1(shell1, shell2);
            if (ov_results[0]) {
                FillSymmetricMatrix(out.overlap, i0, j0, n1, n2, ov_results[0]);
            }

            // Kinetic
            auto kin_results = kinetic_engine.compute1(shell1, shell2);
            if (kin_results[0]) {
                FillSymmetricMatrix(out.hcore, i0, j0, n1, n2, kin_results[0]);
            }

            // Nuclear attraction
            auto nuc_results = nuclear_engine.compute1(shell1, shell2);
            if (nuc_results[0]) {
                // Add nuclear attraction to hcore (hcore = T + V)
                for (int p = 0; p < n1; ++p) {
                    for (int q = 0; q < n2; ++q) {
                        const double val = nuc_results[0][p * n2 + q];
                        out.hcore(i0 + p, j0 + q) += val;
                        if (i0 + p != j0 + q) {
                            out.hcore(j0 + q, i0 + p) += val;
                        }
                    }
                }
            }
        }
    }

    // Two-electron integrals (ERI) using chemist notation (mu nu | lambda sigma)
    libint2::Engine eri_engine(libint2::Operator::coulomb, max_nprim, max_l_i, 0);

    for (size_t a = 0; a < nshell; ++a) {
        for (size_t b = 0; b < nshell; ++b) {
            for (size_t c = 0; c < nshell; ++c) {
                for (size_t d = 0; d < nshell; ++d) {
                    const auto& shell_a = basis[a];
                    const auto& shell_b = basis[b];
                    const auto& shell_c = basis[c];
                    const auto& shell_d = basis[d];

                    // Use template deriv_order matching the engine's runtime deriv order.
                    // libint2 asserts that these must match.
                    // Always compute non-derivative ERIs here.
                    const auto& eri_results =
                        eri_engine.compute2<libint2::Operator::coulomb,
                                            libint2::BraKet::xx_xx,
                                            0>(shell_a, shell_b, shell_c,
                                              shell_d);
                    if (eri_results.empty() || !eri_results[0]) continue;

                    const size_t a0 = shell2bf[a];
                    const size_t b0 = shell2bf[b];
                    const size_t c0 = shell2bf[c];
                    const size_t d0 = shell2bf[d];
                    const int na = static_cast<int>(shell_a.size());
                    const int nb = static_cast<int>(shell_b.size());
                    const int nc = static_cast<int>(shell_c.size());
                    const int nd = static_cast<int>(shell_d.size());
                    const double* buf = eri_results[0];

                    for (int ia = 0; ia < na; ++ia) {
                        for (int ib = 0; ib < nb; ++ib) {
                            for (int ic = 0; ic < nc; ++ic) {
                                for (int id = 0; id < nd; ++id) {
                                    const size_t idx =
                                        (((static_cast<size_t>(ia) * nb + ib) * nc + ic) * nd) +
                                        id;
                                    out.eri(a0 + ia, b0 + ib, c0 + ic, d0 + id) =
                                        buf[idx];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

double Libint2IntegralProvider::ComputeNuclearRepulsionEnergy() const {
    const auto& nuc = impl_->nuclear_charges;
    const int natoms = static_cast<int>(nuc.size());
    double energy = 0.0;
    for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
            const double qi = nuc[i].first;
            const double qj = nuc[j].first;
            const auto& ri = nuc[i].second;
            const auto& rj = nuc[j].second;
            const double dx = ri[0] - rj[0];
            const double dy = ri[1] - rj[1];
            const double dz = ri[2] - rj[2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r > 0.0) {
                energy += qi * qj / r;
            }
        }
    }
    return energy;
}

IntegralDerivatives Libint2IntegralProvider::ComputeFirstDerivatives() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const int natoms = static_cast<int>(impl_->molecule.atoms.size());

    const int ncoords = 3 * natoms;
    IntegralDerivatives out;
    out.d_overlap.assign(ncoords, T2(nbf, nbf));
    out.d_hcore.assign(ncoords, T2(nbf, nbf));
    out.d_eri.assign(ncoords, T4(nbf, nbf, nbf, nbf));

    for (auto& m : out.d_overlap) m.setZero();
    for (auto& m : out.d_hcore) m.setZero();
    for (auto& m : out.d_eri) m.setZero();

    const auto max_nprim = basis.max_nprim();
    const auto max_l = basis.max_l();
    const int max_l_i = static_cast<int>(max_l);

    libint2::Engine overlap_engine(libint2::Operator::overlap, max_nprim,
                                   max_l_i, 1);
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, max_nprim,
                                   max_l_i, 1);
    libint2::Engine nuclear_engine(
        libint2::Operator::nuclear, max_nprim, max_l_i, 1,
        std::numeric_limits<double>::epsilon(), impl_->nuclear_charges);
    libint2::Engine eri_engine(libint2::Operator::coulomb, max_nprim, max_l_i, 1);

    const auto& shell2bf = basis.shell2bf();
    const auto& shell2atom = impl_->shell2atom;
    const auto nshell = basis.size();

    auto process_derivs = [&](const double* buf, int n1, int n2, size_t i0,
                              size_t j0, int atom_index, int coord, auto& target) {
        if (atom_index < 0 || atom_index >= natoms) return;
        auto& mat = target[atom_index * 3 + coord];
        for (int p = 0; p < n1; ++p) {
            for (int q = 0; q < n2; ++q) {
                const double val = buf[p * n2 + q];
                mat(i0 + p, j0 + q) += val;
                if (i0 + p != j0 + q) {
                    mat(j0 + q, i0 + p) += val;
                }
            }
        }
    };

    // One-body integral derivatives.
    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            const int n1 = static_cast<int>(shell1.size());
            const int n2 = static_cast<int>(shell2.size());

            const int atom1 = static_cast<int>(shell2atom[s1]);
            const int atom2 = static_cast<int>(shell2atom[s2]);

            auto handle_one_body = [&](libint2::Engine& engine,
                                       std::vector<T2>& target) {
                auto results = engine.compute1(shell1, shell2);
                const auto nsets = results.size();
                for (size_t idx = 1; idx < nsets; ++idx) {
                    const auto deriv = idx - 1;  // 0-based derivative index
                    const int center = static_cast<int>(deriv / 3);
                    const int coord = static_cast<int>(deriv % 3);
                    const int atom = (center == 0) ? atom1 : atom2;
                    if (results[idx]) {
                        process_derivs(results[idx], n1, n2, i0, j0, atom, coord,
                                       target);
                    }
                }
            };

            handle_one_body(overlap_engine, out.d_overlap);
            handle_one_body(kinetic_engine, out.d_hcore);
            handle_one_body(nuclear_engine, out.d_hcore);
        }
    }

    // Two-electron integral derivatives.
    for (size_t a = 0; a < nshell; ++a) {
        for (size_t b = 0; b < nshell; ++b) {
            for (size_t c = 0; c < nshell; ++c) {
                for (size_t d = 0; d < nshell; ++d) {
                    const auto& shell_a = basis[a];
                    const auto& shell_b = basis[b];
                    const auto& shell_c = basis[c];
                    const auto& shell_d = basis[d];

                    const size_t a0 = shell2bf[a];
                    const size_t b0 = shell2bf[b];
                    const size_t c0 = shell2bf[c];
                    const size_t d0 = shell2bf[d];
                    const int na = static_cast<int>(shell_a.size());
                    const int nb = static_cast<int>(shell_b.size());
                    const int nc = static_cast<int>(shell_c.size());
                    const int nd = static_cast<int>(shell_d.size());

                    const int atom_a = static_cast<int>(shell2atom[a]);
                    const int atom_b = static_cast<int>(shell2atom[b]);
                    const int atom_c = static_cast<int>(shell2atom[c]);
                    const int atom_d = static_cast<int>(shell2atom[d]);

                    auto eri_results =
                        eri_engine.compute2<libint2::Operator::coulomb,
                                            libint2::BraKet::xx_xx,
                                            1>(shell_a, shell_b, shell_c,
                                              shell_d);
                    const auto nsets = eri_results.size();
                    if (nsets == 0) continue;

                    const size_t expected_derivs = 3 * 4;  // 4 centers (a,b,c,d), each has x/y/z
                    const size_t offset =
                        (nsets == expected_derivs + 1 && eri_results[0]) ? 1 : 0;

                    for (size_t idx = offset; idx < nsets; ++idx) {
                        const auto deriv = idx - offset;
                        const int center = static_cast<int>(deriv / 3);
                        const int coord = static_cast<int>(deriv % 3);
                        int atom = -1;
                        switch (center) {
                        case 0:
                            atom = atom_a;
                            break;
                        case 1:
                            atom = atom_b;
                            break;
                        case 2:
                            atom = atom_c;
                            break;
                        case 3:
                            atom = atom_d;
                            break;
                        default:
                            break;
                        }
                        if (atom < 0 || atom >= natoms) continue;

                        const double* buf = eri_results[idx];
                        if (!buf) continue;
                        auto& tensor = out.d_eri[atom * 3 + coord];

                        for (int ia = 0; ia < na; ++ia) {
                            for (int ib = 0; ib < nb; ++ib) {
                                for (int ic = 0; ic < nc; ++ic) {
                                    for (int id = 0; id < nd; ++id) {
                                        const size_t idxx = (((static_cast<size_t>(ia) * nb + ib) * nc + ic) * nd) + id;
                                        tensor(a0 + ia, b0 + ib, c0 + ic, d0 + id) +=
                                            buf[idxx];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

#endif
