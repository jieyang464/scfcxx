// scf.h
#pragma once
#include "IIntegralProvider.h"
#include "types.h"
#include "fock_builders.h"

struct SCFSettings {
    bool break_symmetry {true};  //add one alpha e and remove one beta e in the first iteration 
    int Na {0};
    int Nb {0};

    // DFT / UKS settings
    bool use_dft {false};
    int xc_functional_id {0};
    double exact_exchange_fraction {0.0};
};
//D, K and Vxc are always allocated to simplify logic 
struct JKResults { T2 J; T2 Ka; T2 Kb; };
struct VXCResults { T2 Vxc_a; T2 Vxc_b; double Exc{0.0}; double tr_PVxc{0.0}; };
struct DensityMatrices { T2 Da; T2 Db; };
struct FockMatrices { T2 Fa; T2 Fb; };
struct MOCoefficients { T2 Ca; T2 Cb; };
struct FockEigenvalues { T2 eps_a; T2 eps_b; };

struct SCFResults {
    SCFSettings     settings;
    Integrals       integrals; 
    JKResults       jkResults;
    VXCResults      vxcResults;
    DensityMatrices densityMatrices;
    FockMatrices    fockMatrices;
    MOCoefficients  moCoefficients;
    FockEigenvalues fockEigenvalues;
    int             iteration{0};
    double          energy{0.0};
    double nuclear_repulsion{0.0};
};

void SCFLoop(SCFResults& scfResults, IFockBuilder& fock_builder);

// Utilities
void JKEngine(SCFResults& scfResults); 
void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao,
                    T2& h_mo, T4& eri_mo);