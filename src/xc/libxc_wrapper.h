// Optional LibXC wrapper for local XC evaluation (LDA/GGA)
#pragma once

struct XCInput {
    double rho;      // electron density at point
    double sigma;    // |grad rho|^2 (for GGA); ignored for LDA
};

struct XCOutput {
    double eps_xc;   // energy density per particle
    double v_rho;    // d(epsilon)/d rho (LDA term)
    double v_sigma;  // d(epsilon)/d sigma (GGA term); zero for LDA
};

struct SpinXCInput {
    double rho_a;    // alpha density at point
    double rho_b;    // beta density at point
};

struct SpinXCOutput {
    double eps_xc;   // energy density per particle
    double v_rho_a;  // d(epsilon)/d rho_a
    double v_rho_b;  // d(epsilon)/d rho_b
};

enum class XCFunctionalKind { LDA, GGA };

// Evaluate XC locally using libxc if enabled at build time; otherwise throws.
// functional_id: integer from LibXC functional list (e.g., XC_LDA_X, etc.)
XCOutput evaluate_xc_point(int functional_id, XCFunctionalKind kind, const XCInput& in);

// Evaluate spin-polarized LDA locally using libxc if enabled at build time.
SpinXCOutput evaluate_spin_lda_point(int functional_id, const SpinXCInput& in);
