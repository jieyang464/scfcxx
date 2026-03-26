scfcxx
=====

Simple C++ SCF/quantum-chemistry utilities and examples.

Overview
--------
This repository contains a small C++ project implementing SCF-related code, integral providers, and example tests. It uses Eigen (included under `third_party/eigen`) and has optional hooks for Libint2 and LibXC if available on your system.

Quick start
-----------
Prerequisites:
- A C++17-capable compiler (g++ recommended)
- `make` (optional)
- Optional: `libint2` and `libxc` to enable related features

Build (quick):

Run the provided g++ command used in the workspace (example):

```bash
/usr/bin/g++ -fdiagnostics-color=always -g -std=c++17 \
  src/test.cpp src/scf.cpp src/SzaboHeHIntegral.cpp -o test_bin
```

Or run the existing build task in an editor/IDE that uses the provided task.

Run:

```bash
./test_bin
```

Repository layout (key files)
- `src/` — main sources (`scf.cpp`, `test.cpp`, integral providers)
- `third_party/eigen/` — bundled Eigen headers
- `libint2_basis/` — basis set files used by integral providers
- `build/` — build artifacts (ignored by `.gitignore`)

Notes
- If you want to enable Libint2/LibXC integration, install those libraries on your system and update the build to link against them.
- `install_basis_sets.sh` can help populate basis sets used by the project.

Contribution
------------
See `CONTRIBUTING.md` for guidelines.

License
-------
This project is available under the MIT License. See `LICENSE`.

Contact
-------
Create issues or pull requests on GitHub.
