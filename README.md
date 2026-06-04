# Ogden standalone C++ project

This is a small, GitHub friendly C++17 project for testing Ogden hyperelastic mechanics before integrating anything into SimuCell3D.

It includes two independent validation layers.

1. A true 3D incompressible Ogden material test using prescribed principal stretches.
2. A triangular surface membrane test where a reference triangle is compared with a current triangle, local stretch is computed, Ogden energy is evaluated, and nodal forces are estimated from the energy gradient.

There is no SimuCell3D dependency and no SimuCell3D integration code.

## Build

```bash
cmake -S . -B build
cmake --build build -j8
```

## Run tests

```bash
cd build
ctest --output-on-failure
```

## Run examples

From the project root:

```bash
./build/ogden_3d_uniaxial > uniaxial.csv
./build/ogden_surface_triangle
```

or:

```bash
bash scripts/run_examples.sh
```

## What this project tests

### 3D Ogden material

For an incompressible 3D uniaxial stretch test:

```text
lambda1 = lambda
lambda2 = 1 / sqrt(lambda)
lambda3 = 1 / sqrt(lambda)
```

The volume stays constant because:

```text
lambda1 * lambda2 * lambda3 = 1
```

The Ogden energy density is:

```text
W = sum_p (mu_p / alpha_p) * (lambda1^alpha_p + lambda2^alpha_p + lambda3^alpha_p - 3)
```

At rest:

```text
lambda1 = lambda2 = lambda3 = 1
W = 0
stress = 0
```

### Surface triangle membrane

For a triangular membrane patch:

```text
reference triangle -> current triangle
F = Ds * Dm^-1
C = F^T * F
lambda1, lambda2 = sqrt(eigenvalues of C)
lambda3 = 1 / (lambda1 * lambda2)
E_face = reference_area * W(lambda1, lambda2, lambda3)
force = -dE/dx
```

This is the standalone version of the per face Ogden idea. Each triangular face becomes a tiny hyperelastic membrane element.

## File map

```text
include/ogden/Vec3.hpp                 small 3D vector class
include/ogden/OgdenMaterial.hpp        Ogden energy and uniaxial stress
include/ogden/SurfaceTriangle.hpp      triangle deformation, stretches, energy, finite difference force
src/ogden_3d_uniaxial.cpp              writes lambda, energy, stress CSV to stdout
src/ogden_surface_triangle.cpp         prints per triangle stretch, energy, force
scripts/run_examples.sh                builds example output files
tests/test_ogden_material.cpp          basic material checks
tests/test_surface_triangle.cpp        geometry and energy checks
```

## Important limitation

The surface triangle force in this standalone project is computed by finite difference of the energy. That is simple and useful for validation, but it is not efficient enough for a full simulation. In a production SimuCell3D implementation, this should be replaced by an analytic force assembly using the first Piola Kirchhoff stress.
