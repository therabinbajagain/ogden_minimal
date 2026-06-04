#include "ogden/OgdenMaterial.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

} // namespace

int main() {
    const ogden::OgdenMaterial material({
        {1.0e-3, 2.0},
        {1.0e-4, 5.0}
    });

    const double W_rest = material.energy_density_3d(1.0, 1.0, 1.0);
    const double P_rest = material.uniaxial_nominal_stress_incompressible(1.0);

    assert(close(W_rest, 0.0));
    assert(close(P_rest, 0.0));

    const double lambda = 1.5;
    const double l2 = 1.0 / std::sqrt(lambda);
    const double W_stretched = material.energy_density_3d(lambda, l2, l2);
    const double P_stretched = material.uniaxial_nominal_stress_incompressible(lambda);

    assert(W_stretched > 0.0);
    assert(P_stretched > 0.0);

    const double lambda_compressed = 0.8;
    const double l2c = 1.0 / std::sqrt(lambda_compressed);
    const double W_compressed = material.energy_density_3d(lambda_compressed, l2c, l2c);

    assert(W_compressed > 0.0);

    std::cout << "test_ogden_material passed\n";
    return 0;
}
