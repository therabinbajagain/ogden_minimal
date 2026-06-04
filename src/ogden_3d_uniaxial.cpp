#include "ogden/OgdenMaterial.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    const ogden::OgdenMaterial material({
        {1.0e-3, 2.0},
        {1.0e-4, 5.0}
    });

    std::cout << "lambda,lambda1,lambda2,lambda3,energy_density,nominal_stress\n";

    for (double lambda = 0.60; lambda <= 2.0001; lambda += 0.05) {
        const double lambda1 = lambda;
        const double lambda2 = 1.0 / std::sqrt(lambda);
        const double lambda3 = lambda2;
        const double W = material.energy_density_3d(lambda1, lambda2, lambda3);
        const double P = material.uniaxial_nominal_stress_incompressible(lambda);

        std::cout << std::fixed << std::setprecision(8)
                  << lambda << ","
                  << lambda1 << ","
                  << lambda2 << ","
                  << lambda3 << ","
                  << W << ","
                  << P << "\n";
    }

    return 0;
}
