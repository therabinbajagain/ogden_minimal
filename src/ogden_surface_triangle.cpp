#include "ogden/OgdenMaterial.hpp"
#include "ogden/SurfaceTriangle.hpp"
#include "ogden/Vec3.hpp"

#include <iomanip>
#include <iostream>

int main() {
    const ogden::OgdenMaterial material({
        {1.0e-3, 2.0},
        {1.0e-4, 5.0}
    });

    const ogden::Vec3 X1(0.0, 0.0, 0.0);
    const ogden::Vec3 X2(1.0, 0.0, 0.0);
    const ogden::Vec3 X3(0.0, 1.0, 0.0);

    const ogden::SurfaceTriangle tri(X1, X2, X3);

    // Current triangle is stretched by 1.4 in x and 0.9 in y.
    const ogden::Vec3 x1(0.0, 0.0, 0.0);
    const ogden::Vec3 x2(1.4, 0.0, 0.0);
    const ogden::Vec3 x3(0.0, 0.9, 0.0);

    const auto stretches = tri.principal_stretches(x1, x2, x3);
    const double energy = tri.energy(x1, x2, x3, material);
    const auto forces = tri.finite_difference_forces(x1, x2, x3, material);

    std::cout << std::setprecision(10);
    std::cout << "reference_area = " << tri.reference_area() << "\n";
    std::cout << "lambda1 = " << stretches.lambda1 << "\n";
    std::cout << "lambda2 = " << stretches.lambda2 << "\n";
    std::cout << "lambda3 = " << 1.0 / (stretches.lambda1 * stretches.lambda2) << "\n";
    std::cout << "face_energy = " << energy << "\n";
    std::cout << "force_node_1 = " << forces.f1 << "\n";
    std::cout << "force_node_2 = " << forces.f2 << "\n";
    std::cout << "force_node_3 = " << forces.f3 << "\n";

    const ogden::Vec3 total_force = forces.f1 + forces.f2 + forces.f3;
    std::cout << "total_force = " << total_force << "\n";

    return 0;
}
