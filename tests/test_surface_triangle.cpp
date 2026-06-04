#include "ogden/OgdenMaterial.hpp"
#include "ogden/SurfaceTriangle.hpp"
#include "ogden/Vec3.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double tol = 1e-8) {
    return std::abs(a - b) < tol;
}

double force_norm_sum(const ogden::TriangleForces& f) {
    return f.f1.norm() + f.f2.norm() + f.f3.norm();
}

} // namespace

int main() {
    const ogden::OgdenMaterial material({
        {1.0e-3, 2.0},
        {1.0e-4, 5.0}
    });

    const ogden::Vec3 X1(0.0, 0.0, 0.0);
    const ogden::Vec3 X2(1.0, 0.0, 0.0);
    const ogden::Vec3 X3(0.0, 1.0, 0.0);
    const ogden::SurfaceTriangle tri(X1, X2, X3);

    assert(close(tri.reference_area(), 0.5));

    // Rest state.
    {
        const auto stretches = tri.principal_stretches(X1, X2, X3);
        const double E = tri.energy(X1, X2, X3, material);
        const auto forces = tri.finite_difference_forces(X1, X2, X3, material);

        assert(close(stretches.lambda1, 1.0));
        assert(close(stretches.lambda2, 1.0));
        assert(close(E, 0.0));
        assert(force_norm_sum(forces) < 1e-7);
    }

    // Pure rotation in 3D should not create strain energy.
    {
        const ogden::Vec3 r1(0.0, 0.0, 0.0);
        const ogden::Vec3 r2(0.0, 1.0, 0.0);
        const ogden::Vec3 r3(-1.0, 0.0, 0.0);
        const auto stretches = tri.principal_stretches(r1, r2, r3);
        const double E = tri.energy(r1, r2, r3, material);

        assert(close(stretches.lambda1, 1.0));
        assert(close(stretches.lambda2, 1.0));
        assert(close(E, 0.0));
    }

    // Known biaxial stretch.
    {
        const ogden::Vec3 x1(0.0, 0.0, 0.0);
        const ogden::Vec3 x2(2.0, 0.0, 0.0);
        const ogden::Vec3 x3(0.0, 1.5, 0.0);
        const auto stretches = tri.principal_stretches(x1, x2, x3);
        const double E = tri.energy(x1, x2, x3, material);
        const auto forces = tri.finite_difference_forces(x1, x2, x3, material);
        const ogden::Vec3 total_force = forces.f1 + forces.f2 + forces.f3;

        assert(close(stretches.lambda1, 2.0));
        assert(close(stretches.lambda2, 1.5));
        assert(E > 0.0);
        assert(total_force.norm() < 1e-6);
    }

    std::cout << "test_surface_triangle passed\n";
    return 0;
}
