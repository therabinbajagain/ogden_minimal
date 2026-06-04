#pragma once

#include "ogden/OgdenMaterial.hpp"
#include "ogden/Vec3.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ogden {

struct PrincipalStretches2D {
    double lambda1 = 1.0;
    double lambda2 = 1.0;
};

struct TriangleForces {
    Vec3 f1;
    Vec3 f2;
    Vec3 f3;
};

class SurfaceTriangle {
public:
    SurfaceTriangle(Vec3 ref1, Vec3 ref2, Vec3 ref3)
        : ref1_(ref1), ref2_(ref2), ref3_(ref3) {
        initialize_reference_geometry();
    }

    double reference_area() const {
        return ref_area_;
    }

    std::array<double, 6> deformation_gradient(const Vec3& cur1,
                                               const Vec3& cur2,
                                               const Vec3& cur3) const {
        // F is stored column-major as a 3x2 matrix:
        // [ F00 F01
        //   F10 F11
        //   F20 F21 ]
        const Vec3 d1 = cur2 - cur1;
        const Vec3 d2 = cur3 - cur1;

        return {
            d1.x * Dm_inv_[0] + d2.x * Dm_inv_[2],
            d1.y * Dm_inv_[0] + d2.y * Dm_inv_[2],
            d1.z * Dm_inv_[0] + d2.z * Dm_inv_[2],
            d1.x * Dm_inv_[1] + d2.x * Dm_inv_[3],
            d1.y * Dm_inv_[1] + d2.y * Dm_inv_[3],
            d1.z * Dm_inv_[1] + d2.z * Dm_inv_[3]
        };
    }

    std::array<double, 3> right_cauchy_green(const Vec3& cur1,
                                             const Vec3& cur2,
                                             const Vec3& cur3) const {
        const auto F = deformation_gradient(cur1, cur2, cur3);

        const double C00 = F[0] * F[0] + F[1] * F[1] + F[2] * F[2];
        const double C01 = F[0] * F[3] + F[1] * F[4] + F[2] * F[5];
        const double C11 = F[3] * F[3] + F[4] * F[4] + F[5] * F[5];

        return {C00, C01, C11};
    }

    PrincipalStretches2D principal_stretches(const Vec3& cur1,
                                             const Vec3& cur2,
                                             const Vec3& cur3) const {
        const auto C = right_cauchy_green(cur1, cur2, cur3);
        const double C00 = C[0];
        const double C01 = C[1];
        const double C11 = C[2];

        const double half_trace = 0.5 * (C00 + C11);
        const double half_diff = 0.5 * (C00 - C11);
        const double discriminant = std::sqrt(std::max(0.0, half_diff * half_diff + C01 * C01));

        const double eig1 = std::max(0.0, half_trace + discriminant);
        const double eig2 = std::max(0.0, half_trace - discriminant);

        PrincipalStretches2D stretches;
        stretches.lambda1 = std::sqrt(std::max(eig1, eig2));
        stretches.lambda2 = std::sqrt(std::min(eig1, eig2));

        if (stretches.lambda2 < 1e-12) {
            stretches.lambda2 = 1e-12;
        }
        if (stretches.lambda1 < 1e-12) {
            stretches.lambda1 = 1e-12;
        }

        return stretches;
    }

    double energy(const Vec3& cur1,
                  const Vec3& cur2,
                  const Vec3& cur3,
                  const OgdenMaterial& material) const {
        const auto stretches = principal_stretches(cur1, cur2, cur3);
        const double lambda1 = stretches.lambda1;
        const double lambda2 = stretches.lambda2;
        const double lambda3 = 1.0 / (lambda1 * lambda2);
        return ref_area_ * material.energy_density_3d(lambda1, lambda2, lambda3);
    }

    TriangleForces finite_difference_forces(const Vec3& cur1,
                                            const Vec3& cur2,
                                            const Vec3& cur3,
                                            const OgdenMaterial& material,
                                            double eps = 1e-6) const {
        if (!(eps > 0.0)) {
            throw std::runtime_error("finite difference eps must be positive");
        }

        TriangleForces forces;
        forces.f1 = negative_energy_gradient(cur1, cur2, cur3, material, 0, eps);
        forces.f2 = negative_energy_gradient(cur1, cur2, cur3, material, 1, eps);
        forces.f3 = negative_energy_gradient(cur1, cur2, cur3, material, 2, eps);
        return forces;
    }

private:
    Vec3 ref1_;
    Vec3 ref2_;
    Vec3 ref3_;
    double ref_area_ = 0.0;

    // Row-major 2x2 inverse:
    // [ d00 d01
    //   d10 d11 ]
    std::array<double, 4> Dm_inv_ {0.0, 0.0, 0.0, 0.0};

    void initialize_reference_geometry() {
        const Vec3 e1_3d = ref2_ - ref1_;
        const double edge_length = e1_3d.norm();
        if (edge_length < 1e-12) {
            throw std::runtime_error("Degenerate reference triangle: edge length near zero");
        }

        const Vec3 e1 = e1_3d / edge_length;
        const Vec3 v3 = ref3_ - ref1_;
        const Vec3 normal = cross(e1_3d, v3);
        ref_area_ = 0.5 * normal.norm();
        if (ref_area_ < 1e-18) {
            throw std::runtime_error("Degenerate reference triangle: area near zero");
        }

        const Vec3 e2 = normalized(cross(normalized(normal), e1));

        const double q2x = edge_length;
        const double q2y = 0.0;
        const double q3x = dot(e1, v3);
        const double q3y = dot(e2, v3);

        // Dm columns are q2 - q1 and q3 - q1:
        // [ q2x q3x
        //   q2y q3y ]
        const double det = q2x * q3y - q3x * q2y;
        if (std::abs(det) < 1e-18) {
            throw std::runtime_error("Reference Dm matrix is singular");
        }

        Dm_inv_[0] =  q3y / det;
        Dm_inv_[1] = -q3x / det;
        Dm_inv_[2] = -q2y / det;
        Dm_inv_[3] =  q2x / det;
    }

    Vec3 negative_energy_gradient(const Vec3& cur1,
                                  const Vec3& cur2,
                                  const Vec3& cur3,
                                  const OgdenMaterial& material,
                                  int node_index,
                                  double eps) const {
        Vec3 gradient;
        for (int component = 0; component < 3; ++component) {
            Vec3 p1_plus = cur1;
            Vec3 p2_plus = cur2;
            Vec3 p3_plus = cur3;
            Vec3 p1_minus = cur1;
            Vec3 p2_minus = cur2;
            Vec3 p3_minus = cur3;

            perturb(p1_plus, p2_plus, p3_plus, node_index, component, eps);
            perturb(p1_minus, p2_minus, p3_minus, node_index, component, -eps);

            const double e_plus = energy(p1_plus, p2_plus, p3_plus, material);
            const double e_minus = energy(p1_minus, p2_minus, p3_minus, material);
            const double dE_dx = (e_plus - e_minus) / (2.0 * eps);

            set_component(gradient, component, -dE_dx);
        }
        return gradient;
    }

    static void perturb(Vec3& p1, Vec3& p2, Vec3& p3,
                        int node_index, int component, double amount) {
        if (node_index == 0) {
            add_component(p1, component, amount);
        } else if (node_index == 1) {
            add_component(p2, component, amount);
        } else if (node_index == 2) {
            add_component(p3, component, amount);
        } else {
            throw std::runtime_error("Invalid node index");
        }
    }

    static void add_component(Vec3& v, int component, double amount) {
        if (component == 0) v.x += amount;
        else if (component == 1) v.y += amount;
        else if (component == 2) v.z += amount;
        else throw std::runtime_error("Invalid vector component");
    }

    static void set_component(Vec3& v, int component, double value) {
        if (component == 0) v.x = value;
        else if (component == 1) v.y = value;
        else if (component == 2) v.z = value;
        else throw std::runtime_error("Invalid vector component");
    }
};

} // namespace ogden
