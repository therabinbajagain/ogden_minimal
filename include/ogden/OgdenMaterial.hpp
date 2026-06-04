#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace ogden {

struct OgdenTerm {
    double mu = 0.0;
    double alpha = 0.0;
};

class OgdenMaterial {
public:
    explicit OgdenMaterial(std::vector<OgdenTerm> terms)
        : terms_(std::move(terms)) {
        if (terms_.empty()) {
            throw std::runtime_error("OgdenMaterial needs at least one term");
        }
        for (const auto& term : terms_) {
            if (std::abs(term.alpha) < 1e-12) {
                throw std::runtime_error("Ogden alpha must be non-zero");
            }
        }
    }

    double energy_density_3d(double lambda1, double lambda2, double lambda3) const {
        validate_positive(lambda1, "lambda1");
        validate_positive(lambda2, "lambda2");
        validate_positive(lambda3, "lambda3");

        double W = 0.0;
        for (const auto& term : terms_) {
            W += (term.mu / term.alpha) *
                 (std::pow(lambda1, term.alpha)
                + std::pow(lambda2, term.alpha)
                + std::pow(lambda3, term.alpha)
                - 3.0);
        }
        return W;
    }

    double dW_dlambda(double lambda) const {
        validate_positive(lambda, "lambda");

        double value = 0.0;
        for (const auto& term : terms_) {
            value += term.mu * std::pow(lambda, term.alpha - 1.0);
        }
        return value;
    }

    double uniaxial_nominal_stress_incompressible(double lambda) const {
        validate_positive(lambda, "lambda");

        const double lambda1 = lambda;
        const double lambda2 = 1.0 / std::sqrt(lambda);
        const double dW1 = dW_dlambda(lambda1);
        const double dW2 = dW_dlambda(lambda2);

        // For incompressibility, P_i = dW/dlambda_i - p/lambda_i.
        // The pressure p is chosen so the lateral nominal stress is zero.
        const double pressure = lambda2 * dW2;
        return dW1 - pressure / lambda1;
    }

private:
    std::vector<OgdenTerm> terms_;

    static void validate_positive(double value, const char* name) {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::runtime_error(std::string(name) + " must be positive and finite");
        }
    }
};

} // namespace ogden
