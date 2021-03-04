#include <vector>
#include <utility>
#include "json.hh"
#include "FMHardConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

FMHardConstraint::FMHardConstraint() {}

FMHardConstraint::FMHardConstraint(json &hc_json) {
    this->weight = hc_json["weight"];
    this->threshold = hc_json["threshold"];
    this->type = hc_json["type"];
    this->latest_cost = 0.0;
}

double FMHardConstraint::calculate_cost(std::vector<double> &struct_dose) {
    double struct_cost = 0.0;
    // LOWER LIMIT: No dose can be less than X Gy
    if (this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] < this->threshold) {
                double dose_diff = this->threshold - struct_dose[i];
                struct_cost += dose_diff * dose_diff;
            }
        }
    } else {  // UPPER LIMIT: No dose can be more than X Gy
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] > this->threshold) {
                double dose_diff = struct_dose[i] - this->threshold;
                struct_cost += dose_diff * dose_diff;
            }
        }
    }

    struct_cost = struct_cost * this->weight / struct_dose.size();
    this->latest_cost = struct_cost;
    return struct_cost;
}

double FMHardConstraint::calculate_gradient(std::vector<double> &struct_dose, std::vector<double> &cpt_dose) {
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] < this->threshold) {
                double thresholded_dose = this->threshold - struct_dose[i];
                gradient += -thresholded_dose * cpt_dose[i];
            }
        }
    } else {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] > this->threshold) {
                double thresholded_dose = struct_dose[i] - this->threshold;
                gradient += (thresholded_dose * cpt_dose[i]);
            }
        }
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double FMHardConstraint::calculate_hessian(std::vector<double> &struct_dose,
                                         std::vector<double> &cpt_one,
                                         std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] <= this->threshold) {
                hessian_contribution += cpt_one[i] * cpt_two[i];
            }
        }
    } else {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] >= this->threshold) {
                hessian_contribution += cpt_one[i] * cpt_two[i];
            }
        }
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}