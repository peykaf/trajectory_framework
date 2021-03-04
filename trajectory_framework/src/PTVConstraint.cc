#include <vector>
#include <utility>
#include "json.hh"
#include "PTVConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

PTVConstraint::PTVConstraint() {}

PTVConstraint::PTVConstraint(json &hc_json) {
    this->weight = hc_json["weight"];
    this->threshold = hc_json["threshold"];
    this->latest_cost = 0.0;
}

double PTVConstraint::calculate_cost(std::vector<std::pair<int, double> > &struct_dose) {
    // Assumes dose array sorted in ascending dose values.
    double struct_cost = 0.0;

    for (int i = 0; i < struct_dose.size(); i++) {
        double dose_diff = struct_dose[i].second - this->threshold;
        struct_cost += dose_diff * dose_diff;
    }

    struct_cost = struct_cost * this->weight / struct_dose.size();
    this->latest_cost = struct_cost;
    return struct_cost;
}

double PTVConstraint::calculate_gradient(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_dose) {
    // Assumes dose array sorted in ascending dose values.
    double gradient = 0.0;
    for (int i = 0; i < struct_dose.size(); i++) {
        double thresholded_dose = struct_dose[i].second - this->threshold;
        gradient += thresholded_dose * cpt_dose[struct_dose[i].first];
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double PTVConstraint::calculate_gradient(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_dose) {
    // Assumes dose array sorted in ascending dose values.
    double gradient = 0.0;
    for (int i = 0; i < struct_dose.size(); i++) {
        double thresholded_dose = struct_dose[i].second - this->threshold;
        gradient += thresholded_dose * cpt_dose.grid[struct_dose[i].first];
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double PTVConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    for (int i = 0; i < struct_dose.size(); i++) {
        hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}

double PTVConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_one, SparseDose &cpt_two) {
    double hessian_contribution = 0.0;

    for (int i = 0; i < struct_dose.size(); i++) {
        hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}


void PTVConstraint::calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                                       std::vector<double> &lag_mults) {
    // Want to maximize -gradient
    for (int i = 0; i < struct_dose.size(); i++) {
        double thresholded_dose = struct_dose[i].second - this->threshold;
        lag_mults[struct_dose[i].first] += -2 * this->weight * thresholded_dose / struct_dose.size();
    }
}

nlohmann::json PTVConstraint::to_json() {
    nlohmann::json data;
    data["threshold"] = this->threshold;
    data["latest_cost"] = this->latest_cost;

    return data;
}