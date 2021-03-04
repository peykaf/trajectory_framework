#include <vector>
#include <utility>
#include "json.hh"
#include "HardConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

HardConstraint::HardConstraint() {}

HardConstraint::HardConstraint(json &hc_json) {
    this->weight = hc_json["weight"];
    this->threshold = hc_json["threshold"];
    this->type = hc_json["type"];
    this->latest_cost = 0.0;
}

double HardConstraint::calculate_cost(std::vector<std::pair<int, double> > &struct_dose) {
    // Assumes dose array sorted in ascending dose values.
    double struct_cost = 0.0;

    // LOWER LIMIT: No dose can be less than X Gy
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double dose_diff = this->threshold - struct_dose[i].second;
            struct_cost += dose_diff * dose_diff;
        }
    } else {  // UPPER LIMIT: No dose can be more than X Gy
        for (int i = struct_dose.size() - 1; i >= 0 && struct_dose[i].second > this->threshold; i--) {
            double dose_diff = struct_dose[i].second - this->threshold;
            struct_cost += dose_diff * dose_diff;
        }
    }

    struct_cost = struct_cost * this->weight / struct_dose.size();
    this->latest_cost = struct_cost;
    return struct_cost;
}

double HardConstraint::calculate_unsorted_cost(std::vector<std::pair<int, double> > &struct_dose) {
    double struct_cost = 0.0;
    // LOWER LIMIT: No dose can be less than X Gy
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                double dose_diff = this->threshold - struct_dose[i].second;
                struct_cost += dose_diff * dose_diff;
            }
        }
    } else {  // UPPER LIMIT: No dose can be more than X Gy
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                double dose_diff = struct_dose[i].second - this->threshold;
                struct_cost += dose_diff * dose_diff;
            }
        }
    }

    struct_cost = struct_cost * this->weight / struct_dose.size();
    this->latest_cost = struct_cost;
    return struct_cost;
}

double HardConstraint::calculate_gradient(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_dose) {
    // Assumes dose array sorted in ascending dose values.
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            gradient += -thresholded_dose * cpt_dose[struct_dose[i].first];
        }
    } else {
        for (int i = struct_dose.size() - 1; i >=0 && struct_dose[i].second > this->threshold; i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            gradient += (thresholded_dose * cpt_dose[struct_dose[i].first]);
        }
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double HardConstraint::calculate_gradient(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_dose) {
    // Assumes dose array sorted in ascending dose values.
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            gradient += -thresholded_dose * cpt_dose.grid[struct_dose[i].first];
        }
    } else {
        for (int i = struct_dose.size() - 1; i >=0 && struct_dose[i].second > this->threshold; i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            gradient += (thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
        }
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double HardConstraint::calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_dose) {
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                double thresholded_dose = this->threshold - struct_dose[i].second;
                gradient += -thresholded_dose * cpt_dose[struct_dose[i].first];
            }
        }
    } else {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                double thresholded_dose = struct_dose[i].second - this->threshold;
                gradient += (thresholded_dose * cpt_dose[struct_dose[i].first]);
            }
        }
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double HardConstraint::calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_dose) {
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                double thresholded_dose = this->threshold - struct_dose[i].second;
                gradient += -thresholded_dose * cpt_dose.grid[struct_dose[i].first];
            }
        }
    } else {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                double thresholded_dose = struct_dose[i].second - this->threshold;
                gradient += (thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
            }
        }
    }

    return 2.0 * gradient * this->weight / struct_dose.size();
}

double HardConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second <= this->threshold; i++) {
            hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
        }
    } else {
        for (int i = struct_dose.size() - 1; i >= 0 && struct_dose[i].second >= this->threshold; i--) {
            hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
        }
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}

double HardConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_one, SparseDose &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second <= this->threshold; i++) {
            hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
        }
    } else {
        for (int i = struct_dose.size() - 1; i >= 0 && struct_dose[i].second >= this->threshold; i--) {
            hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
        }
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}


double HardConstraint::calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second <= this->threshold) {
                hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
            }
        }
    } else {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second >= this->threshold) {
                hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
            }
        }
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}

double HardConstraint::calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose, SparseDose &cpt_one, SparseDose &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second <= this->threshold) {
                hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
            }
        }
    } else {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second >= this->threshold) {
                hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
            }
        }
    }

    return 2.0 * hessian_contribution * this->weight / struct_dose.size();
}

void HardConstraint::calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                                       std::vector<double> &lag_mults) {
    // Want to maximize -gradient
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            lag_mults[struct_dose[i].first] += 2 * this->weight * thresholded_dose / struct_dose.size();
        }
    } else {
        for (int i = struct_dose.size() - 1; i >= 0 && struct_dose[i].second > this->threshold; i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            lag_mults[struct_dose[i].first] += -2 * this->weight * thresholded_dose / struct_dose.size();
        }
    }
}

void HardConstraint::calculate_unsorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                                       std::vector<double> &lag_mults) {
    if (this->type == LOWER_LIMIT) {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                double thresholded_dose = this->threshold - struct_dose[i].second;
                lag_mults[struct_dose[i].first] += 2 * this->weight * thresholded_dose / struct_dose.size();
            }
        }
    } else {
        for (int i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                double thresholded_dose = struct_dose[i].second - this->threshold;
                lag_mults[struct_dose[i].first] += -2 * this->weight * thresholded_dose / struct_dose.size();
            }
        }
    }
}

nlohmann::json HardConstraint::to_json() {
    nlohmann::json data;
    data["type"] = this->type;
    data["threshold"] = this->threshold;
    data["latest_cost"] = this->latest_cost;

    return data;
}