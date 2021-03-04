#include <vector>
#include <utility>
#include "json.hh"
#include "FMDVConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

FMDVConstraint::FMDVConstraint() {}

FMDVConstraint::FMDVConstraint(json &dv_json) {
    this->weight = dv_json["weight"];
    this->threshold = dv_json["threshold"];
    this->type = dv_json["type"];
    this->percent_volume = dv_json["volume"];
    this->latest_cost = 0.0;
}

void FMDVConstraint::set_volume(int total_voxels) {
    this->start = (int) ((100.0 - this->percent_volume) / 100.0 * total_voxels);
    if (this->start > total_voxels) this->start = total_voxels;
    //this->voxels_considered = total_voxels - this->start;
    this->voxels_considered = total_voxels;
    if (this->type == UPPER_LIMIT) {
        this->allowed_violators = (int) (this->percent_volume / 100.0) * total_voxels;
    } else {
        this->allowed_violators = this->start;
    }
}

void FMDVConstraint::set_volume(double pct_volume, int total_voxels) {
    if (pct_volume > 100.0) throw "Dose volume constraint must be <= 100%%";
    this->percent_volume = pct_volume;

    this->start = (int) ((100.0 - percent_volume) / 100.0 * total_voxels);
    if (this->start > total_voxels) this->start = total_voxels;
    //this->voxels_considered = total_voxels - this->start;
    this->voxels_considered = total_voxels;

    if (this->type == UPPER_LIMIT) {
        this->allowed_violators = (int) (this->percent_volume / 100.0) * total_voxels;
    } else {
        this->allowed_violators = this->start;
    }
}

void FMDVConstraint::find_violator_doses(std::vector<double> &struct_dose) {
    std::vector<std::pair<int, double> > dose_values;

    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] < this->threshold) {
                dose_values.push_back(std::make_pair(i, struct_dose[i]));
            }
        }
    } else {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i] > this->threshold) {
                dose_values.push_back(std::make_pair(i, struct_dose[i]));
            }
        }
    }

    std::sort(dose_values.begin(),
              dose_values.end(),
              [](const std::pair<int, double>& lhs,const std::pair<int, double>& rhs){ return lhs.second < rhs.second; });

    this->violator_doses = dose_values;
}

double FMDVConstraint::calculate_cost(
    std::vector<double> &struct_dose) {
    double struct_cost = 0.0;

    int counter = 0;

    find_violator_doses(struct_dose);
    if (this->type == LOWER_LIMIT) {
        if (this->violator_doses.size() > this->allowed_violators) {
            for (int i = this->allowed_violators; i < this->violator_doses.size(); i++) {
                double dose_diff = this->threshold - this->violator_doses[i].second;
                struct_cost += dose_diff * dose_diff;
            }
        }
    } else {
        if (this->violator_doses.size() > this->allowed_violators) {
            size_t start = this->violator_doses.size() - this->allowed_violators - 1;
            for (int i = start; i >= 0; i--) {
                double dose_diff = this->violator_doses[i].second - this->threshold;
                struct_cost += dose_diff * dose_diff;
            }
        }
    }

    struct_cost = struct_cost * this->weight / this->voxels_considered;
    this->latest_cost = struct_cost;
    return struct_cost;

}

double FMDVConstraint::calculate_gradient(
    std::vector<double> &struct_dose,
    std::vector<double> &cpt_dose) {
    double gradient = 0.0;

    int counter = 0;

    if (this->type == LOWER_LIMIT) {
        if (this->violator_doses.size() > this->allowed_violators) {
            for (int i = this->allowed_violators; i < this->violator_doses.size(); i++) {
                double dose_diff = this->threshold - this->violator_doses[i].second;
                gradient += -(dose_diff * cpt_dose[this->violator_doses[i].first]);
            }
        }
    } else {
        if (this->violator_doses.size() > this->allowed_violators) {
            int start = this->violator_doses.size() - this->allowed_violators - 1;
            for (int i = start; i >= 0; i--) {
                double dose_diff = this->violator_doses[i].second - this->threshold;
                gradient += (dose_diff * cpt_dose[this->violator_doses[i].first]);
            }
        }
    }

    return 2.0 * gradient * this->weight / this->voxels_considered;
}

double FMDVConstraint::calculate_hessian(std::vector<double> &struct_dose,
    std::vector<double> &cpt_one,
    std::vector<double> &cpt_two) {

    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        if (this->violator_doses.size() > this->allowed_violators) {
            for (int i = this->allowed_violators; i < this->violator_doses.size(); i++) {
                hessian_contribution += cpt_one[this->violator_doses[i].first] * cpt_two[this->violator_doses[i].first];
            }
        }
    } else {
        if (this->violator_doses.size() > this->allowed_violators) {
            int start = this->violator_doses.size() - this->allowed_violators - 1;
            for (int i = start; i >= 0; i--) {
                hessian_contribution += cpt_one[this->violator_doses[i].first] * cpt_two[this->violator_doses[i].first];
            }
        }
    }

    return 2.0 * hessian_contribution * this->weight / this->voxels_considered;
}
