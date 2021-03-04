#include <vector>
#include <utility>
#include <math.h>
#include "json.hh"
#include "MeanConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

MeanConstraint::MeanConstraint() {}

MeanConstraint::MeanConstraint(json &mean_json) {
    this->weight = mean_json["weight"];
    this->threshold = mean_json["threshold"];
    this->type = mean_json["type"];
    this->percent_volume = mean_json["volume"];
    this->latest_cost = 0.0;
}


void MeanConstraint::set_volume(double pct_volume, int total_voxels) {
   if (pct_volume > 100.0) throw "Dose volume constraint must be <= 100%%";
   this->percent_volume = pct_volume;

   this->voxels_considered = (int) ((100.0 - percent_volume) / 100.0 * total_voxels);
   this->fraction = ((100.0 - percent_volume) / 100.0 * total_voxels);

   if (this->type == LOWER_LIMIT) {
      this->start = 0;
      this->end = this->voxels_considered;
   } else {
      this->start = (int) (total_voxels * (1.0 - this->percent_volume / 100.0));
      if (this->start > total_voxels) this->start = total_voxels;
      this->end = total_voxels;
      this->voxels_considered = this->end - this->start;
      this->fraction = this->voxels_considered;
   }
}

void MeanConstraint::set_volume(int total_voxels) {
   this->voxels_considered = (int) ((100.0 - this->percent_volume) / 100.0 * total_voxels);
   this->fraction = ((100.0 - this->percent_volume) / 100.0 * total_voxels);

   if (this->type == LOWER_LIMIT) {
      this->start = 0;
      this->end = this->voxels_considered;
   } else {
      this->start = (int) (total_voxels * (1.0 - this->percent_volume / 100.0));
      if (this->start > total_voxels) this->start = total_voxels;
      this->end = total_voxels;
      this->voxels_considered = this->end - this->start;
      this->fraction = this->voxels_considered;
   }
}

double MeanConstraint::calculate_mean_dose(
    std::vector<std::pair<int, double> > &struct_dose) {

    double mean_dose = 0.0;
    for (int i = this->start; i < this->end; i++) {
        mean_dose += struct_dose[i].second;
    }

    return mean_dose / this->fraction;
}

double MeanConstraint::calculate_cost(
    std::vector<std::pair<int, double> > &struct_dose) {
    double struct_cost = 0.0;

    double mean_dose = this->calculate_mean_dose(struct_dose);
    double difference = mean_dose - this->threshold;

    if(this->type == LOWER_LIMIT) {
        if (mean_dose < this->threshold) {
            struct_cost += difference * difference;
        }
    }
    else {
        if (mean_dose > this->threshold) {
            struct_cost += difference * difference;
        }
    }

    struct_cost = struct_cost * this->weight;
    this->latest_cost = struct_cost;
    return struct_cost;
}

double MeanConstraint::sum_struct_aperture_dose(
    std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_dose) {

    double sum_dose = 0.0;
    for (int i = this->start; i < this->end; i++) {
        sum_dose += cpt_dose.grid[struct_dose[i].first];
    }

    return sum_dose / this->fraction;
}

double MeanConstraint::sum_struct_aperture_dose(
    std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_dose) {

    double sum_dose = 0.0;
    for (int i = this->start; i < this->end; i++) {
        sum_dose += cpt_dose[struct_dose[i].first];
    }

    return sum_dose / this->fraction;
}

double MeanConstraint::calculate_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_dose) {

    double mean_dose = this->calculate_mean_dose(struct_dose);
    double difference = mean_dose - this->threshold;

    double gradient = 0.0;
    if(this->type == LOWER_LIMIT) {
        if (mean_dose < this->threshold) {
            double ap_dose = this->sum_struct_aperture_dose(struct_dose, cpt_dose);
            gradient += 2 * difference * ap_dose;
        }
    }
    else {
        if (mean_dose > this->threshold) {
            double ap_dose = this->sum_struct_aperture_dose(struct_dose, cpt_dose);
            gradient += 2 * difference * ap_dose;
        }
    }

    return gradient * this->weight;
}

double MeanConstraint::calculate_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_dose) {

    double mean_dose = this->calculate_mean_dose(struct_dose);
    double difference = mean_dose - this->threshold;

    double gradient = 0.0;
    if(this->type == LOWER_LIMIT) {
        if (mean_dose < this->threshold) {
            double ap_dose = this->sum_struct_aperture_dose(struct_dose, cpt_dose);
            gradient += 2 * difference * ap_dose;
        }
    }
    else {
        if (mean_dose > this->threshold) {
            double ap_dose = this->sum_struct_aperture_dose(struct_dose, cpt_dose);
            gradient += 2 * difference * ap_dose;
        }
    }

    return gradient * this->weight;
}



double MeanConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_one,
    SparseDose &cpt_two) {
    double hessian_contribution = 0.0;

    double mean_dose = this->calculate_mean_dose(struct_dose);

    // LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if (this->type == LOWER_LIMIT) {
        if (mean_dose <= this->threshold) {
            hessian_contribution += 2 * this->sum_struct_aperture_dose(struct_dose, cpt_one)
                                      * this->sum_struct_aperture_dose(struct_dose, cpt_two);
        }
    } else {
        if (mean_dose >= this->threshold) {
            hessian_contribution += 2 * this->sum_struct_aperture_dose(struct_dose, cpt_one)
                                      * this->sum_struct_aperture_dose(struct_dose, cpt_two);
        }
    }

    return hessian_contribution * this->weight;
}

void MeanConstraint::calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &price_map) {
    double mean_dose = this->calculate_mean_dose(struct_dose);
    double difference = mean_dose - this->threshold;
    double constant_lagmult = -2 * difference * this->weight / this->fraction;

    if(this->type == LOWER_LIMIT) {
        if (mean_dose < this->threshold) {
            for (int i = this->start; i < this->end; i++) {
                price_map[struct_dose[i].first] += constant_lagmult;
            }
        }
    } else {
        if (mean_dose > this->threshold) {
            for (int i = this->start; i < this->end; i++) {
                price_map[struct_dose[i].first] += constant_lagmult;
            }
        }
    }
}

double MeanConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_one,
    std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    double mean_dose = this->calculate_mean_dose(struct_dose);

    // LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if (this->type == LOWER_LIMIT) {
        if (mean_dose <= this->threshold) {
            hessian_contribution += 2 * this->sum_struct_aperture_dose(struct_dose, cpt_one)
                                      * this->sum_struct_aperture_dose(struct_dose, cpt_two);
        }
    } else {
        if (mean_dose >= this->threshold) {
            hessian_contribution += 2 * this->sum_struct_aperture_dose(struct_dose, cpt_one)
                                      * this->sum_struct_aperture_dose(struct_dose, cpt_two);
        }
    }

    return hessian_contribution * this->weight;
}

nlohmann::json MeanConstraint::to_json() {
    nlohmann::json data;
    data["type"] = this->type;
    data["threshold"] = this->threshold;
    data["latest_cost"] = this->latest_cost;
    data["percent_volume"] = this->percent_volume;

    return data;
}