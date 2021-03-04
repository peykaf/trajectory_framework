#include <vector>
#include <utility>
#include "json.hh"
#include "DVConstraint.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

DVConstraint::DVConstraint() {}

DVConstraint::DVConstraint(json &dv_json) {
    this->weight = dv_json["weight"];
    this->threshold = dv_json["threshold"];
    this->type = dv_json["type"];
    this->percent_volume = dv_json["volume"];
    this->latest_cost = 0.0;
}

void DVConstraint::set_volume(int total_voxels) {
   this->start = (int) ((100.0 - this->percent_volume) / 100.0 * total_voxels);
   if (this->start > total_voxels) this->start = total_voxels;
   //this->voxels_considered = total_voxels - this->start;
   this->voxels_considered = total_voxels;
}

void DVConstraint::set_volume(double pct_volume, int total_voxels) {
   if (pct_volume > 100.0) throw "Dose volume constraint must be <= 100%%";
   this->percent_volume = pct_volume;

   this->start = (int) ((100.0 - percent_volume) / 100.0 * total_voxels);
   if (this->start > total_voxels) this->start = total_voxels;
   //this->voxels_considered = total_voxels - this->start;
   this->voxels_considered = total_voxels;
}

double DVConstraint::calculate_unsorted_cost(
    std::vector<std::pair<int, double> > &struct_dose){
    double struct_cost = 0.0;

    int counter = 0;
    int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double dose_diff = this->threshold - struct_dose[i].second;
                    struct_cost += dose_diff * dose_diff;
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double dose_diff = struct_dose[i].second - this->threshold;
                    struct_cost += dose_diff * dose_diff;
                }
            }
        }
    }
    struct_cost = struct_cost * this->weight / this->voxels_considered;
    this->latest_cost = struct_cost;
    return struct_cost;

}

double DVConstraint::calculate_sorted_cost(
    std::vector<std::pair<int, double> > &struct_dose) {
    double struct_cost = 0.0;
    // LOWER LIMIT (no more than Y% of volume can get LESS than X Gy) (reworded: we want Y% of the volume to get at least X Gy)
    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() &&
            struct_dose[i].second < this->threshold;
            i++) {
            double dose_diff = this->threshold - struct_dose[i].second;
            struct_cost += dose_diff * dose_diff;
        }
    } else {  // UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy) (reworded: we want to limit the number of voxels getting X Gy to Y%)
        for (int i = this->start - 1; i >= 0 &&
             struct_dose[i].second > this->threshold;
             i--) {
            double dose_diff = struct_dose[i].second - this->threshold;
            struct_cost += dose_diff * dose_diff;
        }
    }

    struct_cost = struct_cost * this->weight / this->voxels_considered;
    this->latest_cost = struct_cost;
    return struct_cost;
}

double DVConstraint::calculate_cost(
    std::vector<std::pair<int, double> > &struct_dose, bool sorted) {
    if (sorted) {
        return calculate_sorted_cost(struct_dose);
    } else {
        return calculate_unsorted_cost(struct_dose);
    }
}


double DVConstraint::calculate_unsorted_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_dose) {
    double gradient = 0.0;

    int counter = 0; int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = this->threshold - struct_dose[i].second;
                    gradient += (thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = struct_dose[i].second - this->threshold;
                    gradient += (thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
                }
            }
        }
    }

    return 2.0 * gradient * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_sorted_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_dose) {
    // Assumes dose array sorted in ascending dose values, and mask array sorted in the same fashion.
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            gradient += -(thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
        }
    } else {
        for (int i = this->start - 1; i >= 0 &&
             struct_dose[i].second > this->threshold;
             i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            gradient += (thresholded_dose * cpt_dose.grid[struct_dose[i].first]);
        }
    }

    return 2.0 * gradient * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_dose,
    bool sorted) {

    if (sorted) {
        return calculate_sorted_gradient(struct_dose, cpt_dose);
    } else {
        return calculate_unsorted_gradient(struct_dose, cpt_dose);
    }
}

double DVConstraint::calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_one,
    SparseDose &cpt_two) {

    double hessian_contribution = 0.0;

    int counter = 0; int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second <= this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second >= this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
                }
            }
        }
    }

    return 2.0 * hessian_contribution * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_sorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_one,
    SparseDose &cpt_two) {
    double hessian_contribution = 0.0;

    // LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() && struct_dose[i].second <= this->threshold; i++) {
            hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
        }
    } else {
        for (int i = this->start - 1; i >= 0 &&
             struct_dose[i].second >= this->threshold;
             i--) {
            hessian_contribution += cpt_one.grid[struct_dose[i].first] * cpt_two.grid[struct_dose[i].first];
        }
    }

    return 2.0 * hessian_contribution * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
    SparseDose &cpt_one,
    SparseDose &cpt_two,
    bool sorted) {
    if (sorted) {
        return calculate_sorted_hessian(struct_dose, cpt_one, cpt_two);
    } else {
        return calculate_unsorted_hessian(struct_dose, cpt_one, cpt_two);
    }
}

void DVConstraint::calculate_unsorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &price_map) {
    int counter = 0; int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = this->threshold - struct_dose[i].second;
                    price_map[struct_dose[i].first] += 2 * this->weight * thresholded_dose / this->voxels_considered;
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = struct_dose[i].second - this->threshold;
                    price_map[struct_dose[i].first] += -2 * this->weight * thresholded_dose / this->voxels_considered;
                }
            }
        }
    }
}

void DVConstraint::calculate_sorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &price_map) {
    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            price_map[struct_dose[i].first] += 2 * this->weight * thresholded_dose / this->voxels_considered;
        }
    } else {
        for (int i = this->start - 1; i >= 0 &&
             struct_dose[i].second > this->threshold;
             i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            price_map[struct_dose[i].first] += -2 * this->weight * thresholded_dose / this->voxels_considered;
        }
    }
}

void DVConstraint::calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &price_map, bool sorted) {

    if (sorted) {
        calculate_sorted_lag_mults(struct_dose, price_map);
    } else {
        calculate_unsorted_lag_mults(struct_dose, price_map);
    }
}

// vector<double> dose array functions below here

double DVConstraint::calculate_unsorted_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_dose) {

    double gradient = 0.0;

    int counter = 0; int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = this->threshold - struct_dose[i].second;
                    gradient += (thresholded_dose * cpt_dose[struct_dose[i].first]);
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    double thresholded_dose = struct_dose[i].second - this->threshold;
                    gradient += (thresholded_dose * cpt_dose[struct_dose[i].first]);
                }
            }
        }
    }

    return 2.0 * gradient * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_sorted_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_dose) {
    // Assumes dose array sorted in ascending dose values, and mask array sorted in the same fashion.
    double gradient = 0.0;
    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            double thresholded_dose = this->threshold - struct_dose[i].second;
            gradient += -(thresholded_dose * cpt_dose[struct_dose[i].first]);
        }
    } else {
        for (int i = this->start - 1; i >= 0 && struct_dose[i].second > this->threshold; i--) {
            double thresholded_dose = struct_dose[i].second - this->threshold;
            gradient += (thresholded_dose * cpt_dose[struct_dose[i].first]);
        }
    }

    return 2.0 * gradient * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_gradient(
    std::vector<std::pair<int, double> > &struct_dose,
    std::vector<double> &cpt_dose,
    bool sorted) {

    if (sorted) {
        return calculate_sorted_gradient(struct_dose, cpt_dose);
    } else {
        return calculate_unsorted_gradient(struct_dose, cpt_dose);
    }
}

double DVConstraint::calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    int counter = 0;
    int min_voxels = this->start;
    //LOWER LIMIT (no more than Y% of volume can get LESS than X Gy)
    if(this->type == LOWER_LIMIT) {
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second < this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
                }
            }
        }
    } else { //UPPER LIMIT: (no more than Y% of volume can get MORE than X Gy)
        for (size_t i = 0; i < struct_dose.size(); i++) {
            if (struct_dose[i].second > this->threshold) {
                if (counter < min_voxels) {
                    counter++;
                } else {
                    hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
                }
            }
        }
    }
    return 2.0 * hessian_contribution * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_sorted_hessian(std::vector<std::pair<int, double> > &struct_dose, std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian_contribution = 0.0;

    if (this->type == LOWER_LIMIT) {
        for (int i = this->start + 1; i < struct_dose.size() && struct_dose[i].second < this->threshold; i++) {
            hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
        }
    } else {
        for (int i = this->start - 1; i >= 0 &&
             struct_dose[i].second > this->threshold;
             i--) {
            hessian_contribution += cpt_one[struct_dose[i].first] * cpt_two[struct_dose[i].first];
        }
    }
    return 2.0 * hessian_contribution * this->weight / this->voxels_considered;
}

double DVConstraint::calculate_hessian(std::vector<std::pair<int, double> > &struct_dose, std:: vector<double> &cpt_one, std::vector<double> &cpt_two, bool sorted) {
    if (sorted) {
        return calculate_sorted_hessian(struct_dose, cpt_one, cpt_two);
    } else {
        return calculate_unsorted_hessian(struct_dose, cpt_one, cpt_two);
    }
}

nlohmann::json DVConstraint::to_json() {
    nlohmann::json data;
    data["type"] = this->type;
    data["threshold"] = this->threshold;
    data["latest_cost"] = this->latest_cost;
    data["percent_volume"] = this->percent_volume;

    return data;
}