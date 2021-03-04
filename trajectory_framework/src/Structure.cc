#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include "json.hh"

#include "string_utils.hh"
#include "Structure.hh"
#include "PTVConstraint.hh"
#include "HardConstraint.hh"
#include "DVConstraint.hh"

using json = nlohmann::json;

Structure::Structure() {
    this->sorted = true;
    this->voxel_volume = 1.0;
    this->min_resolution = 2000;
    this->target = false;
}

Structure::Structure(json &json_struct) : Structure() {
    this->name = json_struct["name"];
    this->target = json_struct.value("target", false);

    std::vector<int> mask = json_struct["mask"].get<std::vector<int> >();

    // Allow sparse sampling of structure voxels to make optimisation go faster.
    bool downsample = false;
    double frac_volume = 1.0;
    if (json_struct.find("downsample") != json_struct.end()) {
        downsample = json_struct["downsample"];
        frac_volume = json_struct["pct_volume"].get<double>() / 100.0;
    }

    if (downsample) {
        int min_voxels = std::min(int(mask.size()), this->min_resolution);
        int final_num_voxels = int(frac_volume * mask.size());
        final_num_voxels = std::max(final_num_voxels, min_voxels);

        // Shuffle the elements of the mask vector to randomly sample a subset.
        // Always use the same (arbitrary) seed for reproducibility.
        std::mt19937 g(100);

        std::shuffle(mask.begin(), mask.end(), g);

        for (size_t i = 0; i < final_num_voxels; i ++) {
            this->masked_dose.push_back(std::make_pair(mask[i], 0.0));
        }

        std::cout << "Downsampling " << this->name << std::endl;
        std::cout << "Using " << this->masked_dose.size() << " voxels instead of the full " << mask.size() << std::endl;
    } else {
        for (auto voxel_num: mask) {
            this->masked_dose.push_back(std::make_pair(voxel_num, 0.0));
        }
        std::cout << "Using " << this->masked_dose.size() << " voxels" << std::endl;
    }

    this->total_voxels = this->masked_dose.size();

    if (json_struct.find("ptv_constraints") != json_struct.end()) {
        for (auto& ptvc_json : json_struct["ptv_constraints"]) {
            PTVConstraint ptv_constraint(ptvc_json);
            this->ptv_constraints.push_back(ptv_constraint);
        }
    }

    if (json_struct.find("hard_constraints") != json_struct.end()) {
        for (auto& hc_json : json_struct["hard_constraints"]) {
            HardConstraint hard_constraint(hc_json);
            this->hard_constraints.push_back(hard_constraint);
        }
    }

    if (json_struct.find("dv_constraints") != json_struct.end()) {
        for (auto& dv_json : json_struct["dv_constraints"]) {
            DVConstraint dv_constraint(dv_json);
            dv_constraint.set_volume(this->total_voxels);
            this->dv_constraints.push_back(dv_constraint);
        }
    }

    if (json_struct.find("mean_constraints") != json_struct.end()) {
        for (auto& mean_json : json_struct["mean_constraints"]) {
            MeanConstraint mean_constraint(mean_json);
            mean_constraint.set_volume(this->total_voxels);
            this->mean_constraints.push_back(mean_constraint);
        }
    }

    if (json_struct.find("voxel_volume") != json_struct.end()) {
        this->voxel_volume = json_struct["voxel_volume"];
    }

}

void Structure::scale_weights(double scale) {
    for (auto &pc : this->ptv_constraints) {
        pc.weight *= scale;
    }

    for (auto &hc : this->hard_constraints) {
        hc.weight *= scale;
    }

    for (auto &dc : this->dv_constraints) {
        dc.weight *= scale;
    }

    for (auto &mc : this->mean_constraints) {
        mc.weight *= scale;
    }
}

double Structure::get_min_dose() {
    double min_dose = 0.0;
    for (auto &hc : this->hard_constraints) {
        if (hc.type == LOWER_LIMIT) min_dose = hc.threshold;
    }

    return min_dose;
}

void Structure::output_constraints() {
    std::cout << this->name << std::endl;
    if (this->ptv_constraints.size() > 0) {
        std::cout << "PTV constraints" << std::endl;
        for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
            std::cout << this->ptv_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Weight: " << this->ptv_constraints[i].weight << std::endl;
        }
    }

    if (this->hard_constraints.size() > 0) {
        std::cout << "Hard constraints" << std::endl;
        for (size_t i = 0; i < this->hard_constraints.size(); i++) {
            if (this->hard_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->hard_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Weight: " << this->hard_constraints[i].weight << std::endl;
        }
    }

    if (this->dv_constraints.size() > 0) {
        std::cout << "DV constraints" << std::endl;
        for (size_t i = 0; i < this->dv_constraints.size(); i++) {
            if (this->dv_constraints[i].type == UPPER_LIMIT) {
                std:: cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->dv_constraints[i].threshold << "Gy @ "
                << this->dv_constraints[i].percent_volume << "%%" << std::endl;
            std::cout << "Weight: " << this->dv_constraints[i].weight << std::endl;
        }
    }

    if (this->mean_constraints.size() > 0) {
        std::cout << "Mean constraints" << std::endl;
        for (size_t i = 0; i < this->mean_constraints.size(); i++) {
            if (this->mean_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->mean_constraints[i].threshold << "Gy" << std::endl;
            std::cout << this->mean_constraints[i].threshold << "Gy @ "
                << this->mean_constraints[i].percent_volume << "%%" << std::endl;

            std::cout << "Weight: " << this->mean_constraints[i].weight << std::endl;
        }
    }
}

json Structure::to_json() {
    json struct_json;

    struct_json["name"] = this->name;

    if (this->ptv_constraints.size() > 0) {
        json ptvc_constraints;
        for (auto &hc : this->ptv_constraints) {
            json ptv_json;
            ptv_json["threshold"] = hc.threshold;
            ptv_json["weight"] = hc.weight;
            ptvc_constraints.push_back(ptv_json);
        }
        struct_json["ptv_constraints"] = ptvc_constraints;
    }

    if (this->hard_constraints.size() > 0) {
        json hc_constraints;
        for (auto &hc : this->hard_constraints) {
            json hc_json;
            hc_json["threshold"] = hc.threshold;
            hc_json["weight"] = hc.weight;
            hc_json["constraint_type"] = hc.type;
            hc_constraints.push_back(hc_json);
        }
        struct_json["hard_constraints"] = hc_constraints;
    }

    if (this->dv_constraints.size() > 0) {
        json dv_json_constraints;
        for (auto &dv : this->dv_constraints) {
            json dv_json;
            dv_json["threshold"] = dv.threshold;
            dv_json["weight"] = dv.weight;
            dv_json["constraint_type"] = dv.type;
            dv_json["volume"] = dv.percent_volume;

            dv_json_constraints.push_back(dv_json);
        }
        struct_json["dv_constraints"] = dv_json_constraints;
    }

    if (this->mean_constraints.size() > 0) {
        json mean_json_constraints;
        for (auto &mean : this->mean_constraints) {
            json mean_json;
            mean_json["threshold"] = mean.threshold;
            mean_json["weight"] = mean.weight;
            mean_json["constraint_type"] = mean.type;
            mean_json["volume"] = mean.percent_volume;
            mean_json_constraints.push_back(mean_json);
        }
        struct_json["mean_constraints"] = mean_json_constraints;
    }

    return struct_json;

}

void Structure::write_constraints(std::ofstream &outfile) {
    outfile << "Structure" << std::endl;
    outfile << "Name = " << this->name << std::endl;
    outfile << "Num ptv constraints = " << this->ptv_constraints.size() << std::endl;
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        outfile << "Threshold = " << this->ptv_constraints[i].threshold << std::endl;
        outfile << "Weight = " << this->ptv_constraints[i].weight << std::endl;
    }

    outfile << "Num hard constraints = " << this->hard_constraints.size() << std::endl;
    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        outfile << "constraint type = ";
        if (this->hard_constraints[i].type == UPPER_LIMIT) {
            outfile << "Upper limit";
        } else {
            outfile << "Lower limit";
        }
        outfile << std::endl;
        outfile << "Threshold = " << this->hard_constraints[i].threshold << std::endl;
        outfile << "Weight = " << this->hard_constraints[i].weight << std::endl;
    }

    outfile << "Num dv constraints = " << this->dv_constraints.size() << std::endl;
    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        outfile << "constraint type = ";
        if (this->dv_constraints[i].type == UPPER_LIMIT) {
            outfile << "Upper limit";
        } else {
            outfile << "Lower limit";
        }
        outfile << std::endl;
        outfile << "Threshold = " << this->dv_constraints[i].threshold << std::endl;
        outfile << "Percent volume = " << this->dv_constraints[i].percent_volume << std::endl;
        outfile << "Weight = " << this->dv_constraints[i].weight << std::endl;
    }

    outfile << "Num mean constraints = " << this->mean_constraints.size() << std::endl;
    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        outfile << "constraint type = ";
        if (this->mean_constraints[i].type == UPPER_LIMIT) {
            outfile << "Upper limit";
        } else {
            outfile << "Lower limit";
        }
        outfile << std::endl;
        outfile << "Threshold = " << this->mean_constraints[i].threshold << std::endl;
        outfile << "Percent volume = " << this->mean_constraints[i].percent_volume << std::endl;
        outfile << "Weight = " << this->mean_constraints[i].weight << std::endl;
    }
}

void Structure::output_cost() {
    std::cout << this->name << std::endl;
    // Dose inside structures is sorted.
    std::cout << "Min dose: " << this->masked_dose[0].second << std::endl;
    std::cout << "Max dose: " << this->masked_dose[this->masked_dose.size() - 1].second << std::endl;
    if (this->ptv_constraints.size() > 0) {
        std::cout << "PTV constraints" << std::endl;
        for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
            std::cout << this->ptv_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Cost: " << this->ptv_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }

    if (this->hard_constraints.size() > 0) {
        std::cout << "Hard constraints" << std::endl;
        for (size_t i = 0; i < this->hard_constraints.size(); i++) {
            if (this->hard_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->hard_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Cost: " << this->hard_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }

    if (this->dv_constraints.size() > 0) {
        std::cout << "DV constraints" << std::endl;
        for (size_t i = 0; i < this->dv_constraints.size(); i++) {
            if (this->dv_constraints[i].type == UPPER_LIMIT) {
                std:: cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->dv_constraints[i].threshold << "Gy @ "
                << this->dv_constraints[i].percent_volume << "%" << std::endl;
            std::cout << "Cost: " << this->dv_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }

    if (this->mean_constraints.size() > 0) {
        std::cout << "Mean constraints" << std::endl;
        for (size_t i = 0; i < this->mean_constraints.size(); i++) {
            if (this->mean_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->mean_constraints[i].threshold << "Gy" << std::endl;
            std::cout << this->mean_constraints[i].threshold << "Gy @ "
                << this->mean_constraints[i].percent_volume << "%" << std::endl;

            std::cout << "Cost: " << this->mean_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }
}

nlohmann::json Structure::update_data() {
    nlohmann::json iter_data;
    iter_data["ptv_constraints"] = std::vector<json>();
    iter_data["hard_constraints"] = std::vector<json>();
    iter_data["dv_constraints"] = std::vector<json>();
    iter_data["mean_constraints"] = std::vector<json>();
    iter_data["name"] = this->name;

    for (auto &constraint : this->ptv_constraints) {
        iter_data["ptv_constraints"].push_back(constraint.to_json());
    }

    for (auto &constraint : this->hard_constraints) {
        iter_data["hard_constraints"].push_back(constraint.to_json());
    }

    for (auto &constraint : this->dv_constraints) {
        iter_data["dv_constraints"].push_back(constraint.to_json());
    }

    for (auto &constraint : this->mean_constraints) {
        iter_data["mean_constraints"].push_back(constraint.to_json());
    }

    return iter_data;
}

double Structure::calculate_cost() {
    double struct_cost = 0.0;
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        struct_cost += this->ptv_constraints[i].calculate_cost(this->masked_dose);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        struct_cost += this->hard_constraints[i].calculate_cost(this->masked_dose);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        struct_cost += this->dv_constraints[i].calculate_cost(this->masked_dose, this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        struct_cost += this->mean_constraints[i].calculate_cost(this->masked_dose);
    }

    return struct_cost;
}


double Structure::calculate_gradient(std::vector<double> &cpt_dose) {
    double gradient = 0.0;

    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        gradient += this->ptv_constraints[i].calculate_gradient(this->masked_dose, cpt_dose);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        gradient += this->hard_constraints[i].calculate_gradient(this->masked_dose, cpt_dose);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        gradient += this->dv_constraints[i].calculate_gradient(this->masked_dose, cpt_dose, this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        gradient += this->mean_constraints[i].calculate_gradient(this->masked_dose, cpt_dose);
    }

    return gradient;
}


double Structure::calculate_gradient(SparseDose &dose) {
    double gradient = 0.0;
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        gradient += this->ptv_constraints[i].calculate_gradient(this->masked_dose, dose);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        gradient += this->hard_constraints[i].calculate_gradient(this->masked_dose, dose);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        gradient += this->dv_constraints[i].calculate_gradient(this->masked_dose, dose, this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        gradient += this->mean_constraints[i].calculate_gradient(this->masked_dose, dose);
    }


    return gradient;
}


double Structure::calculate_hessian(std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian = 0.0;
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        hessian += this->ptv_constraints[i].calculate_hessian(this->masked_dose,
                                                              cpt_one,
                                                              cpt_two);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        hessian += this->hard_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        hessian += this->dv_constraints[i].calculate_hessian(this->masked_dose,
                                                             cpt_one,
                                                             cpt_two,
                                                             this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        hessian += this->mean_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    return hessian;
}

double Structure::calculate_hessian(SparseDose &cpt_one, SparseDose &cpt_two) {
    double hessian = 0.0;
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
       hessian += this->ptv_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
       hessian += this->hard_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        hessian += this->dv_constraints[i].calculate_hessian(this->masked_dose,
                                                             cpt_one,
                                                             cpt_two,
                                                             this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        hessian += this->mean_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    return hessian;
}


void Structure::calculate_lag_mults(std::vector<double> &lag_mults) {
    for (size_t i = 0; i < this->ptv_constraints.size(); i++) {
        this->ptv_constraints[i].calculate_lag_mults(this->masked_dose, lag_mults);
    }

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        this->hard_constraints[i].calculate_lag_mults(this->masked_dose, lag_mults);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        this->dv_constraints[i].calculate_lag_mults(this->masked_dose, lag_mults, this->sorted);
    }

    for (size_t i = 0; i < this->mean_constraints.size(); i++) {
        this->mean_constraints[i].calculate_lag_mults(this->masked_dose, lag_mults);
    }
}

double Structure::dose_for_volume(double volume, bool percent) {
    // D90%, D1cc, etc
    // Volume in cc or percent depending on boolean "percent" variable

    if (this->cdvh.size() == 0) this->calculate_dvh();

    double dose = 0.0;
    for (int i = 0; i < this->cdvh.size(); i++) {
        double bin_volume = this->cdvh[i];
        if (percent) bin_volume = bin_volume / this->cdvh[0] * 100.0;
        if (bin_volume < volume) {
            if (i == 0) return 0.0;
            double vol_difference, full_difference;
            if (percent) {
                vol_difference = volume - bin_volume;
                full_difference = this->cdvh[i - 1] / this->cdvh[0] * 100.0 - bin_volume;
            } else {
                vol_difference = volume - bin_volume;
                full_difference = this->cdvh[i - 1] - bin_volume;
            }

            double fractional_difference = vol_difference / full_difference;
            double bin_size = 1.0;
            dose = (i + fractional_difference * bin_size) / 100.0;

            break;

        }
    }

    if (dose == 0.0) dose = this->cdvh.size();

    return dose;
}

double Structure::volume_for_dose(double dose, bool percent) {
    // Dose is input in Gy.
    dose = dose * 100.0;
    int bin = static_cast<int>(dose);
    double difference = dose - bin;
    double volume = this->cdvh[bin] + difference * (this->cdvh[bin + 1] - this->cdvh[bin]);

    if (percent) volume = volume / this->cdvh[0] * 100;

    return volume;
}

void Structure::set_voxel_volume(double vol) {
    this->voxel_volume = vol;
}

std::vector<double> Structure::calculate_dvh(bool relative) {
    if (masked_dose.size() == 0) return std::vector<double>();
    int max_dose = static_cast<int>(masked_dose[masked_dose.size() - 1].second) + 1;
    if (max_dose == 0) return std::vector<double>();
    // Want 1 dGy bins
    int num_bins = max_dose * 10;

    std::vector<int> l_ddvh;
    std::vector<double> l_cdvh;

    l_ddvh.resize(num_bins);
    std::fill(l_ddvh.begin(), l_ddvh.end(), 0.0);

    for (int i = 0; i < this->masked_dose.size(); i++) {
        int bin = this->masked_dose[i].second / 0.1;
        if (bin > l_ddvh.size() -1) std::cout << "WARNING" << std::endl;
        l_ddvh[bin] += 1;
    }

    l_cdvh.resize(num_bins);
    std::fill(l_cdvh.begin(), l_cdvh.end(), 0.0);

    double total = 0.0;
    for (int i = l_ddvh.size() - 1; i >= 0; i--) {
        total += l_ddvh[i];
        l_cdvh[i] = total * this->voxel_volume;
    }

    if (relative) {
        double norm_value = l_cdvh[0];
        if (norm_value > 0.0) {
            for (size_t i = 0; i < l_cdvh.size(); i++) {
                l_cdvh[i] = l_cdvh[i] / norm_value * 100;
            }
        }
    }

    this->cdvh = l_cdvh;
    this->ddvh = l_ddvh;

    return l_cdvh;

}

json Structure::get_dvh() {
    double prescription_dose = 0.0;
    for (int i = 0; i < this->hard_constraints.size(); i++) {
        if (this->hard_constraints[i].type == 0) prescription_dose = this->hard_constraints[i].threshold;
        break;
    }
    std::cout << this->name << std::endl;
    json dvh_dump;
    dvh_dump["name"] = this->name;

    if (prescription_dose > 0.0) {
        double d90 = this->dose_for_volume(90, true);
        double v200 = this->volume_for_dose(2.0 * prescription_dose, true);
        double v150 = this->volume_for_dose(1.5 * prescription_dose, true);
        double v100 = this->volume_for_dose(prescription_dose, true);
        std::cout << "D90: " << d90 << std::endl;
        std::cout << "V100: " << v100 << std::endl;
        std::cout << "V150: " << v150 << std::endl;
        std::cout << "V200: " << v200 << std::endl;
        dvh_dump["d90"] = d90;
        dvh_dump["v200"] = v200;
        dvh_dump["v150"] = v150;
        dvh_dump["v100"] = v100;
    }

    double d01cc = this->dose_for_volume(0.1, false);
    double d1cc = this->dose_for_volume(1.0, false);
    double d2cc = this->dose_for_volume(2.0, false);
    dvh_dump["d01cc"] = d01cc;
    dvh_dump["d1cc"] = d1cc;
    dvh_dump["d2cc"] = d2cc;
    std::cout << "D0.1 cm3: " << d01cc << std::endl;
    std::cout << "D1.0 cm3: " << d1cc << std::endl;
    std::cout << "D2.0 cm3: " << d2cc << std::endl;

    return dvh_dump;
}

json Structure::get_json_dvh() {
    std::vector<double> l_cdvh = this->calculate_dvh(true);
    json struct_dvh;
    struct_dvh["name"] = this->name;
    struct_dvh["dvh"] = std::vector<json>();
    for (int i = 0; i < l_cdvh.size(); i++) {
        json dose;
        dose["dose"] = i * 0.1;
        dose["vol"] = l_cdvh[i];
        struct_dvh["dvh"].push_back(dose);
    }
    return struct_dvh;
}
