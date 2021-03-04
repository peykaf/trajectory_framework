// Copyright 2015 Marc-André Renaud
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "json.hh"

#include "string_utils.hh"
#include "KVATOptimisation.hh"
#include "KVATControlPoint.hh"
#include "Phantom.hh"
#include "KVATBeamlet.hh"

using json = nlohmann::json;

KVATOptimisation::KVATOptimisation() {}

KVATOptimisation::KVATOptimisation(json &json_input) {
    std::cout << "Mixed optimisation!" << std::endl;
    this->input_filename = json_input["name"];
    this->max_apertures = 0;
    this->from_json(json_input);
    this->latest_cost = 0.0;
    this->iteration_number = 0;
}

void KVATOptimisation::from_json(json &json_input) {
    this->total_voxels = json_input["total_voxels"];
    this->lag_mults.resize(this->total_voxels);

    int num_cpts = json_input["cpts"].size();
    std::cout << "Number of control points: " << num_cpts << std::endl;

    for (auto& cpt_json : json_input["cpts"]) {
        KVATControlPoint *cpt = new KVATControlPoint(cpt_json);
        this->cpts.push_back(cpt);
    }

    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    for (auto& struct_json : json_input["structures"]) {
        Structure *new_struct = new Structure(struct_json);
        this->structures.push_back(new_struct);
    }

    // Output constraints at the end of the input file reading.
    for (auto& structure: this->structures) {
        structure->output_constraints();
    }

    if (json_input.find("max_apertures") != json_input.end()) {
        this->max_apertures = json_input["max_apertures"];
    }
}

void KVATOptimisation::price_unused_cpts() {
    // For every control point not included in the plan yet,
    // price an aperture and select the control point with the best
    // aperture (highest price). The pricing routine requires knowledge
    // of the previous and next control point to account for machine
    // limitations.
    for (size_t i = 0; i < this->cpts.size(); i++) {
        this->cpts[i]->price_beamlets(this->lag_mults, this->structures);
    }
    std::cout << std::endl;
}

void KVATOptimisation::prune_cpts() {
    double avg_weight = 0.0;
    int num_zero_weights = 0;
    for (size_t i = 0; i < this->active_weights.size(); i++) {
        avg_weight += this->active_weights[i]->weight;
    }

    avg_weight /= this->active_weights.size();

    for (size_t i = 0; i < this->active_weights.size(); i++) {
        if (this->active_weights[i]->weight < avg_weight / 1000.0) {
            num_zero_weights += 1;
        }
    }

    if (num_zero_weights > 3) {
        for (size_t i = 0; i < this->active_weights.size(); i++) {
            if (this->active_weights[i]->weight < avg_weight / 1000.0) {
                this->active_weights[i]->deactivate();
                this->active_weights.erase(this->active_weights.begin() + i);
            }
        }
    }
}

void KVATOptimisation::update_data() {
    this->iteration_number += 1;
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->active_weights.size();
    iter_data["cost_function_value"] = this->latest_cost;

    json p_cpt_stats;
    for (size_t i = 0; i < this->cpts.size(); i++) {
        p_cpt_stats.push_back(this->cpts[i]->write_statistics());
    }

    iter_data["photon_cpt_stats"] = p_cpt_stats;

    this->plan_data[this->iteration_number] = iter_data;

    // Cost data per iteration
    json iter_cost;
    iter_cost["cost_function_value"] = this->latest_cost;
    iter_cost["structures"] = std::vector<nlohmann::json>();

    for (auto &structure : this->structures) {
        iter_cost["structures"].push_back(structure->update_data());
    }

    this->cost_data.push_back(iter_cost);

    // DVH data per iteration
    json iter_dvh;
    iter_dvh["structures"] = std::vector<nlohmann::json>();
    iter_dvh["dose_step"] = 0.1;
    for (auto &structure : this->structures) {
        iter_dvh["structures"].push_back(structure->get_json_dvh());
    }
    this->dvh_data.push_back(iter_dvh);
}

void KVATOptimisation::output_plan_statistics() {
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << this->cpts.size() << " photon control points" << std::endl;
    for (auto &cpt : this->cpts) {
        cpt->output_statistics();
    }
}

bool KVATOptimisation::add_new_weight() {
    size_t max_index = -1;
    double max_price = 0.0;

    this->refresh_doses();

    if (this->active_weights.size() > 0) {
        this->update_data();
        if (this->iteration_number % 5 == 0) {
            this->export_data();
        }
    }

    this->output_cost();

    if (this->max_apertures > 0 && this->active_weights.size() >= this->max_apertures) {
        return false;
    }

    // Export intermediate results every 10 control points.
    if (this->active_weights.size() > 0 && this->active_weights.size() % 10 == 0) {
      this->export_final_dose(true);
      this->export_final_weights(true);
      this->output_plan_statistics();
    }

    this->calculate_lag_mults();
    this->price_unused_cpts();

    bool added_cpt = false;

    for (size_t i = 0; i < this->cpts.size(); i++) {
        if (this->cpts[i]->latest_price > max_price) {
            max_index = i;
            max_price = this->cpts[i]->latest_price;
        }
    }

    std::cout << "Max photon price: " << max_price << std::endl;

    if (max_price > 0.0) {
        this->add_beamlet(max_index);
        added_cpt = true;
    }

    // If max price is less than or equal to 0, then adding more control points
    // does not improve the dose distribution, we can terminate the algorithm.

    return added_cpt;
}

void KVATOptimisation::add_beamlet(int index) {
    this->cpts[index]->included = true;
    KVATBeamlet *ap_ptr = this->cpts[index]->add_best_beamlet();

    this->active_weights.push_back(ap_ptr);
    this->num_weights = this->active_weights.size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;
}

void KVATOptimisation::export_final_dose(bool interm) {
    std::string buffer;
    std::string dose_filename;

    if (interm) {
      dose_filename = "combined_doses/" + remove_extension(input_filename) + "_interm.3ddose";
    } else {
      dose_filename = "combined_doses/" + remove_extension(input_filename) + ".3ddose";
    }

    std::cout << "Exporting final dose to " << dose_filename << "\n";
    std::ofstream dose_file(dose_filename);

    size_t total_voxels = this->total_voxels;

    int *num_voxels;
    double *topleft;
    double *voxel_size;

    num_voxels = this->cpts[0]->num_voxels;
    topleft = this->cpts[0]->topleft;
    voxel_size = this->cpts[0]->voxel_size;

    dose_file << "  " << num_voxels[0] << "  " << num_voxels[1]
      << "  " << num_voxels[2] << "\n";

    for (int i = 0; i < num_voxels[0] + 1; i++) {
      if (i != 0) dose_file << " ";
      dose_file << i * voxel_size[0] + topleft[0];
    }
    dose_file << "\n";
    for (int i = 0; i < num_voxels[1] + 1; i++) {
      if (i != 0) dose_file << " ";
      dose_file << i * voxel_size[1] + topleft[1];
    }
    dose_file << "\n";
    for (int i = 0; i < num_voxels[2] + 1; i++) {
      if (i != 0) dose_file << " ";
      dose_file << i * voxel_size[2] + topleft[2];
    }
    dose_file << "\n";

    std::vector<double> total_dose = this->calculate_total_dose();
    for (size_t i = 0; i < total_voxels; i++) {
      if (i != 0) dose_file << " ";
      if (total_dose[i] > 0) {
        dose_file << total_dose[i];
      } else {
        dose_file << "0";
      }
    }
    dose_file << "\n";

    for (size_t i = 0; i < total_voxels; i++)
        dose_file << "0 ";
    dose_file << "\n";

    dose_file.close();
}

void KVATOptimisation::export_final_weights(bool interm) {
    std::string buffer;
    std::string plan_filename;

    if (interm) {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + "_interm.weights";
    } else {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + ".weights";
    }

    std::cout << "Exporting final weights to " << plan_filename << "\n";
    std::ofstream plan_file(plan_filename);

    json plan_json;
    json structs_json;
    json photon_cpts_json;

    for (auto &roi : this->structures) {
        structs_json.push_back(roi->to_json());
    }

    for (auto &cpt : this->cpts) {
        if (cpt->included) {
            photon_cpts_json.push_back(cpt->to_json());
        }
    }

    plan_json["structures"] = structs_json;
    plan_json["cpts"] = photon_cpts_json;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();
}
