// Copyright 2015 Marc-André Renaud
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "json.hh"

#include "string_utils.hh"
#include "BrachyOptimisation.hh"
#include "DwellPosition.hh"
#include "Phantom.hh"

using json = nlohmann::json;

BrachyOptimisation::BrachyOptimisation() {
}

BrachyOptimisation::BrachyOptimisation(json &input_json) {
    this->input_filename = input_json["name"];
    this->latest_cost = 0.0;
    this->iteration_number = 0;
    this->max_apertures = 0;
    this->from_json(input_json);
}

void BrachyOptimisation::from_json(json &input_json) {
    this->total_voxels = input_json["total_voxels"];
    this->lag_mults.resize(this->total_voxels);

    int num_dwells = input_json["dwells"].size();
    std::cout << "Number of dwells: " << num_dwells << std::endl;

    int dwell_counter = 0;
    for (auto& dwell_json : input_json["dwells"]) {
      std::cout << "Loading dwell: " << dwell_counter+1 << "/" <<
           num_dwells << "\r" << std::flush;

        DwellPosition *dwell = new DwellPosition(dwell_json);
        this->dwells.push_back(dwell);
        dwell_counter += 1;
    }

    if (input_json.find("max_apertures") != input_json.end()) {
        this->max_apertures = input_json["max_apertures"];
    }

    std::cout << "Number of structures: " << input_json["structures"].size() << std::endl;
    for (auto& struct_json : input_json["structures"]) {
        Structure *new_struct = new Structure(struct_json);
        this->structures.push_back(new_struct);
    }

    // Output constraints at the end of the input file reading.
    for (auto& structure: this->structures) {
        structure->output_constraints();
    }
}

void BrachyOptimisation::update_data() {
    this->iteration_number += 1;
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->active_weights.size();
    iter_data["cost_function_value"] = this->latest_cost;

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

void BrachyOptimisation::price_unused_dwells() {
    this->calculate_lag_mults();

    // For every control point not included in the plan yet,
    // price an aperture and select the control point with the best
    // aperture (highest price). The pricing routine requires knowledge
    // of the previous and next control point to account for machine
    // limitations.
    for (size_t i = 0; i < this->dwells.size(); i++) {
        if (this->dwells[i]->included == false) {
            std::cout << "pricing dwell " << i+1 << "/" <<
              this->dwells.size() << "\r" << std::flush;
            this->dwells[i]->price_dwell(this->lag_mults);
        }
    }
    std::cout << std::endl;
}

bool BrachyOptimisation::add_new_weight() {
    size_t max_index = 0;
    double max_price = 0.0;

    if (this->active_weights.size() > 0) {
        this->update_data();
        if (this->iteration_number % 5 == 0) {
            this->export_data();
        }
    }

    // If all available control points are included in the treatment plan,
    // then we can terminate the algorithm.
    if (this->max_apertures > 0 && this->active_weights.size() >= this->max_apertures) {
        return false;
    }

    if (this->active_weights.size() == this->dwells.size()) return 0;

    if (this->active_weights.size() > 0 && this->active_weights.size() % 10 == 0) {
      this->export_final_dose(true);
      this->export_final_weights(true);
    }

    this->price_unused_dwells();

    for (size_t i = 0; i < this->dwells.size(); i++) {
        if (this->dwells[i]->included == false && this->dwells[i]->latest_price > max_price) {
            max_index = i;
            max_price = this->dwells[i]->latest_price;
        }
    }

    if (max_price > 0.0) {
        this->dwells[max_index]->included = true;
        this->active_weights.push_back(this->dwells[max_index]);
        this->num_weights = this->active_weights.size();
        this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;
        return true;
    }

    // If max price is less than or equal to 0, then adding more control points
    // does not improve the dose distribution, we can terminate the algorithm.

    return false;
}

void BrachyOptimisation::export_final_dose(bool interm) {
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
    int *num_voxels = this->active_weights[0]->num_voxels;
    double *topleft = this->active_weights[0]->topleft;
    double *voxel_size = this->active_weights[0]->voxel_size;

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

void BrachyOptimisation::export_final_weights(bool interm) {
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
    json dwells_json;

    for (auto &roi : this->structures) {
        structs_json.push_back(roi->to_json());
    }

    for (auto &cpt : this->active_weights) {
        dwells_json.push_back(cpt->to_json());
    }

    plan_json["structures"] = structs_json;
    plan_json["dwells"] = dwells_json;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();

}
