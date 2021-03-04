// Copyright 2017 Marc-André Renaud
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "json.hh"
#include "string_utils.hh"
#include "FluenceMapOptimisation.hh"
#include "FMStructure.hh"
#include "Aperture.hh"
#include "DenseDose.hh"
#include "SparseDose.hh"

using boost::numeric::ublas::compressed_vector;
using json = nlohmann::json;

FluenceMapOptimisation::FluenceMapOptimisation(json &json_input) {
    std::cout << "Fluence map optimisation!" << std::endl;
    this->input_filename = json_input["name"];
    this->from_json(json_input);
}

FluenceMapOptimisation::FluenceMapOptimisation() {}
FluenceMapOptimisation::~FluenceMapOptimisation() {}

void FluenceMapOptimisation::from_json(json &json_input) {
    this->beamlet_filenames = json_input["filenames"].get<std::vector<std::string> >();;

    this->num_weights = this->beamlet_filenames.size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;

    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    for (auto& struct_json : json_input["structures"]) {
        FMStructure new_struct(struct_json);
        this->structures.push_back(new_struct);
    }

    // Output constraints at the end of the input file reading.
    for (auto& structure: this->structures) {
        structure.output_constraints();
    }


    this->read_header(this->beamlet_filenames[0]);
    this->load_structure_beamlets(this->beamlet_filenames);
}

void FluenceMapOptimisation::read_header(std::string filename) {
    std::string extension;
    std::vector<std::string> buffer;

    extension = find_extension(filename);
    if (extension == "3ddose") {
        this->read_3ddose_header(filename);
    } else if(extension == "minidos") {
        this->read_binary_header(filename);
    } else {
        throw "Dose file extension not recognised";
    }
}

void FluenceMapOptimisation::read_3ddose_header(std::string filename) {
    /*
      Populates the class attributes related to phantom dimensions.
    */

    std::string buffer;
    std::string trimmed_buffer;
    std::vector<std::string> voxel_buffer;

    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer);  // number of voxels
    std::stringstream convertor(buffer);
    convertor >> this->num_voxels[0] >> this->num_voxels[1] >> this->num_voxels[2];
    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
                                       << this->num_voxels[1] << ","
                                       << this->num_voxels[2] << ")" << std::endl;

    double second_voxel;
    std::getline(beamlet_file, buffer);  // x voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[0];
    convertor >> second_voxel;
    std::cout << trimmed_buffer << std::endl;
    this->voxel_size[0] = second_voxel - this->topleft[0];

    std::getline(beamlet_file, buffer);  // y voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[1];
    convertor >> second_voxel;
    this->voxel_size[1] = second_voxel - this->topleft[1];

    std::getline(beamlet_file, buffer);  // z voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[2];
    convertor >> second_voxel;
    this->voxel_size[2] = second_voxel - this->topleft[2];

    std::cout << "Voxel size: (" << this->voxel_size[0] << ","
                                 << this->voxel_size[1] << ","
                                 << this->voxel_size[2] << ")" << std::endl;

    std::cout << "Topleft: (" << this->topleft[0] << ","
                                 << this->topleft[1] << ","
                                 << this->topleft[2] << ")" << std::endl;

    beamlet_file.close();

}

void FluenceMapOptimisation::read_binary_header(std::string filename) {
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    float b_topleft[3], b_voxel_size[3];
    beamlet_file.read((char*)this->num_voxels, 3 * sizeof(int));
    beamlet_file.read((char*)b_voxel_size, 3 * sizeof(float));
    beamlet_file.read((char*)b_topleft, 3 * sizeof(float));
    this->voxel_size[0] = b_voxel_size[0];
    this->voxel_size[1] = b_voxel_size[1];
    this->voxel_size[2] = b_voxel_size[2];

    this->topleft[0] = b_topleft[0];
    this->topleft[1] = b_topleft[1];
    this->topleft[2] = b_topleft[2];

    this->total_voxels = this->num_voxels[0]
                         * this->num_voxels[1]
                         * this->num_voxels[2];
    beamlet_file.close();
    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
                                       << this->num_voxels[1] << ","
                                       << this->num_voxels[2] << ")" << std::endl;


    std::cout << "Voxel size: (" << this->voxel_size[0] << ","
                                 << this->voxel_size[1] << ","
                                 << this->voxel_size[2] << ")" << std::endl;

    std::cout << "Topleft: (" << this->topleft[0] << ","
                                 << this->topleft[1] << ","
                                 << this->topleft[2] << ")" << std::endl;


}

void FluenceMapOptimisation::load_structure_beamlets(std::vector<std::string> beamlet_filenames) {
    size_t num_beamlets = beamlet_filenames.size();
    for (size_t i = 0; i < num_beamlets; i++) {
        std::cout << "loading beamlet " << i+1 << "/" << num_beamlets << "\r" << std::flush;
        SparseDose dose;
        dose.from_file(beamlet_filenames[i], true);
        for (FMStructure &structure: this->structures) {
            std::vector<double> struct_beamlet_dose;
            struct_beamlet_dose.resize(structure.mask.size());
            for (size_t i = 0; i < structure.mask.size(); i++) {
                struct_beamlet_dose[i] = dose.grid[structure.mask[i]];
            }
            structure.beamlets.push_back(struct_beamlet_dose);
        }
    }
    std::cout << std::endl;
}

void FluenceMapOptimisation::load_beamlets(std::vector<std::string> beamlet_filenames) {
    size_t num_beamlets = beamlet_filenames.size();
    for (size_t i = 0; i < num_beamlets; i++) {
        std::cout << "loading beamlet " << i+1 << "/" << num_beamlets << "\r" << std::flush;
        SparseDose dose;
        dose.from_file(beamlet_filenames[i], true);
        this->beamlets.push_back(dose);
    }

    std::cout << std::endl;
}

std::vector<double> FluenceMapOptimisation::calculate_total_dose(const double* weights) {
    std::vector<double> dose;
    dose.resize(this->total_voxels);
    std::fill(dose.begin(), dose.end(), 0.0);

    size_t num_beamlets = this->beamlet_filenames.size();
    for (size_t bbi = 0; bbi < num_beamlets; bbi++) {
        std::cout << "loading beamlet " << bbi+1 << "/" << num_beamlets << "\r" << std::flush;
        DenseDose beamlet_dose;
        beamlet_dose.from_file(this->beamlet_filenames[bbi], this->total_voxels, true);
        for (size_t i = 0; i < this->total_voxels; i++) {
            dose[i] += beamlet_dose.grid[i] * weights[bbi];
        }
    }

    return dose;
}

void FluenceMapOptimisation::refresh_doses(const double* weights) {
    for (FMStructure &structure: this->structures) {
        std::fill(structure.masked_dose.begin(), structure.masked_dose.end(), 0.0);
        for (size_t bbi = 0; bbi < structure.beamlets.size(); bbi++) {
            for (size_t i = 0; i < structure.masked_dose.size(); i++) {
                structure.masked_dose[i] += structure.beamlets[bbi][i] * weights[bbi];
            }
        }
    }
}

// Total cost function is the sum of individual ROI cost functions.
double FluenceMapOptimisation::calculate_cost_function(const double *weights, bool new_weights) {
    if (new_weights) {
        refresh_doses(weights);
    }

    double total_cost_function = 0.0;
    for (size_t i = 0; i < this->structures.size(); i++) {
        total_cost_function += this->structures[i].calculate_cost();
    }

    this->latest_cost = total_cost_function;
    return total_cost_function;
}


void FluenceMapOptimisation::calculate_gradient(double *grad_values, const double *weights, bool new_weights) {
    if (new_weights) {
        refresh_doses(weights);
    }

    for (size_t i = 0; i < this->num_weights; i++) {
        double derivative = 0.0;
        for (FMStructure &structure: this->structures) {
            derivative += structure.calculate_gradient(structure.beamlets[i]);
        }
        grad_values[i] = derivative;
    }
}

void FluenceMapOptimisation::calculate_hessian(double *hess_values, const double* weights, bool new_weights) {
    if (new_weights) {
      refresh_doses(weights);
    }

    std::fill_n(hess_values, num_hess_ele, 0.0);

    int idx = 0;
    for (size_t row = 0; row < this->num_weights; row++) {
        for (size_t col = 0; col <= row; col++) {
            for (size_t i = 0; i < this->structures.size(); i++) {
                hess_values[idx] += this->structures[i].calculate_hessian(this->structures[i].beamlets[row],
                                                                              this->structures[i].beamlets[col]);
            }
            idx++;
        }
    }
}

void FluenceMapOptimisation::export_final_dose(const double *weights) {
    std::string buffer;
    std::string dose_filename;

    dose_filename = "combined_doses/" + split(input_filename, '.')[0] + ".3ddose";
    std::cout << "Exporting final dose to " << dose_filename << "\n";
    std::ofstream dose_file(dose_filename);

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

    std::vector<double> final_dose = calculate_total_dose(weights);
    for (size_t i = 0; i < total_voxels; i++) {
      if (i != 0) dose_file << " ";
      if (final_dose[i] > 0) {
        dose_file << final_dose[i];
      } else {
        dose_file << "0";
      }
    }
    dose_file << "\n";

    // No uncertainties yet.
    for (size_t i = 0; i < total_voxels; i++) {
      if (i != 0) dose_file << " ";
      dose_file << 0.0;
    }
    dose_file << "\n";

    dose_file.close();

    this->export_final_weights(weights);
}

void FluenceMapOptimisation::export_final_weights(const double *weights) {
    std::string buffer;
    std::string plan_filename;

    plan_filename = "combined_doses/" + remove_extension(this->input_filename) + ".weights";

    std::cout << "Exporting final weights to " << plan_filename << "\n";
    std::ofstream plan_file(plan_filename);

    json plan_json;
    json structs_json;

    std::vector<double> vec_weights;
    for (size_t i = 0; i < this->num_weights; i++) {
        vec_weights.push_back(weights[i]);
    }

    for (auto &roi : this->structures) {
        structs_json.push_back(roi.to_json());
    }

    plan_json["structures"] = structs_json;
    plan_json["weights"] = vec_weights;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();
}

