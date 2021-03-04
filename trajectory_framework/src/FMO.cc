// Copyright 2015 Marc-André Renaud
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "json.hh"
#include "string_utils.hh"
#include "FMO.hh"
#include "Structure.hh"
#include "Aperture.hh"
#include "DenseDose.hh"
#include "ControlPoint.hh"

using boost::numeric::ublas::compressed_vector;
using json = nlohmann::json;

FMO::FMO() {
    this->output_modalities = false;
}
FMO::~FMO() {}


FMO::FMO(json &json_input) : FMO()
{
    std::cout << "Dense Fluence map optimisation!" << std::endl;
    this->input_filename = json_input["name"];
    this->from_json(json_input);
}

void FMO::from_json(json &json_input)
{
    this->output_modalities = json_input.value("output_modalities", false);

    int num_cpts = json_input["cpts"].size();
    std::cout << "Number of control points: " << num_cpts << std::endl;

    this->control_points = json_input["cpts"];
    this->load_beamlets(this->control_points);
    this->num_weights = this->beamlets.size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;

    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    for (auto &struct_json : json_input["structures"])
    {
        Structure *new_struct = new Structure(struct_json);
        this->structures.push_back(new_struct);
    }

    // Output constraints at the end of the input file reading.
    for (auto &structure : this->structures)
    {
        structure->output_constraints();
    }
}

void FMO::read_header(std::string filename)
{
    std::string extension;

    extension = find_extension(filename);
    if (extension == "3ddose")
    {
        this->read_3ddose_header(filename);
    }
    else if (extension == "bindos")
    {
        this->read_bindos_header(filename);
    }
    else
    {
        throw "Dose file extension not recognised";
    }
}

void FMO::read_bindos_header(std::string filename)
{
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    beamlet_file.read((char *)this->num_voxels, 3 * sizeof(int));

    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
              << this->num_voxels[1] << ","
              << this->num_voxels[2] << ")" << std::endl;

    float *x_voxels = new float[this->num_voxels[0] + 1];
    float *y_voxels = new float[this->num_voxels[1] + 1];
    float *z_voxels = new float[this->num_voxels[2] + 1];

    beamlet_file.read((char *)x_voxels, (this->num_voxels[0] + 1) * sizeof(float));
    beamlet_file.read((char *)y_voxels, (this->num_voxels[1] + 1) * sizeof(float));
    beamlet_file.read((char *)z_voxels, (this->num_voxels[2] + 1) * sizeof(float));

    beamlet_file.close();

    this->voxel_size[0] = x_voxels[1] - x_voxels[0];
    this->voxel_size[1] = y_voxels[1] - y_voxels[0];
    this->voxel_size[2] = z_voxels[1] - z_voxels[0];

    this->topleft[0] = x_voxels[0];
    this->topleft[1] = y_voxels[0];
    this->topleft[2] = z_voxels[0];

    std::cout << "Voxel size: (" << this->voxel_size[0] << ","
              << this->voxel_size[1] << ","
              << this->voxel_size[2] << ")" << std::endl;

    std::cout << "Topleft: (" << this->topleft[0] << ","
              << this->topleft[1] << ","
              << this->topleft[2] << ")" << std::endl;

    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    delete[] x_voxels;
    delete[] y_voxels;
    delete[] z_voxels;
}

void FMO::read_3ddose_header(std::string filename)
{
    /*
      Populates the class attributes related to phantom dimensions.
    */

    std::string buffer;
    std::string trimmed_buffer;
    std::vector<std::string> voxel_buffer;

    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer); // number of voxels
    std::stringstream convertor(buffer);
    convertor >> this->num_voxels[0] >> this->num_voxels[1] >> this->num_voxels[2];
    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
              << this->num_voxels[1] << ","
              << this->num_voxels[2] << ")" << std::endl;

    double second_voxel;
    std::getline(beamlet_file, buffer); // x voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[0];
    convertor >> second_voxel;
    std::cout << trimmed_buffer << std::endl;
    this->voxel_size[0] = second_voxel - this->topleft[0];

    std::getline(beamlet_file, buffer); // y voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[1];
    convertor >> second_voxel;
    this->voxel_size[1] = second_voxel - this->topleft[1];

    std::getline(beamlet_file, buffer); // z voxel coordinates
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

void FMO::load_beamlets(json &cpts)
{
    bool header_read = false;
    int cpt_index = 0;
    for (auto &cpt : cpts)
    {
        // peyman added cpt_index
        std::cout << cpt_index << std::endl;
        cpt_index++;
        for (std::string filename : cpt["beamlets"])
        {
            // peyman added this if block just to check
            if (filename != "Inactive")
            {
                DenseDose dose;
                dose.from_file(filename, this->total_voxels, false);
                this->beamlets.push_back(dose);

                if (!header_read) {
                    this->read_header(filename);
                    header_read = true;
                }
            }
        }
    }
}

std::vector<double> FMO::calculate_total_dose(const double *weights)
{
    std::vector<double> dose;
    dose.resize(this->total_voxels);
    dose.assign(this->total_voxels, 0.0);
    for (size_t beamlet_i = 0; beamlet_i < this->beamlets.size(); beamlet_i++)
    {
        for (size_t i = 0; i < this->total_voxels; i++)
        {
            dose[i] += this->beamlets[beamlet_i].grid[i] * weights[beamlet_i];
        }
    }

    return dose;
}

void FMO::refresh_doses(const double *weights)
{
    for (size_t st_i = 0; st_i < this->structures.size(); st_i++)
    {
        Structure *structure = this->structures[st_i];
        for (size_t i = 0; i < structure->masked_dose.size(); i++)
        {
            structure->masked_dose[i].second = 0.0;
            int vox_num = structure->masked_dose[i].first;
            for (size_t beamlet_i = 0; beamlet_i < this->beamlets.size(); beamlet_i++)
            {
                structure->masked_dose[i].second += this->beamlets[beamlet_i].grid[vox_num] *
                                                    weights[beamlet_i];
            }
        }
        std::sort(structure->masked_dose.begin(), structure->masked_dose.end(), [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
    }
}

// Total cost function is the sum of individual ROI cost functions.
double FMO::calculate_cost_function(const double *weights, bool new_weights)
{
    if (new_weights)
    {
        refresh_doses(weights);
    }

    double total_cost_function = 0.0;
    for (size_t i = 0; i < this->structures.size(); i++)
    {
        total_cost_function += this->structures[i]->calculate_cost();
    }

    this->latest_cost = total_cost_function;
    return total_cost_function;
}

void FMO::calculate_gradient(double *grad_values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        refresh_doses(weights);
    }

    for (size_t i = 0; i < this->num_weights; i++)
    {
        double derivative = 0.0;
        for (Structure *structure : this->structures)
        {
            derivative += structure->calculate_gradient(this->beamlets[i].grid);
        }
        grad_values[i] = derivative;
    }
}

void FMO::calculate_hessian(double *hess_values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        refresh_doses(weights);
    }

    int idx = 0;
    for (size_t row = 0; row < this->num_weights; row++)
    {
        for (size_t col = 0; col <= row; col++)
        {
            double hess = 0;
            for (size_t i = 0; i < this->structures.size(); i++)
            {
                hess += this->structures[i]->calculate_hessian(this->beamlets[row].grid,
                                                               this->beamlets[col].grid);
            }
            hess_values[idx] = hess;
            idx++;
        }
    }
}

void FMO::update_data() {
    // Cost data per iteration
    json iter_cost;
    iter_cost["structures"] = std::vector<nlohmann::json>();

    for (auto &structure : this->structures)
    {
        iter_cost["structures"].push_back(structure->update_data());
    }

    this->cost_data.push_back(iter_cost);
}

void FMO::export_data() {
    std::string filename3 = "combined_doses/" + remove_extension(this->input_filename) + ".cost";
    std::ofstream myfile3(filename3);
    myfile3 << nlohmann::json(this->cost_data).dump();
    myfile3.close();
}

void FMO::export_final_dose(const double *weights)
{
    
    std::vector<double> beamlet_dose;
    std::string buffer;
    std::string dose_filename;

    dose_filename = "combined_doses/" + split(input_filename, '.')[0] + ".3ddose";
    std::cout << "Exporting final dose to " << dose_filename << "\n";

    std::vector<double> final_dose = calculate_total_dose(weights);
    this->write_dose(dose_filename, final_dose);

    double sum_weights = 0.0;
    for (size_t dwell_i = 0; dwell_i < this->num_weights; dwell_i++)
    {
        sum_weights += weights[dwell_i];
    }

    std::cout << "Sum of weights: " << sum_weights << std::endl;

    this->export_final_weights(weights);
    this->update_data();
    this->export_data();

    if (this->output_modalities) {
        std::cout << "Outputting individual modalities" << std::endl;
        this->export_modality_doses(weights);
    }
}

void FMO::export_final_weights(const double *weights)
{
    std::string weights_filename = "combined_doses/" + remove_extension(this->input_filename) + ".weights";
    std::cout << "Exporting final weights to " << weights_filename << "\n";
    std::ofstream weights_file(weights_filename);

    json plan_json;
    json structs_json;
    json weights_json;

    for (auto roi : this->structures)
    {
        structs_json.push_back(roi->to_json());
    }

    for (size_t i = 0; i < this->num_weights; i++)
    {
        weights_json.push_back(weights[i]);
    }


    size_t idx = 0;
    for (auto &cpt : this->control_points)
    {
        std::vector<double> cpt_weights;
        for (size_t b_i = 0; b_i < cpt["beamlets"].size(); b_i++)
        {
            if (cpt["beamlets"][b_i] == "Inactive")
            {
                cpt_weights.push_back(0.0);
            }
            else
            {
                cpt_weights.push_back(weights[idx]);
                idx += 1;
            }
        }
        cpt["beamlet_weights"] = cpt_weights;
    }

    // size_t idx = 0;
    // for (auto &cpt : this->control_points)
    // {
    //     std::vector<double> cpt_weights;
    //     for (size_t b_i = 0; b_i < cpt["beamlets"].size(); b_i++) {
    //         cpt_weights.push_back(weights[idx]);
    //         idx += 1;
    //     }
    //     cpt["beamlet_weights"] = cpt_weights;
    // }
    plan_json["control_points"] = this->control_points;

    plan_json["structures"] = structs_json;
    plan_json["weights"] = weights_json;
    plan_json["cost_function_value"] = this->latest_cost;

    weights_file << plan_json.dump();

    weights_file.close();
}

std::vector<double> FMO::calculate_modality_dose(std::vector<int> apertures, const double *weights) {
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto ap_i : apertures)
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += this->beamlets[ap_i].grid[i] * weights[ap_i];
        }
    }

    return total_dose;
}

void FMO::export_modality_doses(const double *weights)
{
    std::vector<int> electron_energies;
    std::vector<int> photon_energies;
    std::map<int, std::vector<int> > electron_apertures;
    std::map<int, std::vector<int> > photon_apertures;

    std::vector<int> all_electrons;
    std::vector<int> all_photons;

    int idx = 0;
    for (int i = 0; i < this->control_points.size(); i++) {
        auto cpt = this->control_points[i];

        std::string particle = cpt.value("particle", "electron");
        int energy = cpt["energy"];

        if (particle == "electron") {
            if (electron_apertures.count(energy) == 0)
            {
                electron_apertures[energy] = std::vector<int>();
                electron_energies.push_back(energy);
            }

            for (int b_i = 0; b_i < cpt["beamlets"].size(); b_i++) {
                electron_apertures[energy].push_back(idx);
                all_electrons.push_back(idx);
                idx += 1;
            }
        } else {
            if (photon_apertures.count(energy) == 0)
            {
                photon_apertures[energy] = std::vector<int>();
                photon_energies.push_back(energy);
            }

            for (int b_i = 0; b_i < cpt["beamlets"].size(); b_i++) {
                photon_apertures[energy].push_back(idx);
                all_photons.push_back(idx);
                idx += 1;
            }
        }
    }

    for (auto energy : electron_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MeV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(electron_apertures[energy], weights);
        this->write_dose(dose_filename, modality_dose);
    }

    for (auto energy : photon_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(photon_apertures[energy], weights);
        this->write_dose(dose_filename, modality_dose);
    }

    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_electrons.3ddose";
        std::cout << "Exporting electron dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(all_electrons, weights);
        this->write_dose(dose_filename, modality_dose);
    }

    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_photons.3ddose";
        std::cout << "Exporting photon dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(all_photons, weights);
        this->write_dose(dose_filename, modality_dose);
    }
}

void FMO::write_dose(const std::string filename, const std::vector<double> dose) {
    std::ofstream dose_file(filename);

    dose_file << "  " << num_voxels[0] << "  " << num_voxels[1]
                << "  " << num_voxels[2] << "\n";

    for (int i = 0; i < num_voxels[0] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * voxel_size[0] + topleft[0];
    }
    dose_file << "\n";

    for (int i = 0; i < num_voxels[1] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * voxel_size[1] + topleft[1];
    }
    dose_file << "\n";

    for (int i = 0; i < num_voxels[2] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * voxel_size[2] + topleft[2];
    }
    dose_file << "\n";

    for (size_t i = 0; i < dose.size(); i++)
    {
        if (i != 0)
            dose_file << " ";
        if (dose[i] > 0)
        {
            dose_file << dose[i];
        }
        else
        {
            dose_file << "0";
        }
    }
    dose_file << "\n";

    for (size_t i = 0; i < dose.size(); i++)
    {
        dose_file << "0 ";
    }
    dose_file << "\n";

    dose_file.close();
}
