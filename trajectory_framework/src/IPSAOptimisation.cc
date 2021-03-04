// Copyright 2015 Marc-André Renaud
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>

#include "json.hh"
#include "string_utils.hh"
#include "Structure.hh"
#include "DwellPosition.hh"
#include "IPSAOptimisation.hh"
#include "Phantom.hh"

using json = nlohmann::json;

IPSAOptimisation::IPSAOptimisation() : rnd_gen(876423178)
{
    std::cout << "Brachytherapy IPSA" << std::endl;
    this->latest_cost = 0.0;
    this->iteration_number = 1;
    this->max_iterations = 100000000;
    this->norm_constant = 1.0;
    this->latest_cost = std::numeric_limits<double>::max();
    this->original_cost = std::numeric_limits<double>::max();

    // SA cooling constants
    this->T_0 = 0.0;   // To be initialized
    this->alpha = 0.6; // Based on Lessard/Pouliot paper (2001)
}

IPSAOptimisation::IPSAOptimisation(json &json_input) : IPSAOptimisation()
{
    this->input_filename = json_input["name"];
    this->from_json(json_input);
    this->dwell_distribution = std::uniform_int_distribution<size_t>(0, this->dwell_positions.size()-1);
    this->real_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
}

void IPSAOptimisation::from_json(json &json_input)
{
    this->total_voxels = json_input["total_voxels"];

    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    for (auto &struct_json : json_input["structures"])
    {
        Structure *new_struct = new Structure(struct_json);
        this->structures.push_back(new_struct);
    }

    // Output constraints at the end of the input file reading.
    double min_dose = 0.0;
    for (auto &structure : this->structures)
    {
        structure->output_constraints();
        double struct_min_dose = structure->get_min_dose();
        min_dose = std::max(min_dose, struct_min_dose);
    }
    this->norm_constant = std::max(1.0, min_dose);

    int num_cpts = json_input["dwells"].size();
    std::cout << "Number of dwell positions: " << num_cpts << std::endl;

    int dwell_counter = 0;
    int num_dwells = json_input["dwells"].size();
    for (auto &dwell_json : json_input["dwells"])
    {
        std::cout << "Loading dwell: " << dwell_counter+1 << "/" <<
           num_dwells << "\r" << std::flush;

        DwellPosition *dwell = new DwellPosition(dwell_json);
        this->dwell_positions.push_back(dwell);
        this->active_weights.push_back(dwell);
        dwell_counter += 1;
    }


    this->start_time = std::chrono::steady_clock::now();
}

void IPSAOptimisation::initialize()
{
    size_t num_init_iterations = 500;
    double avg_cost_diff = 0.0;

    for (size_t i = 0; i < num_init_iterations; i++)
    {
        double cost_diff = this->do_iteration();
        avg_cost_diff += std::abs(cost_diff);
    }

    this->T_0 = avg_cost_diff / num_init_iterations;
    this->T_0 = 0.0;
}

void IPSAOptimisation::launch()
{
    this->refresh_doses();
    this->latest_cost = this->calculate_cost();
    this->original_cost = this->latest_cost;
    this->lowest_cost = this->latest_cost;
    this->initialize();

    while (this->iteration_number < this->max_iterations)
    {
        if (this->iteration_number % 1000 == 0 && this->iteration_number % 50000 != 0)
        {
            this->output_cost(false);
            //this->update_data();
            //this->export_data();
        }
        else if (this->iteration_number % 50000 == 0)
        {
            this->output_cost(true);
            this->export_dose(true);
            this->export_weights(true);

            // Get rid of floating point error accumulation
            this->refresh_doses();
        }

        this->do_iteration();
        this->iteration_number += 1;
        if (this->latest_cost < this->lowest_cost)
            this->lowest_cost = this->latest_cost;
    }

    this->export_dose(false);
    this->export_weights(false);
}

double IPSAOptimisation::do_iteration()
{
    size_t dwell_index = this->dwell_distribution(this->rnd_gen);
    DwellPosition *dwell = this->dwell_positions[dwell_index];

    double old_weight = dwell->weight;
    double dwell_diff_sigma = std::max(0.001, this->dwell_positions.size() / std::sqrt(this->iteration_number)) * norm_constant;
    std::normal_distribution<double> norm_dist = std::normal_distribution<double>(0.0, dwell_diff_sigma);

    double new_weight = dwell->weight + norm_dist(this->rnd_gen);
    if (new_weight < 0.0) {
        new_weight = 0.0;
    }

    dwell->weight = new_weight;
    this->refresh_doses(dwell, old_weight);

    double new_cost = this->calculate_cost();
    double cost_diff = new_cost - this->latest_cost;

    if (new_cost >= this->latest_cost)
    {
        // New distribution is worse, accept with some probability
        double temp = this->T_0 / std::pow(this->iteration_number, this->alpha);
        double p_accept = 0.0;
        if (temp > 1e-10) {
            p_accept = std::exp(-cost_diff / temp);
        }

        double rndnum = this->real_distribution(this->rnd_gen);
        if (rndnum < p_accept)
        {
            this->latest_cost = new_cost;
        }
        else
        {
            // Reject this solution
            dwell->weight = old_weight;
            // Reset back to old weight
            this->refresh_doses(dwell, new_weight);
            cost_diff = 0.0;
        }
    }
    else
    {
        this->latest_cost = new_cost;
    }

    return cost_diff;
}

void IPSAOptimisation::update_structure_dose(Structure *structure, DwellPosition *dwell, double old_weight) {
    double weight_diff = dwell->weight - old_weight;
    for (auto &vox : structure->masked_dose) {
        size_t vox_num = vox.first;
        vox.second = vox.second + dwell->dose[vox_num] * weight_diff;
    }

    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
}

void IPSAOptimisation::calculate_structure_dose(Structure *structure)
{
    for (auto &vox : structure->masked_dose) {
        vox.second = 0.0;
        size_t vox_num = vox.first;
        for (const auto &dwell : this->active_weights) {
            vox.second += dwell->dose[vox_num] * dwell->weight;
        }
    }

    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
}

void IPSAOptimisation::refresh_doses()
{
    for (auto &structure : this->structures)
    {
        this->calculate_structure_dose(structure);
    }
}

void IPSAOptimisation::refresh_doses(DwellPosition *dwell, double old_weight)
{
    for (auto &structure : this->structures)
    {
        this->update_structure_dose(structure, dwell, old_weight);
    }
}

double IPSAOptimisation::calculate_cost()
{
    double cost = 0.0;
    for (auto &structure : this->structures)
    {
        cost += structure->calculate_cost();
    }

    return cost;
}

void IPSAOptimisation::update_data()
{
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->active_weights.size();
    iter_data["cost_function_value"] = this->latest_cost;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - this->start_time).count();

    iter_data["running_time"] = elapsed;

    json dwell_stats;
    for (size_t i = 0; i < this->dwell_positions.size(); i++)
    {
        dwell_stats.push_back(this->dwell_positions[i]->write_statistics());
    }

    iter_data["dwell_stats"] = dwell_stats;

    this->plan_data[this->iteration_number] = iter_data;
}

void IPSAOptimisation::output_plan_statistics()
{
    std::cout << this->dwell_positions.size() << " dwell positions" << std::endl;
    for (auto &dwell : this->dwell_positions)
    {
        dwell->output_statistics();
    }
}

size_t IPSAOptimisation::calculate_active_dwells() {
    size_t active_dwells = 0;
    for (auto &dwell : this->dwell_positions)
    {
        if (dwell->weight > 0.0) active_dwells += 1;
    }

    return active_dwells;
}

void IPSAOptimisation::output_cost(bool verbose)
{
    std::cout << "Original cost function value: " << this->original_cost << std::endl;
    std::cout << "Original temperature value: " << this->T_0 << std::endl;
    std::cout << "Lowest cost function value: " << this->lowest_cost << std::endl;
    std::cout << "Latest cost function value: " << this->latest_cost << std::endl;
    std::cout << "Number of active dwells: " << this->calculate_active_dwells() << std::endl;
    std::cout << "Dwell time adjustments: " << std::max(0.001, 1.0 / std::sqrt(this->iteration_number)) * this->norm_constant << std::endl;
    std::cout << "Iteration " << this->iteration_number << "/" << this->max_iterations << std::endl;

    double temp = this->T_0 / std::pow(static_cast<double>(this->iteration_number), this->alpha);
    double p_accept = 0.0;
    if (temp > 1e-10) {
        p_accept = std::exp(-this->T_0 / temp);
    }
    std::cout << "Estimated probability of accepting a worse solution: " << p_accept << std::endl;

    if (verbose)
    {
        for (size_t i = 0; i < this->structures.size(); i++)
        {
            this->structures[i]->output_cost();
        }
    }
}

std::vector<double> IPSAOptimisation::calculate_total_dose()
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &ap : this->active_weights)
    {
        if (ap->weight > 0.0) {
            for (size_t i = 0; i < total_dose.size(); i++)
            {
                total_dose[i] += ap->dose[i] * ap->weight;
            }
        }

    }

    return total_dose;
}

void IPSAOptimisation::export_data()
{
    std::string filename = "combined_doses/" + remove_extension(this->input_filename) + ".data";
    std::ofstream myfile(filename);
    myfile << this->plan_data.dump();
    myfile.close();

    std::string filename2 = "combined_doses/" + remove_extension(this->input_filename) + ".progress";
    std::ofstream myfile2(filename2);
    myfile2 << this->iteration_number << "/" << this->max_iterations << std::endl;
    myfile2.close();
}

void IPSAOptimisation::export_dose(bool interm)
{
    std::string buffer;
    std::string dose_filename;

    if (interm)
    {
        dose_filename = "combined_doses/" + remove_extension(input_filename) + "_interm.3ddose";
    }
    else
    {
        dose_filename = "combined_doses/" + remove_extension(input_filename) + ".3ddose";
    }

    std::cout << "Exporting dose to " << dose_filename << "\n";
    std::ofstream dose_file(dose_filename);

    int *num_voxels;
    double *topleft;
    double *voxel_size;

    num_voxels = this->dwell_positions[0]->num_voxels;
    topleft = this->dwell_positions[0]->topleft;
    voxel_size = this->dwell_positions[0]->voxel_size;

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

    std::vector<double> total_dose = this->calculate_total_dose();
    for (size_t i = 0; i < total_dose.size(); i++)
    {
        if (i != 0)
            dose_file << " ";

        if (total_dose[i] > 0)
        {
            dose_file << total_dose[i];
        }
        else
        {
            dose_file << "0";
        }
    }
    dose_file << "\n";

    for (size_t i = 0; i < total_dose.size(); i++)
    {
        dose_file << "0 ";
    }

    dose_file << "\n";

    dose_file.close();
}

void IPSAOptimisation::export_weights(bool interm)
{
    std::string buffer;
    std::string plan_filename;

    if (interm)
    {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + "_interm.weights";
    }
    else
    {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + ".weights";
    }

    std::cout << "Exporting final weights to " << plan_filename << "\n";
    std::ofstream plan_file(plan_filename);

    json plan_json;
    json structs_json;
    json dwells_json;
    json photon_cpts_json;

    for (auto &roi : this->structures)
    {
        structs_json.push_back(roi->to_json());
    }

    for (auto &cpt : this->dwell_positions)
    {
        dwells_json.push_back(cpt->to_json());
    }

    plan_json["structures"] = structs_json;
    plan_json["dwell_positions"] = dwells_json;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();
}
