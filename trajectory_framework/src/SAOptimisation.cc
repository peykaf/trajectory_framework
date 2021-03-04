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
#include "MixedOptimisation.hh"
#include "MixedArcOptimisation.hh"
#include "SAOptimisation.hh"
#include "ControlPoint.hh"
#include "Phantom.hh"

using json = nlohmann::json;

SAOptimisation::SAOptimisation(MixedOptimisation *mixed_opt) : rnd_gen(100)
{
    this->input_filename = mixed_opt->input_filename + "_SA";
    this->total_voxels = mixed_opt->get_total_voxels();
    this->structures = mixed_opt->get_structures();
    this->photon_cpts = mixed_opt->get_photon_cpts();
    this->electron_cpts = mixed_opt->get_electron_cpts();
    this->N_L = 0;
    this->constrained = false;

    this->initialize();
}

SAOptimisation::SAOptimisation(MixedArcOptimisation *mixed_opt) : rnd_gen(100)
{
    this->input_filename = mixed_opt->input_filename + "_SA";
    this->total_voxels = mixed_opt->get_total_voxels();
    this->structures = mixed_opt->get_structures();
    this->photon_cpts = mixed_opt->get_photon_cpts();
    this->electron_cpts = mixed_opt->get_electron_cpts();
    this->N_L = 0;
    this->constrained = true;

    this->initialize();
}

SAOptimisation::SAOptimisation(json &json_input)
{
    std::cout << "Simulated Annealing" << std::endl;
    this->input_filename = json_input["name"];
    this->max_iterations = 1000000;
    this->from_json(json_input);
    this->latest_cost = 0.0;
    this->iteration_number = 0;
}

void SAOptimisation::initialize()
{
    for (auto &cpt : this->photon_cpts)
    {
        cpt->remove_inactive_apertures();
        this->all_cpts.push_back(cpt);

        for (auto &ap : cpt->apertures)
        {
            this->active_weights.push_back(ap);
        }

        if (cpt->apertures.size() > 0 && this->N_L == 0)
        {
            this->N_L = cpt->apertures[0]->rows.size();
        }
        this->cpt_weights.push_back(cpt->apertures.size());
    }

    for (auto &cpt : this->electron_cpts)
    {
        cpt->remove_inactive_apertures();
        this->all_cpts.push_back(cpt);

        for (auto &ap : cpt->apertures)
        {
            this->active_weights.push_back(ap);
        }

        if (cpt->apertures.size() > 0 && this->N_L == 0)
        {
            this->N_L = cpt->apertures[0]->rows.size();
        }
        this->cpt_weights.push_back(cpt->apertures.size());
    }

    this->cpt_distribution = std::discrete_distribution<int>(this->cpt_weights.begin(), this->cpt_weights.end());
    this->max_iterations = 1000000;
    this->iteration_number = 0;
    this->n_W = 0;
    this->n_S = 0;
    this->N_A = this->active_weights.size();
    this->latest_cost = std::numeric_limits<double>::max();
    this->original_cost = std::numeric_limits<double>::max();

    // SA cooling constants
    this->P_shape = 0.9;
    this->T_P = 3;
    this->P_0 = 0.0;
    //this->P_0 = 0.0;
    this->T_S = 2;
    this->sigma_S0 = 3;
    this->T_W = 3;
    this->sigma_W0 = 0.1;
}

void SAOptimisation::from_json(json &json_input)
{
    this->total_voxels = json_input["total_voxels"];

    int num_cpts = json_input["electron_cpts"].size();
    std::cout << "Number of Electron control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["electron_cpts"])
    {
        ControlPoint *cpt = new ControlPoint(cpt_json);
        this->electron_cpts.push_back(cpt);
    }

    num_cpts = json_input["photon_cpts"].size();
    std::cout << "Number of Photon control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["photon_cpts"])
    {
        ControlPoint *cpt = new ControlPoint(cpt_json);
        this->photon_cpts.push_back(cpt);
    }

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

    this->start_time = std::chrono::steady_clock::now();
}

void SAOptimisation::launch()
{
    this->refresh_doses();
    this->latest_cost = this->calculate_cost();
    this->original_cost = this->latest_cost;
    this->lowest_cost = this->latest_cost;

    while (this->iteration_number < this->max_iterations)
    {
        if (this->iteration_number % 100 == 0 && this->iteration_number % 1000 != 0)
        {
            this->output_cost(false);
            this->update_data();
            this->export_data();
        }
        else if (this->iteration_number % 1000 == 0)
        {
            this->output_cost(true);
            this->export_dose(true);
            this->export_weights(true);
        }

        this->do_iteration();
        this->iteration_number += 1;
        if (this->latest_cost < this->lowest_cost)
            this->lowest_cost = this->latest_cost;
    }

    this->export_dose(false);
    this->export_weights(false);
}

void SAOptimisation::do_iteration()
{
    // Sample from control points, weighted by number of apertures in cpt.
    int cpt_index = this->cpt_distribution(this->rnd_gen);
    ControlPoint *cpt = this->all_cpts[cpt_index];

    // Decide between shape change and weight change.
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand1 = distribution(this->rnd_gen);
    if (rand1 < this->P_shape)
    {
        if (this->constrained && cpt_index < this->photon_cpts.size())
        {
            ControlPoint *prev_cpt = this->photon_cpts[0];
            ControlPoint *next_cpt = this->photon_cpts[0];

            for (int i = cpt_index - 1; i > 0; i--)
            {
                if (this->photon_cpts[i]->included == true)
                {
                    prev_cpt = this->photon_cpts[i];
                    break;
                }
            }

            for (int i = cpt_index + 1; i < this->photon_cpts.size(); i++)
            {
                if (this->photon_cpts[i]->included == true)
                {
                    next_cpt = this->photon_cpts[i];
                    break;
                }
            }

            this->do_constrained_shape_change(cpt, prev_cpt, next_cpt);
        }
        else
        {
            this->do_shape_change(cpt);
        }
    }
    else
    {
        this->do_weight_change(cpt);
    }
}

void SAOptimisation::do_shape_change(ControlPoint *cpt)
{
    std::uniform_int_distribution<int> uniform_apt(0, cpt->apertures.size() - 1);
    Aperture *aperture = cpt->apertures[uniform_apt(this->rnd_gen)];

    std::uniform_int_distribution<int> uniform_leaf(0, 2 * this->N_L - 1);
    int leaf_num = uniform_leaf(this->rnd_gen);

    //double exp_arg = -log(static_cast<float>(this->n_S) + 1) / this->T_S;
    double exp_arg = -log(static_cast<float>(this->n_S) / this->N_L + 1) / this->T_S;
    double sigma_S = 1 + (this->sigma_S0 - 1) * exp(exp_arg);
    std::normal_distribution<double> gaussian(0.0, sigma_S);
    // Coerce sampled leaf change to int to get an integer (beamlet) leaf change.
    int leaf_change = static_cast<int>(gaussian(this->rnd_gen));

    size_t row_num = leaf_num % this->N_L;
    size_t old_l_bound = aperture->rows[row_num].l_bound;
    size_t old_r_bound = aperture->rows[row_num].r_bound;
    size_t new_l_bound = old_l_bound;
    size_t new_r_bound = old_r_bound;
    // Left leaf
    if (leaf_num < this->N_L)
    {
        int new_pos = static_cast<int>(aperture->rows[row_num].l_bound) + leaf_change;
        if (new_pos < 0 || new_pos > aperture->rows[row_num].r_bound)
            return;

        new_l_bound = new_pos;
    }
    else
    {
        int new_pos = static_cast<int>(aperture->rows[row_num].r_bound) + leaf_change;
        if (new_pos < aperture->rows[row_num].l_bound || new_pos > cpt->num_beamlet_columns)
            return;
        new_r_bound = new_pos;
    }

    if (new_l_bound == old_l_bound && new_r_bound == old_r_bound)
        return;

    cpt->update_aperture_dose(aperture, row_num, new_l_bound, new_r_bound);
    this->refresh_doses();
    double new_cost = this->calculate_cost();
    if (new_cost < this->latest_cost)
    {
        this->latest_cost = new_cost;
        this->n_S += 1;
    }
    else
    {
        //double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) + 1) / this->T_P;
        double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) / this->N_A + 1) / this->T_P;
        double p_accept = 2.0 * this->P_0 / (1.0 + exp(exp_arg2));
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double rand1 = uniform(this->rnd_gen);

        // Accept a worse cost function with probability p_accept
        if (rand1 < p_accept)
        {
            this->latest_cost = new_cost;
            this->n_S += 1;
        }
        else
        {
            cpt->update_aperture_dose(aperture, row_num, old_l_bound, old_r_bound);
        }
    }
}

void SAOptimisation::do_constrained_shape_change(ControlPoint *cpt, ControlPoint *prev_cpt, ControlPoint *next_cpt)
{
    // Only one aperture per control point allowed for constrained opt.
    // cpt is guaranteed to have at least one aperture due to the sampling method.
    Aperture *aperture = cpt->apertures[0];

    std::uniform_int_distribution<int> uniform_leaf(0, 2 * this->N_L - 1);
    int leaf_num = uniform_leaf(this->rnd_gen);

    //double exp_arg = -log(static_cast<float>(this->n_S) + 1) / this->T_S;
    double exp_arg = -log(static_cast<float>(this->n_S) / this->N_L + 1) / this->T_S;
    double sigma_S = 1 + (this->sigma_S0 - 1) * exp(exp_arg);
    std::normal_distribution<double> gaussian(0.0, sigma_S);
    // Coerce sampled leaf change to int to get an integer (beamlet) leaf change.
    int leaf_change = static_cast<int>(gaussian(this->rnd_gen));

    size_t row_num = leaf_num % this->N_L;
    size_t old_l_bound = aperture->rows[row_num].l_bound;
    size_t old_r_bound = aperture->rows[row_num].r_bound;
    size_t new_l_bound = old_l_bound;
    size_t new_r_bound = old_r_bound;

    if (leaf_num < this->N_L)
    {
        // Left leaf
        int new_pos = static_cast<int>(aperture->rows[row_num].l_bound) + leaf_change;
        if (new_pos < 0 || new_pos > aperture->rows[row_num].r_bound)
            return;

        new_l_bound = new_pos;
    }
    else
    {
        int new_pos = static_cast<int>(aperture->rows[row_num].r_bound) + leaf_change;
        if (new_pos < aperture->rows[row_num].l_bound || new_pos > cpt->num_beamlet_columns)
            return;
        new_r_bound = new_pos;
    }

    if (new_l_bound == old_l_bound && new_r_bound == old_r_bound)
        return;

    if (!cpt->is_row_deliverable(prev_cpt, next_cpt, row_num, new_l_bound, new_r_bound))
        return;

    cpt->update_aperture_dose(aperture, row_num, new_l_bound, new_r_bound);
    this->refresh_doses();
    double new_cost = this->calculate_cost();
    if (new_cost < this->latest_cost)
    {
        this->latest_cost = new_cost;
        this->n_S += 1;
    }
    else
    {
        //double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) + 1) / this->T_P;
        double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) / this->N_A + 1) / this->T_P;
        double p_accept = 2.0 * this->P_0 / (1.0 + exp(exp_arg2));
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double rand1 = uniform(this->rnd_gen);

        // Accept a worse cost function with probability p_accept
        if (rand1 < p_accept)
        {
            this->latest_cost = new_cost;
            this->n_S += 1;
        }
        else
        {
            cpt->update_aperture_dose(aperture, row_num, old_l_bound, old_r_bound);
        }
    }
}

void SAOptimisation::do_weight_change(ControlPoint *cpt)
{
    std::uniform_int_distribution<int> uniform_apt(0, cpt->apertures.size() - 1);
    Aperture *aperture = cpt->apertures[uniform_apt(this->rnd_gen)];

    //double exp_arg = -log(static_cast<float>(this->n_W) + 1) / this->T_W;
    double exp_arg = -log(static_cast<float>(this->n_W) / this->N_A + 1) / this->T_W;
    double sigma_W = 0.01 + (this->sigma_W0 - 0.01) * exp(exp_arg);

    std::normal_distribution<double> gaussian(0.0, sigma_W);
    double weight_change = gaussian(this->rnd_gen);

    double old_weight = aperture->get_weight();
    double new_weight = std::max(0.0, old_weight * (1.0 + weight_change));
    aperture->set_weight(new_weight);

    this->refresh_doses();
    double new_cost = this->calculate_cost();
    if (new_cost < this->latest_cost)
    {
        this->latest_cost = new_cost;
        this->n_W += 1;
    }
    else
    {
        //double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) + 1.0) / this->T_P;
        double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) / this->N_A + 1.0) / this->T_P;
        double p_accept = 2.0 * this->P_0 / (1.0 + exp(exp_arg2));
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        double rand1 = uniform(this->rnd_gen);

        // Accept a worse cost function with probability p_accept
        if (rand1 < p_accept)
        {
            this->latest_cost = new_cost;
            this->n_W += 1;
        }
        else
        {
            aperture->set_weight(old_weight);
        }
    }
}

void SAOptimisation::calculate_structure_dose(Structure *structure)
{
    for (size_t i = 0; i < structure->masked_dose.size(); i++)
    {
        structure->masked_dose[i].second = 0.0;
        size_t vox_num = structure->masked_dose[i].first;
        for (auto &ap : this->active_weights)
        {
            structure->masked_dose[i].second += ap->dose[vox_num] * ap->weight;
        }
    }
    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
}

void SAOptimisation::refresh_doses()
{
    for (auto &structure : this->structures)
    {
        this->calculate_structure_dose(structure);
    }
}

double SAOptimisation::calculate_cost()
{
    double cost = 0.0;
    for (auto &structure : this->structures)
    {
        cost += structure->calculate_cost();
    }

    return cost;
}

double SAOptimisation::avg_aperture_size()
{
    double total_size = 0.0;
    double total_apertures = 0;

    for (auto &cpt : this->all_cpts)
    {
        for (auto &aperture : cpt->apertures)
        {
            total_size += aperture->size;
            total_apertures += 1;
        }
    }

    return total_size / total_apertures;
}

void SAOptimisation::update_data()
{
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->active_weights.size();
    iter_data["cost_function_value"] = this->latest_cost;
    iter_data["avg_aperture_size"] = this->avg_aperture_size();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - this->start_time).count();

    iter_data["running_time"] = elapsed;

    json e_cpt_stats;
    for (size_t i = 0; i < this->electron_cpts.size(); i++)
    {
        e_cpt_stats.push_back(this->electron_cpts[i]->write_statistics());
    }

    iter_data["electron_cpt_stats"] = e_cpt_stats;

    json p_cpt_stats;
    for (size_t i = 0; i < this->photon_cpts.size(); i++)
    {
        p_cpt_stats.push_back(this->photon_cpts[i]->write_statistics());
    }

    iter_data["photon_cpt_stats"] = p_cpt_stats;

    this->plan_data[this->iteration_number] = iter_data;
}

void SAOptimisation::output_plan_statistics()
{
    std::cout << this->electron_cpts.size() << " electron control points" << std::endl;
    for (auto &cpt : this->electron_cpts)
    {
        cpt->output_statistics();
    }

    std::cout << "---------------------------------------------------------" << std::endl;

    std::cout << this->photon_cpts.size() << " photon control points" << std::endl;
    for (auto &cpt : this->photon_cpts)
    {
        cpt->output_statistics();
    }
}

void SAOptimisation::output_cost(bool verbose)
{
    std::cout << "Original cost function value: " << this->original_cost << std::endl;
    std::cout << "Lowest cost function value: " << this->lowest_cost << std::endl;
    std::cout << "Latest cost function value: " << this->latest_cost << std::endl;
    std::cout << "Iteration " << this->iteration_number << "/" << this->max_iterations << std::endl;
    std::cout << "Accepted weight changes: " << this->n_W << std::endl;
    std::cout << "Accepted shape changes: " << this->n_S << std::endl;
    //double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) + 1.0) / this->T_P;
    double exp_arg2 = log(static_cast<float>(this->n_S + this->n_W) / N_A + 1.0) / this->T_P;
    double p_accept = 2.0 * this->P_0 / (1.0 + exp(exp_arg2));
    std::cout << "Probably of accepting worse solution: " << p_accept << std::endl;

    //double exp_arg = -log(static_cast<float>(this->n_S) + 1) / this->T_S;
    double exp_arg = -log(static_cast<float>(this->n_S) / this->N_L + 1) / this->T_S;
    double sigma_S = 1 + (this->sigma_S0 - 1) * exp(exp_arg);
    std::cout << "Sigma shape change: " << sigma_S << std::endl;

    //double exp_arg3 = -log(static_cast<float>(this->n_W) + 1) / this->T_W;
    double exp_arg3 = -log(static_cast<float>(this->n_W) / this->N_A + 1) / this->T_W;
    double sigma_W = 0.01 + (this->sigma_W0 - 0.01) * exp(exp_arg3);
    std::cout << "Sigma weight change: " << sigma_W << std::endl;

    if (verbose)
    {
        for (size_t i = 0; i < this->active_weights.size(); i++)
        {
            std::cout << "w[" << i << "] = " << this->active_weights[i]->weight << std::endl;
        }

        for (size_t i = 0; i < this->structures.size(); i++)
        {
            this->structures[i]->output_cost();
        }
    }
}

std::vector<double> SAOptimisation::calculate_total_dose()
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &ap : this->active_weights)
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += ap->dose[i] * ap->weight;
        }
    }

    return total_dose;
}

void SAOptimisation::export_data()
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

void SAOptimisation::export_dose(bool interm)
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

    if (this->electron_cpts.size())
    {
        num_voxels = this->electron_cpts[0]->num_voxels;
        topleft = this->electron_cpts[0]->topleft;
        voxel_size = this->electron_cpts[0]->voxel_size;
    }
    else
    {
        num_voxels = this->photon_cpts[1]->num_voxels;
        topleft = this->photon_cpts[1]->topleft;
        voxel_size = this->photon_cpts[1]->voxel_size;
    }

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

void SAOptimisation::export_weights(bool interm)
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
    json electron_cpts_json;
    json photon_cpts_json;

    for (auto &roi : this->structures)
    {
        structs_json.push_back(roi->to_json());
    }

    for (auto &cpt : this->electron_cpts)
    {
        electron_cpts_json.push_back(cpt->to_json());
    }

    for (auto &cpt : this->photon_cpts)
    {
        if (cpt->included)
        {
            photon_cpts_json.push_back(cpt->to_json());
        }
    }

    plan_json["structures"] = structs_json;
    plan_json["electron_cpts"] = electron_cpts_json;
    plan_json["photon_cpts"] = photon_cpts_json;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();
}
