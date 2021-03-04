// Copyright 2018 Marc-Andre Renaud
#include <vector>
#include "json.hh"
#include <chrono>
#include <random>
#include <map>
#include "RobustControlPoint.hh"
#include "WeightedRobustOptimisation.hh"

using json = nlohmann::json;

WeightedRobustOptimisation::WeightedRobustOptimisation() {
    this->mixing_scheme = mixing_scheme_choices::best_per_modality;
    this->max_apertures = 0;
    this->output_modalities = false;
    this->latest_cost = 0.0;
    this->iteration_number = 0;
    this->pruned_apertures = 0;
    this->enough_electrons = false;
    this->num_weights = 0;
}

WeightedRobustOptimisation::WeightedRobustOptimisation(json &json_input) : WeightedRobustOptimisation()
{
    std::cout << "Weighted robust optimisation!" << std::endl;
    this->input_filename = json_input["name"];

    this->from_json(json_input);
}

void WeightedRobustOptimisation::from_json(json &json_input)
{
    this->total_voxels = json_input["total_voxels"];
    this->num_scenarios = json_input["num_scenarios"];

    // Output dose distributions per modality
    this->output_modalities = json_input.value("output_modalities", false);

    this->aperture_sets.resize(this->num_scenarios);
    this->structure_sets.resize(this->num_scenarios);

    int num_cpts = json_input["electron_cpts"].size();
    std::cout << "Number of Electron control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["electron_cpts"])
    {
        RobustControlPoint *cpt = new RobustControlPoint(cpt_json, particle_type::electron);
        this->electron_cpts.push_back(cpt);
    }

    num_cpts = json_input["photon_cpts"].size();
    std::cout << "Number of Photon control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["photon_cpts"])
    {
        RobustControlPoint *cpt = new RobustControlPoint(cpt_json, particle_type::photon);
        this->photon_cpts.push_back(cpt);
    }


    std::vector<double> sc_weights(this->num_scenarios, 1.0);
    if (json_input.find("scenario_weights") != json_input.end())
    {
        sc_weights = json_input["scenario_weights"].get<std::vector<double>>();
    }

    double sc_sum = 0.0;
    for (size_t i = 0; i < sc_weights.size(); i++) {
        sc_sum += sc_weights[i];
    }
    for (size_t i = 0; i < sc_weights.size(); i++) {
        sc_weights[i] /= sc_sum;
    }

    this->scenario_weights = sc_weights;


    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    this->lag_mults.resize(this->num_scenarios);
    for (size_t i = 0; i < this->num_scenarios; i++)
    {
        this->lag_mults[i].resize(this->total_voxels);

        for (auto &struct_json : json_input["structures"])
        {
            bool include_robust = struct_json.value("robust", true);
            // if include_robust is false, only optimise the nominal scenario
            if (include_robust || i == 0) {
                Structure *new_struct = new Structure(struct_json);
                this->structure_sets[i].push_back(new_struct);
            }
        }
    }

    if (json_input.find("max_apertures") != json_input.end())
    {
        this->max_apertures = json_input["max_apertures"];
    }

    if (json_input.find("mixing_scheme") != json_input.end())
    {
        std::string mix_scheme = json_input["mixing_scheme"];

        std::map<std::string, mixing_scheme_choices> choices;
        choices["best_each_particle"] = mixing_scheme_choices::best_each_particle;
        choices["alternating_particle"] = mixing_scheme_choices::alternating_particle;
        choices["best_directional"] = mixing_scheme_choices::best_directional;
        choices["best_per_modality"] = mixing_scheme_choices::best_per_modality;
        choices["all_apertures"] = mixing_scheme_choices::all_apertures;
        choices["random"] = mixing_scheme_choices::random;

        this->mixing_scheme = choices[mix_scheme];
        std::cout << "Modality mixing scheme is: " << mix_scheme << std::endl;
    }

    std::cout << std::endl
              << "SCENARIO WEIGHTS: " << std::endl;
    for (size_t i = 0; i < this->scenario_weights.size(); i++)
    {
        double wt = this->scenario_weights[i];
        std::cout << wt;
        if (i != this->scenario_weights.size() - 1)
            std::cout << ", ";
    }
    std::cout << std::endl
              << std::endl;

    // Output constraints at the end of the input file reading.
    for (auto &structure : this->structure_sets[0])
    {
        structure->output_constraints();
    }
    std::cout << std::endl;

    this->start_time = std::chrono::steady_clock::now();
}

void WeightedRobustOptimisation::get_starting_point(double *weights)
{
    for (size_t i = 0; i < this->aperture_sets[0].size(); i++)
    {
        weights[i] = this->aperture_sets[0][i]->weight;
    }
}

void WeightedRobustOptimisation::update_weights(const double *weights)
{
    for (auto &apertures : this->aperture_sets)
    {
        for (size_t i = 0; i < apertures.size(); i++)
        {
            apertures[i]->set_weight(weights[i]);
        }
    }
}

void WeightedRobustOptimisation::update_weights(const double *weights, std::vector<WeightClass *> apertures)
{
    for (size_t i = 0; i < apertures.size(); i++)
    {
        apertures[i]->set_weight(weights[i]);
    }
}

void WeightedRobustOptimisation::calculate_structure_dose(Structure *structure, std::vector<WeightClass *> apertures)
{
    for (auto &vox : structure->masked_dose)
    {
        vox.second = 0.0;
        size_t vox_num = vox.first;
        for (auto &aperture : apertures)
        {
            vox.second += aperture->dose[vox_num] * aperture->weight;
        }
    }

    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
}

void WeightedRobustOptimisation::update_doses()
{
    for (size_t ap_index = 0; ap_index < this->aperture_sets.size(); ap_index++)
    {
        for (auto &structure : this->structure_sets[ap_index])
        {
            this->calculate_structure_dose(structure, this->aperture_sets[ap_index]);
        }
    }
}

void WeightedRobustOptimisation::update_doses(const double *weights)
{
    // Every aperture set has an associated structure set where doses are calculated.
    this->update_weights(weights);
    this->update_doses();
}

double WeightedRobustOptimisation::calculate_cost_function(const double *weights,
                                                           bool new_weights)
{
    double cost = 0.0;

    if (new_weights)
    {
        this->update_doses(weights);
    }

    for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
    {
        if (this->scenario_weights[sc_i] > 0.0)
        {
            double sc_cost = 0.0;
            for (auto &structure : this->structure_sets[sc_i])
            {
                sc_cost += structure->calculate_cost();
            }

            cost += sc_cost * this->scenario_weights[sc_i];
        }
    }

    this->latest_cost = cost;
    return cost;
}

void WeightedRobustOptimisation::calculate_gradient(double *values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        this->update_doses(weights);
    }

    for (size_t wt_i = 0; wt_i < this->num_weights; wt_i++)
    {
        double grad_value = 0.0;
        for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
        {
            if (this->scenario_weights[sc_i] > 0.0)
            {
                double scenario_grad = 0.0;
                for (auto &structure : this->structure_sets[sc_i])
                {
                    scenario_grad += structure->calculate_gradient(this->aperture_sets[sc_i][wt_i]->dose);
                }
                grad_value += scenario_grad * this->scenario_weights[sc_i];
            }
        }
        values[wt_i] = grad_value;
    }
}

void WeightedRobustOptimisation::calculate_hessian(double *hess_values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        this->update_doses(weights);
    }

    int idx = 0;
    for (size_t row = 0; row < this->num_weights; row++)
    {
        for (size_t col = 0; col <= row; col++)
        {
            double hess = 0.0;
            for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
            {
                if (this->scenario_weights[sc_i] > 0.0)
                {
                    double sc_hess = 0.0;
                    for (auto &structure : this->structure_sets[sc_i])
                    {
                        sc_hess += structure->calculate_hessian(this->aperture_sets[sc_i][row]->dose,
                                                                this->aperture_sets[sc_i][col]->dose);
                    }

                    hess += this->scenario_weights[sc_i] * sc_hess;
                }
            }

            hess_values[idx] = hess;
            idx++;
        }
    }
}

std::vector<double> WeightedRobustOptimisation::calculate_total_dose()
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : this->aperture_sets[0])
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

std::vector<double> WeightedRobustOptimisation::calculate_scenario_dose(size_t sc_i)
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : this->aperture_sets[sc_i])
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

std::vector<double> WeightedRobustOptimisation::calculate_modality_dose(std::vector<Aperture *> apertures) {
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : apertures)
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

int WeightedRobustOptimisation::get_num_apertures()
{
    return this->num_weights;
}

double WeightedRobustOptimisation::avg_aperture_size()
{
    double total_size = 0.0;
    double total_apertures = 0;

    for (auto &cpt : this->electron_cpts)
    {
        for (auto &aperture : cpt->apertures)
        {
            total_size += aperture->size;
            total_apertures += 1;
        }
    }

    for (auto &cpt : this->photon_cpts)
    {
        for (auto &aperture : cpt->apertures)
        {
            total_size += aperture->size;
            total_apertures += 1;
        }
    }

    return total_size / total_apertures;
}

void WeightedRobustOptimisation::update_data()
{
    this->iteration_number += 1;
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->aperture_sets[0].size();
    iter_data["cost_function_value"] = this->latest_cost;
    iter_data["avg_aperture_size"] = this->avg_aperture_size();
    iter_data["pruned_apertures"] = this->pruned_apertures;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - this->start_time).count();

    iter_data["running_time"] = elapsed;

    json e_cpt_stats;
    for (auto &cpt : this->electron_cpts)
    {
        e_cpt_stats.push_back(cpt->write_statistics());
    }
    iter_data["electron_cpt_stats"] = e_cpt_stats;

    json p_cpt_stats;
    for (auto &cpt : this->photon_cpts)
    {
        p_cpt_stats.push_back(cpt->write_statistics());
    }
    iter_data["photon_cpt_stats"] = p_cpt_stats;

    this->plan_data[this->iteration_number] = iter_data;

    // Cost data per iteration
    json iter_cost;
    iter_cost["cost_function_value"] = this->latest_cost;
    iter_cost["structures"] = std::vector<nlohmann::json>();

    for (auto &structure : this->structure_sets[0])
    {
        iter_cost["structures"].push_back(structure->update_data());
    }

    this->cost_data.push_back(iter_cost);

    // DVH data per iteration
    json iter_dvh;
    iter_dvh["structures"] = std::vector<nlohmann::json>();
    iter_dvh["dose_step"] = 0.1;
    for (auto &structure : this->structure_sets[0])
    {
        iter_dvh["structures"].push_back(structure->get_json_dvh());
    }
    this->dvh_data.push_back(iter_dvh);
}

void WeightedRobustOptimisation::output_plan_statistics()
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

/*
    PRICING PROBLEM RELATED
*/

bool WeightedRobustOptimisation::prune_cpts()
{
    double avg_weight = 0.0;
    int num_zero_weights = 0;
    int max_zero_weights = 0;
    bool pruned_cpts = 0;

    for (auto &aperture : this->aperture_sets[0])
    {
        avg_weight += aperture->weight;
    }

    avg_weight /= this->aperture_sets[0].size();

    for (auto &aperture : this->aperture_sets[0])
    {
        if (aperture->weight < avg_weight / 1000.0)
        {
            num_zero_weights += 1;
        }
    }

    if (num_zero_weights > max_zero_weights)
    {
        for (size_t sc_i = 0; sc_i < this->aperture_sets.size(); sc_i++)
        {
            for (size_t i = 0; i < this->aperture_sets[sc_i].size(); i++)
            {
                if (this->aperture_sets[sc_i][i]->weight < avg_weight / 1000.0)
                {
                    this->aperture_sets[sc_i][i]->deactivate();
                    this->pruned_apertures += 1;
                }
            }

            this->aperture_sets[sc_i].erase(std::remove_if(this->aperture_sets[sc_i].begin(), this->aperture_sets[sc_i].end(),
                                                           [](WeightClass *ap) { return !(ap->active); }),
                                            this->aperture_sets[sc_i].end());
        }

        for (auto &cpt : this->photon_cpts)
        {
            cpt->remove_inactive_apertures();
        }

        for (auto &cpt : this->electron_cpts)
        {
            cpt->remove_inactive_apertures();
        }

        pruned_cpts = true;
    }

    this->num_weights = this->aperture_sets[0].size();
    return pruned_cpts;
}

void WeightedRobustOptimisation::prepare_pricing()
{
    this->max_photon_index = -1;
    this->max_photon_price = 0.0;

    this->max_electron_index = -1;
    this->max_electron_price = 0.0;

    if (this->aperture_sets[0].size() > 0)
    {
        this->update_data();

        if (this->iteration_number % 5 == 0)
        {
            this->export_data();
        }
    }
    this->output_cost();

    // Export intermediate results every 10 aperture.
    if (this->aperture_sets[0].size() > 0 && this->iteration_number % 10 == 0)
    {
        this->export_final_dose(true);
        this->export_final_weights(true);
        this->output_plan_statistics();
    }
}

bool WeightedRobustOptimisation::add_new_weight()
{
    bool pruned_cpts = this->prune_cpts();

    if (this->max_apertures > 0 && this->aperture_sets[0].size() >= this->max_apertures)
    {
        return 0;
    }

    this->prepare_pricing();
    this->calculate_lag_mults();

    for (size_t cpt_i = 0; cpt_i < this->photon_cpts.size(); cpt_i++)
    {
        RobustControlPoint *cpt = this->photon_cpts[cpt_i];
        cpt->price_aperture(this->lag_mults, this->scenario_weights, this->structure_sets);
        if (cpt->latest_price > this->max_photon_price)
        {
            this->max_photon_index = cpt_i;
            this->max_photon_price = cpt->latest_price;
        }
    }

    std::cout << std::endl;

    for (size_t cpt_i = 0; cpt_i < this->electron_cpts.size(); cpt_i++)
    {
        RobustControlPoint *cpt = this->electron_cpts[cpt_i];
        cpt->price_aperture(this->lag_mults, this->scenario_weights, this->structure_sets);
        if (cpt->latest_price > this->max_electron_price)
        {
            this->max_electron_index = cpt_i;
            this->max_electron_price = cpt->latest_price;
        }
    }

    std::cout << "Max photon price: " << this->max_photon_price << std::endl;
    std::cout << "Max electron price: " << this->max_electron_price << std::endl;

    bool added_cpt = 0;
    switch (this->mixing_scheme) {
        case mixing_scheme_choices::best_each_particle: added_cpt = this->best_each_particle_scheme(); break;
        case mixing_scheme_choices::alternating_particle: added_cpt = this->alternating_particle_scheme(); break;
        case mixing_scheme_choices::best_directional: added_cpt = this->best_directional_price_scheme(); break;
        case mixing_scheme_choices::best_per_modality: added_cpt = this->best_per_modality_scheme(); break;
        case mixing_scheme_choices::all_apertures: added_cpt = this->add_all_generated_apertures_scheme(); break;
        case mixing_scheme_choices::random: added_cpt = this->random_scheme(); break;
    }

    return added_cpt;
}

void WeightedRobustOptimisation::calculate_lag_mults()
{
    for (size_t sc_i = 0; sc_i < this->aperture_sets.size(); sc_i++)
    {
        std::fill(this->lag_mults[sc_i].begin(), this->lag_mults[sc_i].end(), 0.0);

        if (this->scenario_weights[sc_i] > 0.0)
        {
            for (auto &structure : this->structure_sets[sc_i])
            {
                this->calculate_structure_dose(structure, this->aperture_sets[sc_i]);
            }

            for (auto &structure : this->structure_sets[sc_i])
            {
                structure->calculate_lag_mults(this->lag_mults[sc_i]);
            }
        }
    }
}

void WeightedRobustOptimisation::add_aperture(RobustControlPoint *cpt)
{
    std::cout << "Adding " << cpt->energy << " MV cpt" << std::endl;
    cpt->included = true;

    std::vector<Aperture *> apertures = cpt->save_robust_aperture();
    apertures[0]->print_aperture();

    for (size_t sc_i = 0; sc_i < this->aperture_sets.size(); sc_i++)
    {
        this->aperture_sets[sc_i].push_back(apertures[sc_i]);
    }

    this->num_weights = this->aperture_sets[0].size();
}

bool WeightedRobustOptimisation::random_scheme()
{
    std::cout << "Random scheme" << std::endl;
    std::vector<int> active_photons;
    std::vector<int> active_electrons;
    for (size_t i = 0; i < this->photon_cpts.size(); i++)
    {
        if (this->photon_cpts[i]->latest_price > 0.0)
        {
            active_photons.push_back(i);
        }
    }

    if (!this->enough_electrons)
    {
        for (size_t i = 0; i < this->electron_cpts.size(); i++)
        {
            if (this->electron_cpts[i]->latest_price > 0.0)
            {
                active_electrons.push_back(i);
            }
        }
    }

    std::random_device random_device;
    std::mt19937 engine{random_device()};

    int total_active = active_photons.size() + active_electrons.size();

    std::uniform_int_distribution<int> dist(0, total_active - 1);
    size_t random_aperture_index = dist(engine);

    bool added_cpt = 0;
    if (random_aperture_index < active_photons.size())
    {
        size_t cpt_index = active_photons[random_aperture_index];
        this->add_aperture(this->photon_cpts[cpt_index]);
        added_cpt = 1;
    }
    else
    {
        size_t cpt_index = active_electrons[random_aperture_index - active_photons.size()];
        this->add_aperture(this->electron_cpts[cpt_index]);
        added_cpt = 1;
    }

    return added_cpt;
}

bool WeightedRobustOptimisation::add_all_generated_apertures_scheme()
{
    // Adds the aperture from every modality and from every beam angle.
    // (# photon beam angles + # electron energies * # electron beam angles)

    bool added_cpt = 0;
    std::cout << "All generated apertures" << std::endl;

    for (auto &cpt : this->photon_cpts)
    {
        if (cpt->latest_price > 0.0)
        {
            this->add_aperture(cpt);
            added_cpt = 1;
        }
    }

    if (!this->enough_electrons)
        return added_cpt;

    for (auto &cpt : this->electron_cpts)
    {
        if (cpt->latest_price > 0.0)
        {
            this->add_aperture(cpt);
            added_cpt = 1;
        }
    }

    return added_cpt;
}

bool WeightedRobustOptimisation::best_each_particle_scheme()
{
    // Adds the best aperture from each particle type.
    // (2 apertures per iteration)

    bool added_cpt = 0;
    std::cout << "Best aperture per particle" << std::endl;
    if (this->max_photon_price > 0.0)
    {
        this->add_aperture(this->photon_cpts[this->max_photon_index]);
        added_cpt = 1;
    }

    if (this->max_electron_price > 0.0 && !this->enough_electrons)
    {
        this->add_aperture(this->electron_cpts[this->max_electron_index]);
        added_cpt = 1;
    }

    return added_cpt;
}

bool WeightedRobustOptimisation::alternating_particle_scheme()
{
    // Alternates between the best photon aperture and the best electron aperture.
    // (1 aperture per iteration)

    bool added_cpt = 0;
    std::cout << "Alternating particle scheme" << std::endl;

    if (this->last_ptype == ELECTRON && this->max_photon_price > 0)
    {
        this->add_aperture(this->photon_cpts[this->max_photon_index]);
        added_cpt = 1;
        this->last_ptype = PHOTON;
    }
    else if (this->last_ptype == PHOTON && this->max_electron_price > 0 && !this->enough_electrons)
    {
        this->add_aperture(this->electron_cpts[max_electron_index]);
        added_cpt = 1;
        this->last_ptype = ELECTRON;
    }

    // This happens only when one of the two modalities has 0 price.
    // We still want to add an aperture even if it's not that modality's turn.
    if (added_cpt == 0)
    {
        if (this->max_photon_price > 0)
        {
            this->add_aperture(this->photon_cpts[this->max_photon_index]);
            added_cpt = 1;
            this->last_ptype = PHOTON;
        }
        else if (this->max_electron_price > 0 && !this->enough_electrons)
        {
            this->add_aperture(this->electron_cpts[this->max_electron_index]);
            added_cpt = 1;
            this->last_ptype = ELECTRON;
        }
    }

    return added_cpt;
}

bool WeightedRobustOptimisation::best_directional_price_scheme()
{
    // Adds the best aperture out of all beam angles and modalities.
    // (1 aperture per iteration)

    bool added_cpt = 0;
    std::cout << "Best directional aperture" << std::endl;
    if (this->max_photon_price > 0.0 && this->max_photon_price >= this->max_electron_price)
    {
        this->add_aperture(this->photon_cpts[this->max_photon_index]);
        added_cpt = 1;
    }
    else if (this->max_electron_price > 0.0 && this->max_electron_price > this->max_photon_price && !this->enough_electrons)
    {
        this->add_aperture(this->electron_cpts[this->max_electron_index]);
        added_cpt = 1;
    }
    else if (this->max_photon_price > 0.0 && this->enough_electrons)
    {
        this->add_aperture(this->photon_cpts[this->max_photon_index]);
        added_cpt = 1;
    }

    return added_cpt;
}

bool WeightedRobustOptimisation::best_per_modality_scheme()
{
    // Adds the best aperture for each modality
    // (# of apertures is equal to number of modalities)
    // (eg. 1 photon aperture + # electron energies apertures)

    bool added_cpt = 0;
    std::vector<int> vints;
    std::cout << "Best per modality" << std::endl;

    this->energy_costs.clear();
    this->energy_indices.clear();

    for (size_t i = 0; i < this->electron_cpts.size(); i++)
    {
        ControlPoint *cpt = this->electron_cpts[i];
        if (energy_costs.count(cpt->energy) == 0)
        {
            energy_costs[cpt->energy] = cpt->latest_price;
            energy_indices[cpt->energy] = i;
        }

        if (cpt->latest_price > energy_costs[cpt->energy])
        {
            energy_costs[cpt->energy] = cpt->latest_price;
            energy_indices[cpt->energy] = i;
        }
    }

    if (this->max_photon_price > 0.0)
    {
        this->add_aperture(this->photon_cpts[this->max_photon_index]);
        added_cpt = 1;
    }

    if (this->enough_electrons)
        return added_cpt;

    for (auto imap : this->energy_costs)
    {
        vints.push_back(imap.first);
    }

    for (size_t i = 0; i < vints.size(); i++)
    {
        if (this->energy_costs[vints[i]] > 0.0)
        {
            size_t cpt_index = this->energy_indices[vints[i]];
            this->add_aperture(this->electron_cpts[cpt_index]);
            added_cpt = 1;
        }
    }

    return added_cpt;
}

/*
    OUTPUT ROUTINES
*/

Coordinates WeightedRobustOptimisation::get_coordinates() {
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
        // The first control point is a dummy control point in arc deliveries
        if (this->photon_cpts.size() > 1)
        {
            num_voxels = this->photon_cpts[1]->num_voxels;
            topleft = this->photon_cpts[1]->topleft;
            voxel_size = this->photon_cpts[1]->voxel_size;
        }
        else
        {
            num_voxels = this->photon_cpts[0]->num_voxels;
            topleft = this->photon_cpts[0]->topleft;
            voxel_size = this->photon_cpts[0]->voxel_size;
        }
    }


    return Coordinates(num_voxels, topleft, voxel_size);
}

void WeightedRobustOptimisation::output_cost()
{
    for (auto &structure : this->structure_sets[0])
    {
        structure->output_cost();
    }
}

void WeightedRobustOptimisation::export_modality_doses()
{
    std::vector<int> electron_energies;
    std::vector<int> photon_energies;
    std::map<int, std::vector<Aperture *> > electron_apertures;
    std::map<int, std::vector<Aperture *> > photon_apertures;

    std::vector<Aperture *> all_electrons;
    std::vector<Aperture *> all_photons;


    for (auto &cpt : this->electron_cpts)
    {
        if (electron_apertures.count(cpt->energy) == 0)
        {
            electron_apertures[cpt->energy] = std::vector<Aperture *>();
            electron_energies.push_back(cpt->energy);
        }

        for (auto &ap : cpt->apertures)
        {
            if (ap->active)
            {
                electron_apertures[cpt->energy].push_back(ap);
                all_electrons.push_back(ap);
            }
        }
    }

    for (auto &cpt : this->photon_cpts)
    {
        if (photon_apertures.count(cpt->energy) == 0)
        {
            photon_apertures[cpt->energy] = std::vector<Aperture *>();
            photon_energies.push_back(cpt->energy);
        }

        for (auto &ap : cpt->apertures)
        {
            if (ap->active)
            {
                photon_apertures[cpt->energy].push_back(ap);
                all_photons.push_back(ap);
            }
        }
    }

    Coordinates coords = this->get_coordinates();

    for (auto energy : electron_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MeV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(electron_apertures[energy]);
        this->write_dose(dose_filename, coords, modality_dose);
    }

    for (auto energy : photon_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(photon_apertures[energy]);
        this->write_dose(dose_filename, coords, modality_dose);
    }

    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_electrons.3ddose";
        std::cout << "Exporting electron dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(all_electrons);
        this->write_dose(dose_filename, coords, modality_dose);
    }

    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_photons.3ddose";
        std::cout << "Exporting photon dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(all_photons);
        this->write_dose(dose_filename, coords, modality_dose);
    }
}

void WeightedRobustOptimisation::export_scenario_doses(bool interm)
{
    std::string buffer;
    std::string dose_filename;

    Coordinates coords = this->get_coordinates();

    for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
    {
        if (sc_i == 0)
        {
            dose_filename = "combined_doses/" + remove_extension(input_filename);
        }
        else
        {
            dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(sc_i);
        }

        if (interm)
        {
            dose_filename += "_interm.3ddose";
        }
        else
        {
            dose_filename += ".3ddose";
        }

        std::cout << "Exporting scenario dose to " << dose_filename << "\n";

        std::vector<double> total_dose = this->calculate_scenario_dose(sc_i);
        this->write_dose(dose_filename, coords, total_dose);
    }
}

void WeightedRobustOptimisation::write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose) {
    std::ofstream dose_file(filename);

    dose_file << "  " << coords.num_voxels[0] << "  " << coords.num_voxels[1]
                << "  " << coords.num_voxels[2] << "\n";

    for (int i = 0; i < coords.num_voxels[0] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[0] + coords.topleft[0];
    }
    dose_file << "\n";

    for (int i = 0; i < coords.num_voxels[1] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[1] + coords.topleft[1];
    }
    dose_file << "\n";

    for (int i = 0; i < coords.num_voxels[2] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[2] + coords.topleft[2];
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

void WeightedRobustOptimisation::export_final_dose(bool interm)
{
    this->export_scenario_doses(interm);
    if (this->output_modalities) this->export_modality_doses();
}

void WeightedRobustOptimisation::export_final_weights(bool interm)
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

    for (auto &roi : this->structure_sets[0])
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

void WeightedRobustOptimisation::export_data()
{
    std::string filename = "combined_doses/" + remove_extension(this->input_filename) + ".data";
    std::ofstream myfile(filename);
    myfile << this->plan_data.dump();
    myfile.close();

    std::string filename3 = "combined_doses/" + remove_extension(this->input_filename) + ".cost";
    std::ofstream myfile3(filename3);
    myfile3 << nlohmann::json(this->cost_data).dump();
    myfile3.close();

    std::string filename4 = "combined_doses/" + remove_extension(this->input_filename) + ".dvh";
    std::ofstream myfile4(filename4);
    myfile4 << nlohmann::json(this->dvh_data).dump();
    myfile4.close();

    this->export_progress();
}

void WeightedRobustOptimisation::export_progress()
{
    std::string filename2 = "combined_doses/" + remove_extension(this->input_filename) + ".progress";
    std::ofstream myfile2(filename2);
    myfile2 << this->aperture_sets[0].size() << "/" << this->max_apertures << std::endl;
    myfile2.close();
}
