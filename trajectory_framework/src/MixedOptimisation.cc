// Copyright 2015 Marc-André Renaud
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>
#include <map>
#include "json.hh"

#include "enums.hh"
#include "string_utils.hh"
#include "MixedOptimisation.hh"
#include "ControlPoint.hh"
#include "Aperture.hh"
#include "Phantom.hh"
#include "Coordinates.hh"

using json = nlohmann::json;

MixedOptimisation::MixedOptimisation() : PricingGenerator() {
    this->mixing_scheme = mixing_scheme_choices::best_per_modality;
    this->num_pruned = 0;
    this->pruned_apertures = 0;
    this->max_zero_weights = 0;
    this->enough_electrons = false;
    this->output_modalities = false;
}

MixedOptimisation::MixedOptimisation(json &json_input) : MixedOptimisation()
{
    std::cout << "Mixed optimisation!" << std::endl;
    this->input_filename = json_input["name"];
    this->from_json(json_input);
}

void MixedOptimisation::from_json(json &json_input)
{
    this->total_voxels = json_input["total_voxels"];
    this->lag_mults.resize(this->total_voxels);

    // Output dose distributions per modality
    this->output_modalities = json_input.value("output_modalities", false);

    int num_cpts = json_input["electron_cpts"].size();
    std::cout << "Number of Electron control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["electron_cpts"])
    {
        ControlPoint *cpt = new ControlPoint(cpt_json, particle_type::electron);
        this->electron_cpts.push_back(cpt);
    }

    num_cpts = json_input["photon_cpts"].size();
    std::cout << "Number of Photon control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["photon_cpts"])
    {
        ControlPoint *cpt = new ControlPoint(cpt_json, particle_type::photon);
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

    if (json_input.find("max_apertures") != json_input.end())
    {
        this->max_apertures = json_input["max_apertures"];
        std::cout << "Max apertures: " << this->max_apertures << std::endl;
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

    this->start_time = std::chrono::steady_clock::now();
}

void MixedOptimisation::price_unused_photon_cpts()
{
    for (auto &cpt : this->photon_cpts)
    {
        cpt->price_aperture(this->lag_mults, this->structures);
    }
    std::cout << std::endl;
}

void MixedOptimisation::price_unused_electron_cpts()
{
    for (auto &cpt : this->electron_cpts)
    {
        cpt->price_aperture(this->lag_mults, this->structures);
    }
}

bool MixedOptimisation::prune_cpts()
{
    double avg_weight = 0.0;
    int num_zero_weights = 0;
    bool pruned = false;

    // If the cost function isn't improving, don't prune.
    double cost_improvement = std::abs((this->latest_cost - this->previous_cost) / this->latest_cost * 100.0);
    std::cout << "Cost improvement was " << cost_improvement << "%" << std::endl;
    if (cost_improvement < 0.1) {
        return false;
    }

    for (auto &aperture : this->active_weights)
    {
        avg_weight += aperture->weight;
    }

    avg_weight /= this->active_weights.size();

    for (auto &aperture : this->active_weights)
    {
        if (aperture->weight < avg_weight / 1000.0)
        {
            num_zero_weights += 1;
        }
    }

    if (num_zero_weights > max_zero_weights)
    {
        for (size_t i = 0; i < this->active_weights.size(); i++)
        {
            if (this->active_weights[i]->weight < avg_weight / 1000.0)
            {
                this->active_weights[i]->deactivate();
                this->pruned_apertures += 1;
                pruned = true;
            }
        }

        this->active_weights.erase(std::remove_if(this->active_weights.begin(), this->active_weights.end(),
                                                  [](WeightClass *ap) { return !(ap->active); }),
                                   this->active_weights.end());

        for (auto &cpt : this->photon_cpts)
        {
            cpt->remove_inactive_apertures();
        }

        for (auto &cpt : this->electron_cpts)
        {
            cpt->remove_inactive_apertures();
        }
    }

    this->num_weights = this->active_weights.size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;

    if (pruned) this->num_pruned += 1;

    return pruned;
}

double MixedOptimisation::avg_aperture_size()
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

void MixedOptimisation::update_data()
{
    this->iteration_number += 1;
    this->previous_cost = this->latest_cost;
    json iter_data;
    iter_data["iteration_number"] = this->iteration_number;
    iter_data["active_apertures"] = this->active_weights.size();
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

    for (auto &structure : this->structures)
    {
        iter_cost["structures"].push_back(structure->update_data());
    }

    this->cost_data.push_back(iter_cost);

    // DVH data per iteration
    json iter_dvh;
    iter_dvh["structures"] = std::vector<nlohmann::json>();
    iter_dvh["dose_step"] = 0.1;
    for (auto &structure : this->structures)
    {
        iter_dvh["structures"].push_back(structure->get_json_dvh());
    }
    this->dvh_data.push_back(iter_dvh);
}

void MixedOptimisation::output_plan_statistics()
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

bool MixedOptimisation::check_convergence()
{
    const int avg_aperture_iterations = 35;
    const double epsilon = 1e-6;

    double avg_apertures = 0.0;
    double avg_cost = 0.0;

    // Check if the current number of apertures is the same as the last 20
    // iterations' number of apertures. If so, it means that the optimiser
    // is probably always pruning the latest aperture added and has entered
    // an infinite loop.
    // TEMPORARY
    return false;

    if (this->last_apertures.size() >= avg_aperture_iterations)
    {
        for (std::list<int>::iterator p = this->last_apertures.begin(); p != this->last_apertures.end(); ++p)
        {
            avg_apertures += *p;
        }
        avg_apertures /= this->last_apertures.size();

        if (std::abs(avg_apertures - this->active_weights.size()) < epsilon)
        {
            std::cout << "Number of apertures converged" << std::endl;
            std::cout << "Avg apertures: " << avg_apertures << std::endl;
            std::cout << "latest num apertures: " << this->active_weights.size() << std::endl;
            return true;
        }

        this->last_apertures.pop_front();
    }
    // Check if the relative improvement of the cost function between this
    // iteration and the last avg_aperture_iterations is less than epsilon.
    // If so, each new aperture is probably contributing very little to
    // reducing the cost function, so we have converged.
    if (this->last_cost.size() >= avg_aperture_iterations)
    {
        for (std::list<double>::iterator p = this->last_cost.begin(); p != this->last_cost.end(); ++p)
        {
            avg_cost += *p;
        }
        avg_cost /= this->last_cost.size();

        if (std::abs(avg_cost - this->latest_cost) / this->latest_cost < epsilon)
        {
            std::cout << "Cost function converged" << std::endl;
            std::cout << "Avg_cost: " << avg_cost << std::endl;
            std::cout << "latest_cost: " << this->latest_cost << std::endl;
            return true;
        }

        this->last_cost.pop_front();
    }

    this->last_apertures.push_back(this->active_weights.size());
    this->last_cost.push_back(this->latest_cost);

    return false;
}

bool MixedOptimisation::add_new_weight()
{
    this->max_photon_index = -1;
    this->max_photon_price = 0.0;

    this->max_electron_index = -1;
    this->max_electron_price = 0.0;

    if (this->num_pruned < 5) {
        this->prune_cpts();
    } else {
        this->num_pruned = 0;
        this->max_zero_weights += 1;
    }

    this->refresh_doses();

    if (this->active_weights.size() > 0)
    {
        this->update_data();

        if (this->iteration_number % 5 == 0)
        {
            this->export_data();
        }
    }

    this->output_cost();

    if (this->check_convergence())
        return false;

    if (this->max_apertures > 0 && this->active_weights.size() >= this->max_apertures)
    {
        return false;
    }

    // Export intermediate results every 10 aperture.
    if (this->active_weights.size() > 0 && this->active_weights.size() % 10 == 0)
    {
        this->export_final_dose(true);
        this->export_final_weights(true);
        this->output_plan_statistics();
    }

    this->calculate_lag_mults();
    this->price_unused_photon_cpts();
    this->price_unused_electron_cpts();

    for (size_t i = 0; i < this->photon_cpts.size(); i++)
    {
        ControlPoint *cpt = this->photon_cpts[i];
        if (cpt->latest_price > max_photon_price)
        {
            max_photon_index = i;
            max_photon_price = cpt->latest_price;
        }
    }

    for (size_t i = 0; i < this->electron_cpts.size(); i++)
    {
        ControlPoint *cpt = this->electron_cpts[i];
        if (cpt->latest_price > max_electron_price)
        {
            max_electron_index = i;
            max_electron_price = cpt->latest_price;
        }
    }

    std::cout << "Max photon price: " << max_photon_price << std::endl;
    std::cout << "Max electron price: " << max_electron_price << std::endl;

    bool added_cpt = 0;

    switch (this->mixing_scheme) {
        case mixing_scheme_choices::best_each_particle: added_cpt = this->best_each_particle_scheme(); break;
        case mixing_scheme_choices::alternating_particle: added_cpt = this->alternating_particle_scheme(); break;
        case mixing_scheme_choices::best_directional: added_cpt = this->best_directional_price_scheme(); break;
        case mixing_scheme_choices::best_per_modality: added_cpt = this->best_per_modality_scheme(); break;
        case mixing_scheme_choices::all_apertures: added_cpt = this->add_all_generated_apertures_scheme(); break;
        case mixing_scheme_choices::random: added_cpt = this->random_scheme(); break;
    }

    // If max price is less than or equal to 0, then adding more control points
    // does not improve the dose distribution, we can terminate the algorithm.

    return added_cpt;
}

bool MixedOptimisation::random_scheme()
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

bool MixedOptimisation::add_all_generated_apertures_scheme()
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

bool MixedOptimisation::best_each_particle_scheme()
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

bool MixedOptimisation::alternating_particle_scheme()
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

bool MixedOptimisation::best_directional_price_scheme()
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

bool MixedOptimisation::best_per_modality_scheme()
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

void MixedOptimisation::add_aperture(ControlPoint *cpt)
{
    std::cout << "Adding " << cpt->energy << " MV cpt" << std::endl;
    cpt->included = true;

    Aperture *ap_ptr = cpt->save_current_aperture();
    ap_ptr->print_aperture();

    this->active_weights.push_back(ap_ptr);
    this->num_weights = this->active_weights.size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;
}

void MixedOptimisation::export_final_dose(bool interm)
{
    std::string buffer;
    std::string dose_filename;

    if (this->output_modalities) this->export_modality_doses();

    if (interm)
    {
        dose_filename = "combined_doses/" + remove_extension(input_filename) + "_interm.3ddose";
    }
    else
    {
        dose_filename = "combined_doses/" + remove_extension(input_filename) + ".3ddose";
    }

    std::cout << "Exporting final dose to " << dose_filename << "\n";

    Coordinates coords = this->get_coordinates();
    std::vector<double> total_dose = this->calculate_total_dose();
    this->write_dose(dose_filename, coords, total_dose);
}

Coordinates MixedOptimisation::get_coordinates() {
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

std::vector<double> MixedOptimisation::calculate_modality_dose(std::vector<Aperture *> apertures) {
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

void MixedOptimisation::export_modality_doses()
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

void MixedOptimisation::write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose) {
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

void MixedOptimisation::export_final_weights(bool interm)
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
        if (cpt->included && !cpt->dummy)
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
