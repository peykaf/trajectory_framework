// Copyright 2017 Marc-André Renaud
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "json.hh"

#include "MixedArcOptimisation.hh"
#include "MixedOptimisation.hh"
#include "ControlPoint.hh"

using json = nlohmann::json;

MixedArcOptimisation::MixedArcOptimisation(json &json_input) : MixedOptimisation()
{
    std::cout << "Mixed arc optimisation!" << std::endl;

    // Add a dummy control point for arcs. It's used in enforcing machine limitations.
    ControlPoint *dummyCpt = new ControlPoint();
    dummyCpt->included = true;
    dummyCpt->dummy = true;
    this->photon_cpts.push_back(dummyCpt);

    this->input_filename = json_input["name"];
    this->from_json(json_input);

    for (auto &cpt : this->photon_cpts) {
        cpt->is_arc = true;
    }
}

void MixedArcOptimisation::price_unused_photon_cpts()
{
    // For every control point not included in the plan yet,
    // price an aperture and select the control point with the best
    // aperture (highest price). The pricing routine requires knowledge
    // of the previous and next control point to account for machine
    // limitations.

    // Could do this in one pass, but I'm doing the less efficient version of this loop
    // to make sure there are no bugs. Might refactor later if it actually affects speed.
    for (size_t cpt_index = 1; cpt_index < this->photon_cpts.size(); cpt_index++)
    {
        ControlPoint *cur_cpt = this->photon_cpts[cpt_index];
        ControlPoint *prev_cpt = this->photon_cpts[0];
        ControlPoint *next_cpt = this->photon_cpts[0];
        if (!cur_cpt->included)
        {
            for (size_t i = cpt_index - 1; i > 0; i--)
            {
                ControlPoint *candidate_cpt = this->photon_cpts[i];
                if (candidate_cpt->included == true && candidate_cpt->arc_index == cur_cpt->arc_index)
                {
                    prev_cpt = candidate_cpt;
                    break;
                }
            }

            for (size_t i = cpt_index + 1; i < this->photon_cpts.size(); i++)
            {
                ControlPoint *candidate_cpt = this->photon_cpts[i];
                if (candidate_cpt->included == true && candidate_cpt->arc_index == cur_cpt->arc_index)
                {
                    next_cpt = candidate_cpt;
                    break;
                }
            }

            cur_cpt->price_constrained_aperture(this->lag_mults, this->structures, prev_cpt, next_cpt);
        }
        else
        {
            cur_cpt->latest_price = 0.0;
        }
    }

    std::cout << std::endl;
}

size_t MixedArcOptimisation::count_electron_apertures()
{
    size_t num_electron_apertures = 0;
    for (auto &cpt : this->electron_cpts)
    {
        num_electron_apertures += cpt->apertures.size();
    }

    return num_electron_apertures;
}

bool MixedArcOptimisation::add_new_weight()
{
    this->max_photon_index = -1;
    this->max_photon_price = 0.0;

    this->max_electron_index = -1;
    this->max_electron_price = 0.0;
    this->energy_costs.clear();
    this->energy_indices.clear();

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

        if (this->iteration_number % 10 == 0)
        {
            this->export_final_dose(true);
            this->export_final_weights(true);
            this->output_plan_statistics();
        }
    }

    this->output_cost();

    if (this->check_convergence())
        return false;

    if (this->max_apertures > 0)
    {
        this->enough_electrons = (this->count_electron_apertures() >= this->max_apertures);
    }

    this->calculate_lag_mults();
    this->price_unused_photon_cpts();

    if (!this->enough_electrons)
        this->price_unused_electron_cpts();

    // Start loop at i = 1 to avoid the dummy control point
    for (size_t i = 1; i < this->photon_cpts.size(); i++)
    {
        ControlPoint *cpt = this->photon_cpts[i];
        if (cpt->included == false && cpt->latest_price > max_photon_price)
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

    std::cout << "Max photon price: " << max_photon_price << std::endl;
    std::cout << "Max electron price: " << max_electron_price << std::endl;

    bool added_cpt = false;
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

    if (!added_cpt)
        this->export_data();

    return added_cpt;
}

void MixedArcOptimisation::export_progress()
{
    std::string filename2 = "combined_doses/" + remove_extension(this->input_filename) + ".progress";
    std::ofstream myfile2(filename2);
    myfile2 << this->active_weights.size() << "/" << (this->max_apertures + this->photon_cpts.size()) << std::endl;
    myfile2.close();
}
