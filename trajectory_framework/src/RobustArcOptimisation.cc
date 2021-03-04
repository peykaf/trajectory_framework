// Copyright 2018 Marc-Andre Renaud
#include <vector>
#include "json.hh"
#include "RobustControlPoint.hh"
#include "RobustArcOptimisation.hh"

using json = nlohmann::json;

RobustArcOptimisation::RobustArcOptimisation(json &json_input) : RobustOptimisation()
{
    std::cout << "Robust arc optimisation!" << std::endl;

    // Add a dummy control point for arcs. It's used in enforcing machine limitations.
    RobustControlPoint *dummyCpt = new RobustControlPoint();
    dummyCpt->included = true;
    this->photon_cpts.push_back(dummyCpt);

    this->input_filename = json_input["name"];
    this->from_json(json_input);

    for (auto &cpt : this->photon_cpts) {
        cpt->is_arc = true;
    }
}

void RobustArcOptimisation::price_unused_photon_cpts()
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
        RobustControlPoint *cur_cpt = this->photon_cpts[cpt_index];
        RobustControlPoint *prev_cpt = this->photon_cpts[0];
        RobustControlPoint *next_cpt = this->photon_cpts[0];
        if (!cur_cpt->included)
        {
            for (size_t i = cpt_index - 1; i > 0; i--)
            {
                RobustControlPoint *candidate_cpt = this->photon_cpts[i];
                if (candidate_cpt->included == true && candidate_cpt->arc_index == cur_cpt->arc_index)
                {
                    prev_cpt = candidate_cpt;
                    break;
                }
            }

            for (size_t i = cpt_index + 1; i < this->photon_cpts.size(); i++)
            {
                RobustControlPoint *candidate_cpt = this->photon_cpts[i];
                if (candidate_cpt->included == true && candidate_cpt->arc_index == cur_cpt->arc_index)
                {
                    next_cpt = candidate_cpt;
                    break;
                }
            }

            cur_cpt->price_constrained_aperture(this->lag_mults, this->kkt_mults, this->structure_sets, prev_cpt, next_cpt);
        }
    }

    std::cout << std::endl;
}

size_t RobustArcOptimisation::count_electron_apertures()
{
    size_t num_electron_apertures = 0;
    for (auto &cpt : this->electron_cpts)
    {
        num_electron_apertures += cpt->apertures.size();
    }

    return num_electron_apertures;
}

bool RobustArcOptimisation::add_new_aperture()
{
    bool pruned_cpts = this->prune_cpts();
    if (pruned_cpts)
        this->update_doses();

    if (this->max_apertures > 0)
    {
        this->enough_electrons = (this->count_electron_apertures() >= this->max_apertures);
    }

    this->prepare_pricing();
    this->calculate_lag_mults();

    this->price_unused_photon_cpts();

    if (!this->enough_electrons)
    {
        for (auto &cpt : this->electron_cpts)
        {
            cpt->price_aperture(this->lag_mults, this->kkt_mults, this->structure_sets);
        }
    }

    for (size_t i = 1; i < this->photon_cpts.size(); i++)
    {
        RobustControlPoint *cpt = this->photon_cpts[i];
        if (cpt->included == false && cpt->latest_price > this->max_photon_price)
        {
            this->max_photon_index = i;
            this->max_photon_price = cpt->latest_price;
        }
    }

    for (size_t i = 0; i < this->electron_cpts.size(); i++)
    {
        RobustControlPoint *cpt = this->electron_cpts[i];
        if (cpt->latest_price > this->max_electron_price)
        {
            this->max_electron_index = i;
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

void RobustArcOptimisation::export_progress()
{
    std::string filename2 = "combined_doses/" + remove_extension(this->input_filename) + ".progress";
    std::ofstream myfile2(filename2);
    myfile2 << this->aperture_sets[0].size() << "/" << (this->max_apertures + this->photon_cpts.size()) << std::endl;
    myfile2.close();
}
