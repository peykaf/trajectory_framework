// Copyright 2015 Marc-André Renaud
#ifndef MixedOptimisation_H
#define MixedOptimisation_H 1

#include <string>
#include <vector>
#include <chrono>
#include "json.hh"
#include "enums.hh"

#include "PricingGenerator.hh"
#include "Phantom.hh"
#include "Structure.hh"
#include "ControlPoint.hh"
#include "Coordinates.hh"

class MixedOptimisation : public PricingGenerator {
    public:
        MixedOptimisation();
        MixedOptimisation(nlohmann::json &json_input);

        virtual void from_json(nlohmann::json &json_input);

        bool add_new_weight();
        virtual bool prune_cpts();

        virtual void export_final_dose(bool interm=false);
        void export_final_weights(bool interm=false);

        // GETTERS
        std::vector<ControlPoint *> get_photon_cpts() {return this->photon_cpts;}
        std::vector<ControlPoint *> get_electron_cpts() {return this->electron_cpts;}

        size_t get_total_voxels() {return this->total_voxels;}
    protected:
        std::vector<double> calculate_modality_dose(std::vector<Aperture *> apertures);
        void export_modality_doses();
        Coordinates get_coordinates();
        void write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose);

        virtual void price_unused_photon_cpts();
        void price_unused_electron_cpts();

        void output_plan_statistics();

        virtual void add_aperture(ControlPoint *);

        double avg_aperture_size();
        void update_data();

        bool check_convergence();

        bool best_each_particle_scheme();
        bool alternating_particle_scheme();
        bool best_directional_price_scheme();
        bool best_per_modality_scheme();
        bool add_all_generated_apertures_scheme();
        bool random_scheme();

        int num_control_points;
        int num_pruned;
        int max_zero_weights;  // Max zero weights apertures needed before pruning
        bool last_ptype;
        bool enough_electrons;
        bool output_modalities;

        int pruned_apertures;

        std::vector<ControlPoint *> photon_cpts;
        std::vector<ControlPoint *> electron_cpts;

        std::chrono::steady_clock::time_point start_time;

        std::map<int, double> energy_costs;
        std::map<int, int> energy_indices;

        double max_photon_price;
        double max_electron_price;
        int max_electron_index;
        int max_photon_index;

        const bool ELECTRON = 0;
        const bool PHOTON = 1;


        mixing_scheme_choices mixing_scheme;
};

#endif
