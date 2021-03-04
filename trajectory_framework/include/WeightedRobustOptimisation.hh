// Copyright 2018 Marc-Andre Renaud
#ifndef WeightedRobustOptimisation_h
#define WeightedRobustOptimisation_h
#include <vector>
#include <chrono>
#include "json.hh"

#include "enums.hh"
#include "PricingGenerator.hh"
#include "WeightClass.hh"
#include "Structure.hh"
#include "RobustControlPoint.hh"
#include "Coordinates.hh"

class WeightedRobustOptimisation {
    public:
        WeightedRobustOptimisation();
        WeightedRobustOptimisation(nlohmann::json &json_input);
        ~WeightedRobustOptimisation(){}

        void get_starting_point(double *weights);
        double calculate_cost_function(const double *weights, bool new_weights);
        void calculate_gradient(double *grad_values, const double *weights, bool new_weights);
        void calculate_hessian(double *hess_values, const double *weights, bool new_weights);

        int get_num_apertures();
        double get_highest_cost();

        bool virtual add_new_weight();

        void update_weights(const double *);
        void update_doses(const double *);
        void update_doses();

        void export_final_dose(bool interm=false);
        void export_final_weights(bool interm=false);
        void export_data();
        virtual void export_progress();

        size_t num_scenarios;
        size_t num_weights;
        std::string input_filename;

    protected:
        void from_json(nlohmann::json &json_input);

        void calculate_structure_dose(Structure *, std::vector<WeightClass *>);
        void update_weights(const double *, std::vector<WeightClass *>);
        void add_aperture(RobustControlPoint *cpt);
        void write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose);
        Coordinates get_coordinates();

        std::vector<double> calculate_total_dose();
        std::vector<double> calculate_scenario_dose(size_t);
        std::vector<double> calculate_modality_dose(std::vector<Aperture *> apertures);
        void export_scenario_doses(bool interm);
        void export_modality_doses();

        void update_data();
        void output_plan_statistics();
        double avg_aperture_size();
        void output_cost();

        void calculate_lag_mults();
        void prepare_pricing();
        bool prune_cpts();

        bool best_each_particle_scheme();
        bool alternating_particle_scheme();
        bool best_directional_price_scheme();
        bool best_per_modality_scheme();
        bool add_all_generated_apertures_scheme();
        bool random_scheme();

        std::vector<std::vector<double> > lag_mults;
        std::vector<double> scenario_weights;
        std::vector<RobustControlPoint *> photon_cpts;
        std::vector<RobustControlPoint *> electron_cpts;
        std::vector<std::vector<WeightClass *> > aperture_sets;
        std::vector<std::vector<Structure *> > structure_sets;

        int num_control_points;
        int max_apertures;
        size_t total_voxels;

        double latest_cost;
        int pruned_apertures;
        int iteration_number;
        bool last_ptype;
        bool enough_electrons;
        bool output_modalities;

        std::chrono::steady_clock::time_point start_time;

        std::map<int, double> energy_costs;
        std::map<int, int> energy_indices;

        nlohmann::json plan_data;
        std::vector<nlohmann::json> cost_data;
        std::vector<nlohmann::json> dvh_data;


        double max_photon_price;
        double max_electron_price;
        int max_electron_index;
        int max_photon_index;

        const bool ELECTRON = 0;
        const bool PHOTON = 1;

        mixing_scheme_choices mixing_scheme;
};

#endif