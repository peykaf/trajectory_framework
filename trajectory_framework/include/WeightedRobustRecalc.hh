// Copyright 2018 Marc-Andre Renaud
#ifndef WeightedRobustRecalc_h
#define WeightedRobustRecalc_h
#include <vector>
#include <chrono>
#include "json.hh"

#include "PricingGenerator.hh"
#include "WeightClass.hh"
#include "Structure.hh"
#include "RobustRecalcControlPoint.hh"
#include "Coordinates.hh"
#include "RecalcAperture.hh"
#include "TrajectoryGenerator.hh"

class WeightedRobustRecalc : TrajectoryGenerator {
    public:
        WeightedRobustRecalc(){}
        WeightedRobustRecalc(nlohmann::json &json_input);
        ~WeightedRobustRecalc(){}

        void get_starting_point(double *weights);
        double calculate_cost_function(const double *weights, bool new_weights);
        void calculate_gradient(double *grad_values, const double *weights, bool new_weights);
        void calculate_hessian(double *hess_values, const double *weights, bool new_weights);

        int get_num_apertures();
        double get_highest_cost();

        void update_weights(const double *);
        void update_doses(const double *);
        void update_doses();

        void export_final_dose(bool interm=false);
        void export_final_weights(bool interm=false);
        void export_data();

        size_t num_scenarios;
        size_t num_weights;
        std::string input_filename;

    protected:
        void from_json(nlohmann::json &json_input);
        void init_apertures();

        void calculate_structure_dose(Structure *, std::vector<RecalcAperture *>);
        void update_weights(const double *, std::vector<RecalcAperture *>);
        void write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose);
        Coordinates get_coordinates();

        std::vector<double> calculate_total_dose();
        std::vector<double> calculate_scenario_dose(size_t);
        std::vector<double> calculate_modality_dose(std::vector<RecalcAperture *> apertures);
        void export_scenario_doses(bool interm);
        void export_modality_doses();

        void update_data();
        void output_plan_statistics();
        double avg_aperture_size();
        void output_cost();

        std::vector<double> scenario_weights;
        std::vector<RobustRecalcControlPoint *> photon_cpts;
        std::vector<RobustRecalcControlPoint *> electron_cpts;
        std::vector<std::vector<RecalcAperture *> > aperture_sets;
        std::vector<std::vector<Structure *> > structure_sets;

        int num_control_points;
        size_t total_voxels;

        double latest_cost;
        bool output_modalities;

        std::chrono::steady_clock::time_point start_time;

        nlohmann::json plan_data;
};

#endif