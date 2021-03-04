// Copyright 2015 Marc-André Renaud
#ifndef SA_OPTIMISATION_H
#define SA_OPTIMISATION_H 1
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include "json.hh"
#include "TrajectoryGenerator.hh"
#include "ControlPoint.hh"
#include "Aperture.hh"
#include "FMStructure.hh"
#include "SparseDose.hh"
#include "MixedOptimisation.hh"
#include "MixedArcOptimisation.hh"

using json = nlohmann::json;

class SAOptimisation : public TrajectoryGenerator
{
  public:
    SAOptimisation(json &json_input);
    SAOptimisation(MixedOptimisation *mixed_opt);
    SAOptimisation(MixedArcOptimisation *mixed_opt);
    ~SAOptimisation(){};
    void read_input_file(std::string);
    double calculate_cost();
    void launch();

    size_t total_voxels;

    int num_voxels[3];
    double voxel_size[3];
    double topleft[3];

  private:
    void initialize();
    void from_json(json &json_input);
    void load_beamlets(std::vector<std::string> beamlet_filenames);
    void load_structure_beamlets(std::vector<std::string> beamlet_filenames);

    // Reads 3ddose header for dose grid information
    void read_header(std::string);
    void read_3ddose_header(std::string);
    void read_binary_header(std::string);

    // Various methods for calculating dose based on aperture weights
    void refresh_doses();
    void calculate_aperture_dose(Aperture &aperture);
    void calculate_structure_dose(Structure *structure);

    std::vector<double> calculate_total_dose();
    void do_iteration();
    void do_shape_change(ControlPoint *);
    void do_constrained_shape_change(ControlPoint *, ControlPoint *, ControlPoint *);
    void do_weight_change(ControlPoint *);

    double avg_aperture_size();
    void update_data();
    void output_plan_statistics();
    void output_cost(bool verbose = false);
    void export_data();
    void export_dose(bool interm);
    void export_weights(bool interm);

    // Containers for apertures and control points in various formats
    // for easy iteration.
    std::vector<ControlPoint *> photon_cpts;
    std::vector<ControlPoint *> electron_cpts;
    std::vector<ControlPoint *> all_cpts;
    std::vector<WeightClass *> active_weights;

    std::vector<Structure *> structures;

    std::default_random_engine rnd_gen;
    std::chrono::steady_clock::time_point start_time;

    std::vector<double> cpt_weights;
    std::discrete_distribution<int> cpt_distribution;

    int iteration_number; // Current iteration number
    int max_iterations;   // Max number of simulated annealing iterations

    nlohmann::json plan_data; // Plan metadata for analysis

    // CONSTANTS
    int N_A; // Total number of apertures
    int N_L; // Number of leaf pairs in MLC

    double P_shape; // Probability of shape change instead of weight change

    double sigma_W0; // sigma cooling schedule weight change
    double T_W;      // Cooling schedule weight change

    double sigma_S0; // sigma cooling schedule shape change
    double T_S;      // Cooling schedule shape change

    double T_P; // Cooling rate for accepting worse cost functions
    double P_0; // Initial probability of accepting worse cost functions

    // VARIABLE PARAMETERS
    int n_W; // Number of accepted weight changes
    int n_S; // Number of accepted shape changes

    bool constrained; // For arc optimisations, aperture shapes are constrained

    double latest_cost;
    double original_cost;
    double lowest_cost;
};

#endif
