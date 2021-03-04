// Copyright 2015 Marc-André Renaud
#ifndef IPSA_OPTIMISATION_H
#define IPSA_OPTIMISATION_H 1
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include "json.hh"
#include "TrajectoryGenerator.hh"
#include "DwellPosition.hh"

using json = nlohmann::json;

class IPSAOptimisation : public TrajectoryGenerator
{
  public:
    IPSAOptimisation();
    IPSAOptimisation(json &json_input);
    ~IPSAOptimisation(){};
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

    // Reads 3ddose header for dose grid information
    void read_header(std::string);
    void read_3ddose_header(std::string);
    void read_binary_header(std::string);

    // Various methods for calculating dose based on aperture weights
    void refresh_doses();
    void refresh_doses(DwellPosition *, double old_weight);
    void calculate_structure_dose(Structure *structure);
    void update_structure_dose(Structure *structure, DwellPosition *dwell, double old_weight);

    std::vector<double> calculate_total_dose();
    double do_iteration();

    void update_data();
    void output_plan_statistics();
    void output_cost(bool verbose = false);
    void export_data();
    void export_dose(bool interm);
    void export_weights(bool interm);

    size_t calculate_active_dwells();

    // Containers for apertures and control points in various formats
    // for easy iteration.
    std::vector<DwellPosition *> dwell_positions;
    std::vector<WeightClass *> active_weights;

    std::vector<Structure *> structures;
    std::vector<std::vector<std::pair<int, double> > > structure_caches;

    std::default_random_engine rnd_gen;
    std::uniform_int_distribution<size_t> dwell_distribution;
    std::uniform_real_distribution<double> real_distribution;
    std::chrono::steady_clock::time_point start_time;

    int iteration_number; // Current iteration number
    int max_iterations;   // Max number of simulated annealing iterations

    nlohmann::json plan_data; // Plan metadata for analysis

    // IPSA CONSTANTS
    double alpha; // Speed constant
    double T_0; // Initial temperature

    double norm_constant;

    // VARIABLE PARAMETERS
    double latest_cost;
    double original_cost;
    double lowest_cost;
};

#endif
