// Copyright 2015 Marc-André Renaud
#ifndef FLUENCE_MAP_OPTIMISATION_H
#define FLUENCE_MAP_OPTIMISATION_H 1
#include <string>
#include <vector>
#include "json.hh"
#include "TrajectoryGenerator.hh"
#include "Aperture.hh"
#include "FMStructure.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

class FluenceMapOptimisation : public TrajectoryGenerator {
  public:
    FluenceMapOptimisation(json &json_input);
    FluenceMapOptimisation();
    ~FluenceMapOptimisation();
    void read_input_file(std::string);
    double calculate_cost_function(const double *weights,
                                   bool new_weights = true);
    void calculate_gradient(double *grad_values,
                            const double *weights,
                            bool new_weights = true);
    void calculate_hessian(double* hess_values,
                           const double *weights,
                           bool new_weights = true);

    void calculate_jacobian() {}
    void calculate_constraints() {}

    size_t total_voxels;
    size_t beamlet_rows;
    size_t beamlet_columns;

    int num_voxels[3];
    double voxel_size[3];
    double topleft[3];

  private:
    void from_json(json &json_input);
    void load_beamlets(std::vector<std::string> beamlet_filenames);
    void load_structure_beamlets(std::vector<std::string> beamlet_filenames);

    // Reads 3ddose header for dose grid information
    void read_header(std::string);
    void read_3ddose_header(std::string);
    void read_binary_header(std::string);

    // Various methods for calculating dose based on aperture weights
    void refresh_doses(const double* weights);
    void calculate_aperture_dose(Aperture& aperture);

    std::vector<double> calculate_total_dose(const double* weights);

    void export_final_dose(const double *);
    void export_final_weights(const double *);

    std::vector<double> optimised_weights;
    std::vector<FMStructure> structures;
    std::vector<SparseDose> beamlets;

    const int LOWER_LIMIT = 0;
    const int UPPER_LIMIT = 1;
    const int TARGET = 0;
    const int OAR = 1;

    double latest_cost;

    std::vector<std::string> beamlet_filenames;

    // double *hessian;
    // bool set_hessian;
};

#endif
