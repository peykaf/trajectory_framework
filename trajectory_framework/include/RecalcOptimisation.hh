// Copyright 2015 Marc-André Renaud
#ifndef RECALC_OPTIMISATION_H
#define RECALC_OPTIMISATION_H 1
#include <string>
#include <vector>
#include "json.hh"

#include "TrajectoryGenerator.hh"
#include "Aperture.hh"
#include "Structure.hh"
#include "DenseDose.hh"

using json = nlohmann::json;

class RecalcOptimisation : public TrajectoryGenerator
{
  public:
    RecalcOptimisation();
    RecalcOptimisation(json &json_input);
    ~RecalcOptimisation();

    void from_json(json &json_input);
    double calculate_cost_function(const double *weights,
                                   bool new_weights = true);
    void calculate_gradient(double *grad_values,
                            const double *weights,
                            bool new_weights = true);
    void calculate_hessian(double *hess_values,
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
    void load_beamlets(std::vector<json> beamlet_sets);
    void load_sparse_beamlets(std::vector<std::string> beamlet_filenames);

    // Reads 3ddose header for dose grid information
    void read_header(std::string);
    void read_3ddose_header(std::string);
    void read_minidos_header(std::string);
    void read_bindos_header(std::string);

    // Various methods for calculating dose based on aperture weights
    void refresh_doses(const double *weights);
    void calculate_aperture_dose(Aperture &aperture);

    std::vector<double> calculate_total_dose(const double *weights);
    std::vector<double> calculate_total_dose(std::vector<double>);

    void export_modality_doses(const double *weights);
    void write_dose(const std::string filename, const std::vector<double> dose);
    std::vector<double> calculate_modality_dose(std::vector<int> apertures, const double *weights);

    void export_final_dose(const double *);
    void export_final_dose(std::vector<double>);

    void export_final_weights(std::vector<double>);
    void export_final_weights(const double *);

    void load_original_dose(std::string);

    void init_hessian();

    std::vector<double> optimised_weights;
    std::vector<Structure> structure_list;
    std::vector<DenseDose> beamlets;
    std::vector<double> hessian;
    DenseDose original_dose;
    json control_points;

    double latest_cost;
    bool spdose;

    bool output_modalities;

    const int LOWER_LIMIT = 0;
    const int UPPER_LIMIT = 1;
    const int TARGET = 0;
    const int OAR = 1;

    // double *hessian;
    // bool set_hessian;
};

#endif
