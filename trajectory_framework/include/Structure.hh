#ifndef STRUCTURE_H
#define STRUCTURE_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"
#include "HardConstraint.hh"
#include "DVConstraint.hh"
#include "MeanConstraint.hh"
#include "PTVConstraint.hh"
#include "SparseDose.hh"

class Structure
{
  public:
    Structure();
    Structure(nlohmann::json &json_struct);
    std::string name;
    int total_voxels;
    int min_resolution;
    bool target;

    std::vector<std::pair<int, double>> masked_dose;

    std::vector<double> cost;
    std::vector<double> cdvh;
    std::vector<int> ddvh;

    std::vector<PTVConstraint> ptv_constraints;
    std::vector<HardConstraint> hard_constraints;
    std::vector<DVConstraint> dv_constraints;
    std::vector<MeanConstraint> mean_constraints;

    void read_struct(std::ifstream &input_file);
    void output_constraints();
    void write_constraints(std::ofstream &outfile);

    double calculate_cost();
    double calculate_gradient(std::vector<double> &cpt_dose);
    double calculate_gradient(SparseDose &dose);
    double calculate_hessian(std::vector<double> &cpt_one, std::vector<double> &cpt_two);
    double calculate_hessian(SparseDose &cpt_one, SparseDose &cpt_two);
    void calculate_lag_mults(std::vector<double> &lag_mults);
    void scale_weights(double scale);
    std::vector<double> calculate_dvh(bool relative = false);

    nlohmann::json get_dvh();
    nlohmann::json update_data();
    nlohmann::json to_json();
    nlohmann::json get_json_dvh();

    double dose_for_volume(double, bool percent = true);
    double volume_for_dose(double, bool percent = true);
    double get_min_dose();

    void set_voxel_volume(double vol);
    double voxel_volume;

    bool sorted;

    void output_cost();

    const int LOWER_LIMIT = 0;
    const int UPPER_LIMIT = 1;
};

#endif
