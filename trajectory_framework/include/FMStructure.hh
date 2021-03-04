#ifndef FMStructure_H
#define FMStructure_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"
#include "FMHardConstraint.hh"
#include "FMDVConstraint.hh"
#include "SparseDose.hh"

class FMStructure {
public:
  FMStructure();
  FMStructure(nlohmann::json &json_struct);
  std::string name;
  int roi_number;
  int total_voxels;
  std::vector<int> mask;
  std::vector<double> masked_dose;
  std::vector<double> cost;
  std::vector<double> cdvh;
  std::vector<int> ddvh;
  std::vector<std::vector<double> > beamlets;
  std::vector<FMHardConstraint> hard_constraints;
  std::vector<FMDVConstraint> dv_constraints;

  void read_struct(std::ifstream &input_file);
  void output_constraints();
  void write_constraints(std::ofstream &outfile);
  double calculate_cost();
  double calculate_gradient(std::vector<double> &dose);
  double calculate_hessian(std::vector<double> &cpt_one, std::vector<double> &cpt_two);
  void calculate_dvh();

  nlohmann::json to_json();

  void set_voxel_volume(double vol);
  double voxel_volume;

  void output_cost();

  const int LOWER_LIMIT = 0;
  const int UPPER_LIMIT = 1;
};

#endif
