// Copyright 2015 Marc-André Renaud
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "json.hh"
#include "DwellPosition.hh"
#include "SparseDose.hh"
#include "DenseDose.hh"
#include "string_utils.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t> > sparse_vector;

DwellPosition::DwellPosition() {
    this->included = false;
    this->weight = 0.0;
    this->position.resize(3);
}

DwellPosition::DwellPosition(json &dwell_json) : DwellPosition() {
    this->catheter_index = dwell_json["catheter_index"];
    this->angle_index = dwell_json["angle_index"];
    this->shield_angle = dwell_json["shield_angle"];
    this->position = dwell_json["position"].get<std::vector<double> >();
    this->dose_filename = dwell_json["dose_filename"];

    this->read_dose(this->dose_filename);
}

void DwellPosition::read_3ddose_header(std::string filename) {

    std::string buffer;
    std::string trimmed_buffer;
    std::vector<std::string> voxel_buffer;

    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer);  // number of voxels
    std::stringstream convertor(buffer);
    convertor >> this->num_voxels[0] >> this->num_voxels[1] >> this->num_voxels[2];
    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
                                       << this->num_voxels[1] << ","
                                       << this->num_voxels[2] << ")" << std::endl;

    double second_voxel;
    std::getline(beamlet_file, buffer);  // x voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[0];
    convertor >> second_voxel;
    std::cout << trimmed_buffer << std::endl;
    this->voxel_size[0] = second_voxel - this->topleft[0];

    std::getline(beamlet_file, buffer);  // y voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[1];
    convertor >> second_voxel;
    this->voxel_size[1] = second_voxel - this->topleft[1];

    std::getline(beamlet_file, buffer);  // z voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[2];
    convertor >> second_voxel;
    this->voxel_size[2] = second_voxel - this->topleft[2];

    std::cout << "Voxel size: (" << this->voxel_size[0] << ","
                                 << this->voxel_size[1] << ","
                                 << this->voxel_size[2] << ")" << std::endl;

    std::cout << "Topleft: (" << this->topleft[0] << ","
                                 << this->topleft[1] << ","
                                 << this->topleft[2] << ")" << std::endl;

    beamlet_file.close();

}

void DwellPosition::read_binary_header(std::string filename) {
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    float b_topleft[3], b_voxel_size[3];
    beamlet_file.read((char*)this->num_voxels, 3 * sizeof(int));
    beamlet_file.read((char*)b_voxel_size, 3 * sizeof(float));
    beamlet_file.read((char*)b_topleft, 3 * sizeof(float));
    this->voxel_size[0] = b_voxel_size[0];
    this->voxel_size[1] = b_voxel_size[1];
    this->voxel_size[2] = b_voxel_size[2];

    this->topleft[0] = b_topleft[0];
    this->topleft[1] = b_topleft[1];
    this->topleft[2] = b_topleft[2];

    this->total_voxels = this->num_voxels[0]
                         * this->num_voxels[1]
                         * this->num_voxels[2];
    beamlet_file.close();
}

void DwellPosition::read_dose(std::string filename) {
  std::string extension;
  std::vector<std::string> buffer;

  buffer = split(filename, '.');
  extension = buffer[buffer.size()-1];
  DenseDose d_dose;
  if (extension == "3ddose") {
      this->read_3ddose_header(filename);
      d_dose.from_3ddose(filename, true);
  } else if (extension == "minidos") {
      this->read_binary_header(filename);
      d_dose.from_minidos(filename, true);
  } else {
      throw "Dose file extension not recognised";
  }

  this->dose = d_dose.grid;
}

void DwellPosition::price_dwell(std::vector<double> &lag_mults) {
    this->latest_price = 0.0;

    for (size_t i = 0; i < lag_mults.size(); i++) {
        this->latest_price += this->dose[i] * lag_mults[i];
    }
}

json DwellPosition::to_json() {
    json dwell;
    dwell["catheter_index"] = this->catheter_index;
    dwell["angle_index"] = this->angle_index;
    dwell["shield_angle"] = this->shield_angle;
    dwell["weight"] = this->weight;
    dwell["position"] = this->position;

    return dwell;
}

json DwellPosition::write_statistics() {
    return this->to_json();
}

void DwellPosition::output_statistics() {
    std::cout << "Dwell position" << std::endl;
    std::cout << "Cathether index: " << this->catheter_index << std::endl;
    std::cout << "Angle index: " << this->angle_index << std::endl;
    std::cout << "Shield angle: " << this->shield_angle << std::endl;
    std::cout << "Position: (" << this->position[0] << ", " << this->position[1] << "," << this->position[2] << ")" << std::endl;
    std::cout << "Dwell time: " << this->weight << std::endl;
}