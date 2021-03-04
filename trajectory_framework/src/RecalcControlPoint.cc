// Copyright 2015 Marc-André Renaud
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include "json.hh"
#include "RecalcControlPoint.hh"
#include "SparseDose.hh"
#include "string_utils.hh"
#include "Structure.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t>> sparse_vector;

RecalcControlPoint::RecalcControlPoint()
{
    this->gantry_angle = -9999;
    this->couch_angle = 0.0;
    this->collimator_angle = 0.0;
    this->arc_index = 1;
    this->particle = this->photon;
    this->arclength_scaling = 1.0;
    this->FFF = false;
}

RecalcControlPoint::RecalcControlPoint(json &cpt_json) : RecalcControlPoint()
{
    this->energy = cpt_json["energy"];

    if (cpt_json.find("arc_index") != cpt_json.end())
    {
        this->arc_index = cpt_json["arc_index"];
    }

    if (cpt_json.find("FFF") != cpt_json.end())
    {
        this->FFF = cpt_json["FFF"];
    }

    this->gantry_angle = cpt_json["gantry_angle"];
    this->couch_angle = cpt_json["couch_angle"];
    this->collimator_angle = cpt_json["col_angle"];

    this->sad = cpt_json["sad"];
    this->iso = cpt_json["iso"].get<std::vector<double>>();

    this->coords = this->read_header(cpt_json["apertures"][0]["dose_filename"]);

    for (auto &ap : cpt_json["apertures"]) {
        this->apertures.push_back(new RecalcAperture(ap));
    }
}

RecalcControlPoint::RecalcControlPoint(json &cpt_json, particle_type ptype) : RecalcControlPoint(cpt_json)
{
    this->particle = ptype;
}

Coordinates RecalcControlPoint::read_3ddose_header(std::string filename)
{
    /*
      Populates the class attributes related to phantom dimensions.
    */

    Coordinates coords;
    std::string buffer;
    std::string trimmed_buffer;
    std::vector<std::string> voxel_buffer;

    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer); // number of voxels
    std::stringstream convertor(buffer);
    convertor >> coords.num_voxels[0] >> coords.num_voxels[1] >> coords.num_voxels[2];

    std::cout << "Number of voxels: (" << coords.num_voxels[0] << ","
              << coords.num_voxels[1] << ","
              << coords.num_voxels[2] << ")" << std::endl;

    double second_voxel;
    std::getline(beamlet_file, buffer); // x voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> coords.topleft[0];
    convertor >> second_voxel;
    std::cout << trimmed_buffer << std::endl;
    coords.voxel_size[0] = second_voxel - coords.topleft[0];

    std::getline(beamlet_file, buffer); // y voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> coords.topleft[1];
    convertor >> second_voxel;
    coords.voxel_size[1] = second_voxel - coords.topleft[1];

    std::getline(beamlet_file, buffer); // z voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> coords.topleft[2];
    convertor >> second_voxel;
    coords.voxel_size[2] = second_voxel - coords.topleft[2];

    std::cout << "Voxel size: (" << coords.voxel_size[0] << ","
              << coords.voxel_size[1] << ","
              << coords.voxel_size[2] << ")" << std::endl;

    std::cout << "Topleft: (" << coords.topleft[0] << ","
              << coords.topleft[1] << ","
              << coords.topleft[2] << ")" << std::endl;

    beamlet_file.close();

    return coords;
}

Coordinates RecalcControlPoint::read_minidos_header(std::string filename)
{
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    float b_topleft[3], b_voxel_size[3];
    int b_num_voxels[3];
    beamlet_file.read((char *)b_num_voxels, 3 * sizeof(int));
    beamlet_file.read((char *)b_voxel_size, 3 * sizeof(float));
    beamlet_file.read((char *)b_topleft, 3 * sizeof(float));

    coords.num_voxels[0] = b_num_voxels[0];
    coords.num_voxels[1] = b_num_voxels[1];
    coords.num_voxels[2] = b_num_voxels[2];

    coords.voxel_size[0] = b_voxel_size[0];
    coords.voxel_size[1] = b_voxel_size[1];
    coords.voxel_size[2] = b_voxel_size[2];

    coords.topleft[0] = b_topleft[0];
    coords.topleft[1] = b_topleft[1];
    coords.topleft[2] = b_topleft[2];

    beamlet_file.close();

    return coords;
}

Coordinates RecalcControlPoint::read_bindos_header(std::string filename)
{
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    int b_num_voxels[3];
    beamlet_file.read((char *)b_num_voxels, 3 * sizeof(int));

    float *x_voxels = new float[b_num_voxels[0] + 1];
    float *y_voxels = new float[b_num_voxels[1] + 1];
    float *z_voxels = new float[b_num_voxels[2] + 1];

    beamlet_file.read((char *)x_voxels, (b_num_voxels[0] + 1) * sizeof(float));
    beamlet_file.read((char *)y_voxels, (b_num_voxels[1] + 1) * sizeof(float));
    beamlet_file.read((char *)z_voxels, (b_num_voxels[2] + 1) * sizeof(float));

    beamlet_file.close();

    Coordinates coords;
    coords.num_voxels[0] = b_num_voxels[0];
    coords.num_voxels[1] = b_num_voxels[1];
    coords.num_voxels[2] = b_num_voxels[2];

    coords.voxel_size[0] = x_voxels[1] - x_voxels[0];
    coords.voxel_size[1] = y_voxels[1] - y_voxels[0];
    coords.voxel_size[2] = z_voxels[1] - z_voxels[0];

    coords.topleft[0] = x_voxels[0];
    coords.topleft[1] = y_voxels[0];
    coords.topleft[2] = z_voxels[0];

    delete[] x_voxels;
    delete[] y_voxels;
    delete[] z_voxels;

    return coords;
}

Coordinates RecalcControlPoint::read_header(std::string filename)
{
    std::string extension;

    Coordinates coords;
    extension = find_extension(filename);
    if (extension == "3ddose")
    {
        coords = this->read_3ddose_header(filename);
    }
    else if (extension == "bindos")
    {
        coords = this->read_bindos_header(filename);
    }
    else if (extension == "minidos")
    {
        coords = this->read_minidos_header(filename);
    }
    else
    {
        throw "Dose file extension not recognised";
    }

    return coords;
}

nlohmann::json RecalcControlPoint::write_statistics()
{
    json stats;
    stats["energy"] = this->energy;
    stats["gantry_angle"] = this->gantry_angle;
    stats["arc_index"] = this->arc_index;

    double weight_sum = 0.0;
    int num_apertures = 0;
    std::vector<double> aperture_sizes;

    for (auto &aperture : this->apertures)
    {
        if (aperture->active)
        {
            num_apertures += 1;
            weight_sum += aperture->weight;
        }
    }

    stats["num_apertures"] = num_apertures;
    stats["weight_sum"] = weight_sum;

    return stats;
}

void RecalcControlPoint::output_statistics()
{
    std::cout << this->energy << " MV" << std::endl;
    std::cout << "Gantry angle: " << this->gantry_angle << std::endl;

    double weight_sum = 0.0;
    int num_apertures = 0;
    for (auto &aperture : this->apertures)
    {
        if (aperture->active)
        {
            num_apertures += 1;
            weight_sum += aperture->weight;
        }
    }

    std::cout << "Number of apertures: " << this->apertures.size() << std::endl;
    std::cout << "Sum of weights: " << weight_sum << std::endl;
    std::cout << std::endl;
}

json RecalcControlPoint::to_json()
{
    json cpt_json;
    cpt_json["energy"] = this->energy;
    cpt_json["FFF"] = this->FFF;

    cpt_json["gantry_angle"] = this->gantry_angle;
    cpt_json["couch_angle"] = this->couch_angle;
    cpt_json["col_angle"] = this->collimator_angle;
    cpt_json["arclength_scaling"] = this->arclength_scaling;

    cpt_json["iso"] = this->iso;
    cpt_json["sad"] = this->sad;

    json cpt_apertures;

    for (auto &aperture : this->apertures) {
        cpt_apertures.push_back(aperture->to_json());
    }

    cpt_json["apertures"] = cpt_apertures;

    if (this->particle == electron)
    {
        cpt_json["particle"] = "electron";
    }
    else if (this->particle == photon)
    {
        cpt_json["particle"] = "photon";
    }

    return cpt_json;
}
