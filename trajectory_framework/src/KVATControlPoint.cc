// Copyright 2015 Marc-André Renaud
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "json.hh"
#include "KVATControlPoint.hh"
#include "KVATBeamlet.hh"
#include "SparseDose.hh"
#include "Structure.hh"
#include "string_utils.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t> > sparse_vector;

KVATControlPoint::KVATControlPoint()
{
    this->included = false;
    this->beamlets_loaded = false;
    this->gantry_angle = -9999;
    this->couch_angle = 0.0;
    this->max_fluence_rate = 1.6666667; // MU/deg calculated from (10 MU/s) / (6 deg/s)
    this->gantry_speed = 0.1;           // deg / s
    this->arclength_scaling = 1.0;
    this->arc_index = 1;
}

KVATControlPoint::KVATControlPoint(json &cpt_json) : KVATControlPoint()
{
    this->price_norm = 1.0;

    this->energy = cpt_json["energy"];
    this->num_beamlet_rows = 1; // Only one row for now

    if (cpt_json.find("arc_index") != cpt_json.end())
    {
        this->arc_index = cpt_json["arc_index"];
    }

    this->gantry_angle = cpt_json["gantry_angle"];
    this->couch_angle = cpt_json["couch_angle"];
    this->collimator_angle = cpt_json["col_angle"];

    this->sad = cpt_json["sad"];
    this->iso = cpt_json["iso"].get<std::vector<double>>();

    this->beamlet_filenames = cpt_json["beamlets"].get<std::vector<std::string>>();

    this->num_beamlet_columns = this->beamlet_filenames.size();
    int first_active = -1;
    for (size_t i = 0; i < this->beamlet_filenames.size(); i++)
    {
        if (this->beamlet_filenames[i] != "Inactive")
        {
            first_active = i;
            break;
        }
    }

    this->read_header(this->beamlet_filenames[first_active]);
}

void KVATControlPoint::normalize_beamlets()
{
    // Converts from dosxyznrc dose to Gy / min
    this->global_weight = 2 * 0.69 * 6.4218e20 * (100.0 / 30.0) * (1.0 / 100.0);
    double norm = this->global_weight;
    for (auto &beamlet : this->beamlets)
    {
        if (beamlet.loaded)
            beamlet.normalize_dose(norm);
    }
}

void KVATControlPoint::read_3ddose_header(std::string filename)
{
    /*
      Populates the class attributes related to phantom dimensions.
    */

    std::string buffer;
    std::string trimmed_buffer;
    std::vector<std::string> voxel_buffer;

    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer); // number of voxels
    std::stringstream convertor(buffer);
    convertor >> this->num_voxels[0] >> this->num_voxels[1] >> this->num_voxels[2];
    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    std::cout << "Number of voxels: (" << this->num_voxels[0] << ","
              << this->num_voxels[1] << ","
              << this->num_voxels[2] << ")" << std::endl;

    double second_voxel;
    std::getline(beamlet_file, buffer); // x voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[0];
    convertor >> second_voxel;
    std::cout << trimmed_buffer << std::endl;
    this->voxel_size[0] = second_voxel - this->topleft[0];

    std::getline(beamlet_file, buffer); // y voxel coordinates
    trimmed_buffer = trim(buffer);
    convertor.str(trimmed_buffer);
    convertor.clear();
    convertor >> this->topleft[1];
    convertor >> second_voxel;
    this->voxel_size[1] = second_voxel - this->topleft[1];

    std::getline(beamlet_file, buffer); // z voxel coordinates
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

void KVATControlPoint::read_minidos_header(std::string filename)
{
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    float b_topleft[3], b_voxel_size[3];
    beamlet_file.read((char *)this->num_voxels, 3 * sizeof(int));
    beamlet_file.read((char *)b_voxel_size, 3 * sizeof(float));
    beamlet_file.read((char *)b_topleft, 3 * sizeof(float));
    this->voxel_size[0] = b_voxel_size[0];
    this->voxel_size[1] = b_voxel_size[1];
    this->voxel_size[2] = b_voxel_size[2];

    this->topleft[0] = b_topleft[0];
    this->topleft[1] = b_topleft[1];
    this->topleft[2] = b_topleft[2];

    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];
    beamlet_file.close();
}

void KVATControlPoint::read_header(std::string filename)
{
    std::string extension;

    extension = find_extension(filename);
    if (extension == "3ddose")
    {
        this->read_3ddose_header(filename);
    }
    else if (extension == "bindos")
    {
        this->read_bindos_header(filename);
    }
    else if (extension == "minidos")
    {
        this->read_minidos_header(filename);
    }
    else
    {
        throw "Dose file extension not recognised";
    }
}

void KVATControlPoint::read_bindos_header(std::string filename)
{
    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    beamlet_file.read((char *)this->num_voxels, 3 * sizeof(int));

    float *x_voxels = new float[this->num_voxels[0] + 1];
    float *y_voxels = new float[this->num_voxels[1] + 1];
    float *z_voxels = new float[this->num_voxels[2] + 1];

    beamlet_file.read((char *)x_voxels, (this->num_voxels[0] + 1) * sizeof(float));
    beamlet_file.read((char *)y_voxels, (this->num_voxels[1] + 1) * sizeof(float));
    beamlet_file.read((char *)z_voxels, (this->num_voxels[2] + 1) * sizeof(float));

    beamlet_file.close();

    this->voxel_size[0] = x_voxels[1] - x_voxels[0];
    this->voxel_size[1] = y_voxels[1] - y_voxels[0];
    this->voxel_size[2] = z_voxels[1] - z_voxels[0];

    this->topleft[0] = x_voxels[0];
    this->topleft[1] = y_voxels[0];
    this->topleft[2] = z_voxels[0];

    this->total_voxels = this->num_voxels[0] * this->num_voxels[1] * this->num_voxels[2];

    delete[] x_voxels;
    delete[] y_voxels;
    delete[] z_voxels;
}

void KVATControlPoint::load_beamlets(std::vector<std::string> beamlet_names)
{
    this->beamlet_filenames = beamlet_names;
    this->load_beamlets();
}

void KVATControlPoint::load_beamlets()
{
    std::cout << "Loading beamlets for " << this->energy << " kV control point." << std::endl;
    size_t num_beamlets = this->beamlet_filenames.size();

    for (size_t i = 0; i < num_beamlets; i++)
    {
        SparseDose dose;
        if (this->beamlet_filenames[i] != "Inactive")
        {
            dose.from_file(this->beamlet_filenames[i]);
        }
        this->beamlets.push_back(dose);
    }
    this->normalize_beamlets();
    this->beamlets_loaded = true;
}

void KVATControlPoint::free_beamlets()
{
    if (this->beamlets_loaded)
    {
        for (size_t i = 0; i < this->beamlets.size(); i++)
        {
            sparse_vector emptyvect;
            this->beamlets[i].grid.swap(emptyvect);
        }
        this->beamlets.clear();
        this->beamlets.shrink_to_fit();
    }
    this->beamlets_loaded = false;
}

double KVATControlPoint::calculate_beamlet_price(size_t row, size_t beamlet_id, std::vector<double> &lag_mults)
{
    double price = 0.0;
    size_t beamlet_offset = row * this->num_beamlet_columns + beamlet_id;
    if (!this->beamlets[beamlet_offset].loaded)
        return price;

    for (sparse_vector::iterator it = this->beamlets[beamlet_offset].grid.begin();
         it != this->beamlets[beamlet_offset].grid.end();
         ++it)
    {
        price += lag_mults[it.index()] * (*it);
    }

    return price;
}

void KVATControlPoint::price_beamlets(std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    this->latest_price = 0.0;
    if (this->beamlets_loaded == false)
        this->load_beamlets();
    int best_bid = -1;
    double best_price = 0.0;

    for (size_t row = 0; row < this->num_beamlet_rows; row++)
    {
        for (size_t beamlet_id = 0; beamlet_id < this->num_beamlet_columns; beamlet_id++)
        {
            double col_price = calculate_beamlet_price(row, beamlet_id, lag_mults);
            double price_norm = directional_price_norm(row, beamlet_id, lag_mults, structures);
            if (price_norm > 0.0) col_price /= price_norm;

            if (col_price > best_price)
            {
                best_bid = beamlet_id;
                best_price = col_price;
            }
        }
    }

    this->latest_price = best_price;
    this->best_beamlet = best_bid;
}

double KVATControlPoint::directional_price_norm(size_t row, size_t beamlet_id, std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    double norm;

    // Directional norm rule
    double inner_dose = 0.0;
    double inner_lag = 0.0;
    size_t beamlet_offset = row * this->num_beamlet_columns + beamlet_id;
    for (auto &structure : structures)
    {
        for (size_t i = 0; i < structure->masked_dose.size(); i++)
        {
            size_t vox_num = structure->masked_dose[i].first;
            if (std::abs(lag_mults[vox_num]) > 0.0)
            {
                inner_lag += lag_mults[vox_num] * lag_mults[vox_num];
            }
        }
    }

    for (sparse_vector::iterator it = this->beamlets[beamlet_offset].grid.begin();
         it != this->beamlets[beamlet_offset].grid.end();
         ++it)
    {
        if (std::abs(lag_mults[it.index()]) > 0.0) {
            inner_dose += (*it) * (*it);
        }
    }

    norm = std::sqrt(inner_dose * inner_lag);
    return norm;
}

KVATBeamlet *KVATControlPoint::add_best_beamlet()
{
    KVATBeamlet *new_beamlet = new KVATBeamlet();
    new_beamlet->dose.resize(this->total_voxels);
    std::fill(new_beamlet->dose.begin(), new_beamlet->dose.end(), 0.0);

    for (sparse_vector::iterator it = this->beamlets[this->best_beamlet].grid.begin();
         it != this->beamlets[this->best_beamlet].grid.end();
         ++it)
    {
        new_beamlet->dose[it.index()] = (*it);
    }

    new_beamlet->weight = 0.0;
    new_beamlet->active = true;
    new_beamlet->index = this->best_beamlet;
    this->active_beamlets.push_back(new_beamlet);
    return new_beamlet;
}

nlohmann::json KVATControlPoint::write_statistics()
{
    json stats;
    stats["energy"] = this->energy;
    stats["gantry_angle"] = this->gantry_angle;

    stats["num_active_beamlets"] = this->active_beamlets.size();

    return stats;
}

void KVATControlPoint::output_statistics()
{
    std::cout << "Photon control point" << std::endl;
    std::cout << this->energy << " MV" << std::endl;
    std::cout << "Gantry angle: " << this->gantry_angle << std::endl;

    double weight_sum = 0.0;
    int num_apertures = 0;
    for (auto &beamlet : this->active_beamlets)
    {
        if (beamlet->active)
        {
            weight_sum += beamlet->weight;
            num_apertures += 1;
        }
    }

    std::cout << "Number of active_beamlets: " << num_apertures << std::endl;
    std::cout << "Sum of weights: " << weight_sum << std::endl;
    std::cout << std::endl;
}

json KVATControlPoint::to_json()
{
    json cpt_json;
    cpt_json["energy"] = this->energy;
    cpt_json["beamlet_columns"] = this->num_beamlet_columns;

    cpt_json["gantry_angle"] = this->gantry_angle;
    cpt_json["couch_angle"] = this->couch_angle;
    cpt_json["col_angle"] = this->collimator_angle;
    cpt_json["arclength_scaling"] = this->arclength_scaling;

    cpt_json["iso"] = this->iso;
    cpt_json["sad"] = this->sad;

    std::vector<int> beamlet_indices;
    std::vector<double> beamlet_weights;
    for (auto &beamlet : this->active_beamlets)
    {
        if (beamlet->active && beamlet->weight > 0.0)
        {
            beamlet_indices.push_back(beamlet->index);
            beamlet_weights.push_back(beamlet->weight);
        }
    }

    cpt_json["active_beamlets"] = beamlet_indices;
    cpt_json["beamlet_weights"] = beamlet_weights;

    return cpt_json;
}
