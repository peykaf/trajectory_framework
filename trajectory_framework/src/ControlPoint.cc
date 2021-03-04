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
#include "ControlPoint.hh"
#include "SparseDose.hh"
#include "string_utils.hh"
#include "Structure.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t>> sparse_vector;

ControlPoint::ControlPoint()
{
    this->included = false;
    this->beamlets_loaded = false;
    this->dummy = false;
    this->gantry_angle = -9999;
    this->couch_angle = 0.0;
    this->max_leaf_speed = 25.0;        // mm/s defined at MLC plane
    this->max_fluence_rate = 1.6666667; // MU/deg calculated from (10 MU/s) / (6 deg/s)
    this->gantry_speed = 0.1;           // deg / s
    this->arc_index = 1;
    this->transmission_factor = 0.016; // MLC transmission factor
    this->do_mlc_transmission = true;
    this->particle = particle_type::photon;
    this->arclength_scaling = 1.0;
    this->FFF = false;
    this->is_arc = false;
    this->price_norm_type = price_norm_types::directional;
    this->max_overtravel = 157.8; // Varian MLC maximum overtravel in mm
}

ControlPoint::ControlPoint(json &cpt_json) : ControlPoint()
{
    this->price_norm = 1.0;

    this->energy = cpt_json["energy"];
    this->num_beamlet_rows = cpt_json["beamlet_rows"];
    this->num_beamlet_columns = cpt_json["beamlet_columns"];

    this->current_aperture.init_aperture(this->num_beamlet_rows, this->num_beamlet_columns);

    this->row_size = cpt_json["iso_row_size"];
    this->col_size = cpt_json["iso_col_size"];

    // Leaves are not allowed to overtravel more than max_overtravel.
    this->max_l_column = this->num_beamlet_columns / 2 + int(this->max_overtravel / this->col_size);
    this->min_r_column = int((200.0 - this->max_overtravel) / this->col_size) + 1;

    // Use default values if not specified
    this->arc_index = cpt_json.value("arc_index", this->arc_index);
    this->FFF = cpt_json.value("FFF", this->FFF);
    this->gantry_speed = cpt_json.value("gantry_speed", this->gantry_speed);
    this->do_mlc_transmission = cpt_json.value("do_mlc_transmission", this->do_mlc_transmission);

    this->gantry_angle = cpt_json.value("gantry_angle", 0);
    this->couch_angle = cpt_json.value("couch_angle", 0);
    this->collimator_angle = cpt_json.value("col_angle", 0);

    this->sad = cpt_json["sad"];
    this->iso = cpt_json["iso"].get<std::vector<double>>();

    this->beamlet_filenames = cpt_json["beamlets"].get<std::vector<std::string>>();

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

ControlPoint::ControlPoint(json &cpt_json, particle_type ptype) : ControlPoint(cpt_json)
{
    this->particle = ptype;

    // MLC transmission is only taken into account for photons
    if (this->particle == particle_type::electron)
    {
        this->do_mlc_transmission = false;
    }
}

void ControlPoint::set_price_norm(price_norm_types new_type)
{
    this->price_norm_type = new_type;
}

double ControlPoint::calculate_aperture_dose(Aperture &aperture)
{
    if (aperture.dose.size() != this->total_voxels)
    {
        aperture.dose.resize(this->total_voxels);
    }

    std::fill(aperture.dose.begin(), aperture.dose.end(), 0.0);
    for (size_t row_i = 0; row_i < aperture.rows.size(); row_i++)
    {
        size_t beamlet_offset = row_i * this->num_beamlet_columns;

        for (size_t beamlet_i = aperture.rows[row_i].l_bound;
             beamlet_i < aperture.rows[row_i].r_bound;
             beamlet_i++)
        {
            size_t index = beamlet_i + beamlet_offset;

            for (sparse_vector::const_iterator it = this->beamlets[index].grid.begin();
                 it != this->beamlets[index].grid.end();
                 ++it)
            {
                aperture.dose[it.index()] += *it;
            }
        }
    }

    double unused = 0.0;
    return unused;
}

double ControlPoint::calculate_transmission_aperture_dose(Aperture &aperture)
{
    if (aperture.dose.size() != this->total_voxels)
    {
        aperture.dose.resize(this->total_voxels);
    }

    std::fill(aperture.dose.begin(), aperture.dose.end(), 0.0);
    for (size_t row_i = 0; row_i < aperture.rows.size(); row_i++)
    {
        size_t beamlet_offset = row_i * num_beamlet_columns;

        for (size_t beamlet_i = 0; beamlet_i < this->num_beamlet_columns; beamlet_i++)
        {
            size_t index = beamlet_i + beamlet_offset;

            // MLC transmission factor for closed beamlets
            double scaling = this->transmission_factor;
            if (beamlet_i >= aperture.rows[row_i].l_bound && beamlet_i < aperture.rows[row_i].r_bound)
            {
                scaling = 1.0;
            }

            for (sparse_vector::const_iterator it = beamlets[index].grid.begin();
                 it != beamlets[index].grid.end();
                 ++it)
            {
                aperture.dose[it.index()] += *it * scaling;
            }
        }
    }

    double unused = 0.0;
    return unused;
}

void ControlPoint::update_aperture_dose(Aperture *aperture, size_t row, size_t l_bound, size_t r_bound)
{
    size_t beamlet_offset = row * this->num_beamlet_columns;

    for (size_t beamlet_i = aperture->rows[row].l_bound;
         beamlet_i < aperture->rows[row].r_bound;
         beamlet_i++)
    {
        size_t index = beamlet_i + beamlet_offset;

        for (sparse_vector::const_iterator it = beamlets[index].grid.begin();
             it != beamlets[index].grid.end();
             ++it)
        {
            aperture->dose[it.index()] -= (*it / aperture->norm);
        }
    }

    for (size_t beamlet_i = l_bound; beamlet_i < r_bound; beamlet_i++)
    {
        size_t index = beamlet_i + beamlet_offset;

        for (sparse_vector::const_iterator it = beamlets[index].grid.begin();
             it != beamlets[index].grid.end();
             ++it)
        {
            aperture->dose[it.index()] += (*it / aperture->norm);
        }
    }

    aperture->rows[row].l_bound = l_bound;
    aperture->rows[row].r_bound = r_bound;
}

void ControlPoint::remove_inactive_apertures()
{
    this->apertures.erase(std::remove_if(this->apertures.begin(), this->apertures.end(),
                                         [](Aperture *ap) { return !ap->active; }),
                          this->apertures.end());

    if (this->apertures.size() == 0)
        this->included = false;
}

void ControlPoint::calculate_aperture_size()
{
    double ap_size = 0.0;
    for (size_t row_i = 0; row_i < this->current_aperture.rows.size(); row_i++)
    {
        for (size_t beamlet_i = this->current_aperture.rows[row_i].l_bound;
             beamlet_i < this->current_aperture.rows[row_i].r_bound;
             beamlet_i++)
        {
            ap_size += this->col_size * this->row_size;
        }
    }

    this->current_aperture.size = ap_size;
}

void ControlPoint::normalize_beamlets()
{
    // All beamlet doses within this control point are normalized by the
    // maximum dose of all beamlets.
    double max_dose = 0;
    for (auto &beamlet : this->beamlets)
    {
        if (beamlet.loaded && beamlet.max_dose > max_dose)
        {
            max_dose = beamlet.max_dose;
        }
    }

    this->global_weight = max_dose;
    double norm = 1.0 / max_dose;
    for (auto &beamlet : this->beamlets)
    {
        if (beamlet.loaded) {
            //beamlet.normalize_dose(norm);
            if (this->particle == particle_type::electron && beamlet.max_dose > 0) {
                beamlet.normalize_dose(norm);
            }
        }
    }
}

void ControlPoint::read_3ddose_header(std::string filename)
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

void ControlPoint::read_minidos_header(std::string filename)
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

void ControlPoint::read_bindos_header(std::string filename)
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

void ControlPoint::read_header(std::string filename)
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

void ControlPoint::load_beamlets(std::vector<std::string> beamlet_names)
{
    this->beamlet_filenames = beamlet_names;
    this->load_beamlets();
}

void ControlPoint::load_beamlets()
{
    if (this->particle == particle_type::electron)
    {
        std::cout << "Loading " << this->energy << " MeV electron beamlets at " << this->gantry_angle << " deg" << std::endl;
    }
    else
    {
        std::cout << "Loading " << this->energy << " MV photon beamlets at " << this->gantry_angle << " deg" << std::endl;
    }

    for (auto beamlet_filename : this->beamlet_filenames)
    {
        SparseDose dose;
        if (beamlet_filename != "Inactive")
        {
            dose.from_file(beamlet_filename);
        }
        this->beamlets.push_back(dose);
    }

    this->normalize_beamlets();
    this->beamlets_loaded = true;
}

void ControlPoint::free_beamlets()
{
    if (this->beamlets_loaded)
    {
        for (auto &beamlet : this->beamlets)
        {
            sparse_vector emptyvect;
            beamlet.grid.swap(emptyvect);
        }
        this->beamlets.clear();
        this->beamlets.shrink_to_fit();
    }
    this->beamlets_loaded = false;
}

double ControlPoint::calculate_col_price(size_t row, size_t col, std::vector<double> &lag_mults)
{
    double price = 0.0;
    size_t beamlet_offset = row * this->num_beamlet_columns + col;
    if (!this->beamlets[beamlet_offset].loaded)
        return price;

    for (sparse_vector::const_iterator it = this->beamlets[beamlet_offset].grid.begin();
         it != this->beamlets[beamlet_offset].grid.end();
         ++it)
    {
        price += lag_mults[it.index()] * (*it);
    }

    return price;
}

bool ControlPoint::is_row_deliverable(ControlPoint *prev_cpt, ControlPoint *next_cpt,
                                      size_t row, size_t l_bound, size_t r_bound)
{
    double t_prev = fabs(this->gantry_angle - prev_cpt->gantry_angle);
    double t_next = fabs(this->gantry_angle - next_cpt->gantry_angle);

    double max_dist_prev = max_leaf_speed * t_prev / this->gantry_speed;
    double max_dist_next = max_leaf_speed * t_next / this->gantry_speed;

    int max_unit_prev = max_dist_prev / this->col_size;
    int max_unit_next = max_dist_next / this->col_size;

    // If the previous/next control points belong to another arc, then
    // there are no machine restrictions on the aperture of the current
    // cpt with respect to the previous/next cpts.
    if (this->arc_index != prev_cpt->arc_index)
        max_unit_prev = 9999;
    if (this->arc_index != next_cpt->arc_index)
        max_unit_next = 9999;

    size_t prev_l, prev_r, next_l, next_r;
    // The dummy control point has no aperture defined
    if (prev_cpt->apertures.size() > 0 && next_cpt->apertures.size() > 0)
    {
        prev_l = prev_cpt->apertures[0]->rows[row].l_bound;
        prev_r = prev_cpt->apertures[0]->rows[row].r_bound;
        next_l = next_cpt->apertures[0]->rows[row].l_bound;
        next_r = next_cpt->apertures[0]->rows[row].r_bound;
    }
    else
    {
        prev_l = 0;
        prev_r = this->num_beamlet_columns;
        next_l = 0;
        next_r = this->num_beamlet_columns;
    }

    return (abs(l_bound - prev_l) <= max_unit_prev && abs(l_bound - next_l) <= max_unit_prev && abs(r_bound - prev_r) <= max_unit_next && abs(r_bound - next_r) <= max_unit_next);
}

void ControlPoint::price_constrained_aperture(std::vector<double> &lag_mults,
                                              std::vector<Structure *> structures,
                                              ControlPoint *prev_cpt,
                                              ControlPoint *next_cpt)
{
    if (this->beamlets_loaded == false)
        this->load_beamlets();

    double t_prev = fabs(this->gantry_angle - prev_cpt->gantry_angle);
    double t_next = fabs(this->gantry_angle - next_cpt->gantry_angle);

    double max_dist_prev = max_leaf_speed * t_prev / this->gantry_speed;
    double max_dist_next = max_leaf_speed * t_next / this->gantry_speed;

    int max_unit_prev = max_dist_prev / this->col_size;
    int max_unit_next = max_dist_next / this->col_size;
    std::cout << "max unit prev " << max_unit_prev << std::endl;
    std::cout << "max unit next " << max_unit_next << std::endl;

    // If the previous/next control points belong to another arc, then
    // there are no machine restrictions on the aperture of the current
    // cpt with respect to the previous/next cpts.
    if (this->arc_index != prev_cpt->arc_index)
        max_unit_prev = 9999;
    if (this->arc_index != next_cpt->arc_index)
        max_unit_next = 9999;

    this->latest_price = 0.0;

    for (size_t row = 0; row < this->current_aperture.rows.size(); row++)
    {
        int prev_l, prev_r, next_l, next_r;
        int unit_prev, unit_next;
        double current_val = 0.0;
        double min_val = 0.0;
        double optimal_val = 0.0;
        int optimal_cl = 0;
        int optimal_cr = 0;
        int cl = 0;

        // The dummy control point has no aperture defined
        if (prev_cpt->current_aperture.rows.size() > 0 && next_cpt->current_aperture.rows.size() > 0)
        {
            prev_l = prev_cpt->current_aperture.rows[row].l_bound;
            prev_r = prev_cpt->current_aperture.rows[row].r_bound;
            next_l = next_cpt->current_aperture.rows[row].l_bound;
            next_r = next_cpt->current_aperture.rows[row].r_bound;
            unit_prev = max_unit_prev;
            unit_next = max_unit_next;

            // If the leaves on both the previous and next control point are closed,
            // place no restriction on leaf movement for this row.
            if (prev_r == prev_l || next_r == next_l) {
                unit_prev = 9999;
                unit_next = 9999;
            }
        }
        else
        {
            prev_l = 0;
            prev_r = this->num_beamlet_columns;
            next_l = 0;
            next_r = this->num_beamlet_columns;
            unit_prev = 9999;
            unit_next = 9999;
        }

        // fuck me
        // Find the allowable range of left and right MLC leaf positions.
        int min_prev_l = std::max(0, prev_l - unit_prev);
        int min_next_l = std::max(0, next_l - unit_next);
        int min_cl = std::max(min_prev_l, min_next_l);

        int max_prev_cl = std::min(prev_l + unit_prev, int(this->num_beamlet_columns) - 1);
        int max_next_cl = std::min(next_l + unit_next, int(this->num_beamlet_columns) - 1);
        int max_cl = std::min(max_prev_cl, max_next_cl);

        int min_prev_r = std::max(1, prev_r - unit_prev);
        int min_next_r = std::max(1, next_r - unit_next);
        int min_cr = std::max(min_prev_r, min_next_r);

        int max_prev_cr = std::min(prev_r + unit_prev, int(this->num_beamlet_columns));
        int max_next_cr = std::min(next_r + unit_next, int(this->num_beamlet_columns));
        int max_cr = std::min(max_prev_cr, max_next_cr);

        // Calculate the prices for the beamlets between min_cl and min_cr to figure
        // out what the initial value of cl is going to be.
        std::vector<double> initial_prices;
        for (int cl = min_cl; cl < min_cr; cl++)
        {
            initial_prices.push_back(calculate_col_price(row, cl, lag_mults));
        }

        // Sometimes max_cl will be less than min_cr, we want to catch this.
        // This will happen if leaf constraints force us to have a sub-optimal open row.
        int opt_max_cl = std::min(max_cl, min_cr);

        int opt_init_cl = min_cl;
        // We may be forced to start with a negative price due to leaf constraints.
        // Don't want to initialize this to 0.
        double opt_init_price = -1e10;

        if (opt_max_cl - min_cl < 0) {
            opt_init_price = 0.0;
        }

        // <= instead of < because we can have completely closed initial state, ie cl = min_cr
        for (int i = 0; i <= (opt_max_cl - min_cl); i++)
        {
            double running_price = 0.0;
            for (int j = i; j < initial_prices.size(); j++)
            {
                running_price += initial_prices[j];
            }

            if (running_price >= opt_init_price)
            {
                opt_init_cl = min_cl + i;
                opt_init_price = running_price;
            }

            if (running_price <= min_val)
            {
                min_val = running_price;
            }
        }

        cl = opt_init_cl;
        optimal_cl = opt_init_cl;
        optimal_cr = min_cr;
        optimal_val = opt_init_price;
        current_val = opt_init_price;
        if (current_val <= min_val)
        {
            min_val = current_val;
        }

        for (int c_right = min_cr + 1; c_right <= max_cr; c_right++)
        {
            // As long as we're allowed to move our right leaf
            current_val += calculate_col_price(row, c_right - 1, lag_mults);

            // If the price was lowered by this beamlet, then it's a bad beamlet, skip over it.
            if (current_val <= min_val)
            {
                min_val = current_val;
                // We would like to move our left leaf forward, but we may not be able to.
                if (abs(cl - prev_l) < unit_prev && abs(cl - next_l) < unit_next)
                {
                    cl = c_right;
                }
            }

            if ((current_val - min_val) > optimal_val)
            {
                optimal_val = current_val - min_val;
                optimal_cl = cl;
                optimal_cr = c_right;
            }
        }

        // If leaves violate the overtravel constraints, then the row must be closed.
        if (optimal_cl > this->max_l_column || optimal_cr < this->min_r_column)
        {
            optimal_val = 0.0;
            optimal_cl = this->num_beamlet_columns / 2;
            optimal_cr = this->num_beamlet_columns / 2;
        }

        if (optimal_cl == optimal_cr) {
            // Reposition closed leaves in the center of the field.
            int average_cl = (prev_l + next_l) / 2;
            int average_cr = (prev_r + next_r) / 2;
            int avg_closed = (average_cl + average_cr) / 2;

            optimal_cl = avg_closed;
            optimal_cr = avg_closed;
        }
        this->latest_price += optimal_val;
        this->current_aperture.rows[row].l_bound = optimal_cl;
        this->current_aperture.rows[row].r_bound = optimal_cr;
    }

    if (this->do_mlc_transmission)
    {
        this->calculate_transmission_aperture_dose(this->current_aperture);
    }
    else
    {
        this->calculate_aperture_dose(this->current_aperture);
    }

    double price_norm = 0.0;
    if (this->price_norm_type == price_norm_types::directional)
    {
        price_norm = this->directional_price_norm(lag_mults, structures);
    }
    else if (this->price_norm_type == price_norm_types::dmax)
    {
        price_norm = this->dmax_price_norm(lag_mults, structures);
    }

    if (price_norm > 0.0)
    {
        this->latest_price = latest_price / price_norm;
    }

    std::cout << this->energy << " MV: " << this->latest_price << std::endl;
    //this->current_aperture.normalize_dose();
}

void ControlPoint::price_aperture(std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    this->latest_price = 0.0;
    if (this->beamlets_loaded == false)
        this->load_beamlets();

    for (size_t row = 0; row < this->current_aperture.rows.size(); row++)
    {
        double current_val = 0.0;
        double min_val = 0.0;
        double optimal_val = 0.0;
        int optimal_cl = 0;
        int optimal_cr = 0;
        int cl = 0;

        for (size_t c_right = 1; c_right <= this->num_beamlet_columns; c_right++)
        {
            current_val += calculate_col_price(row, c_right - 1, lag_mults);
            if (current_val <= min_val)
            {
                min_val = current_val;
                cl = c_right;
            }

            if ((current_val - min_val) > optimal_val)
            {
                optimal_val = current_val - min_val;
                optimal_cl = cl;
                optimal_cr = c_right;
            }
        }

        // If leaves violate the overtravel constraints, then the row must be closed.
        if (optimal_cl > this->max_l_column || optimal_cr < this->min_r_column)
        {
            optimal_val = 0.0;
            optimal_cl = 0;
            optimal_cr = 0;
        }
        this->latest_price += optimal_val;
        this->current_aperture.rows[row].l_bound = optimal_cl;
        this->current_aperture.rows[row].r_bound = optimal_cr;
    }

    if (this->do_mlc_transmission)
    {
        this->calculate_transmission_aperture_dose(this->current_aperture);
    }
    else
    {
        this->calculate_aperture_dose(this->current_aperture);
    }

    double price_norm = 0.0;
    if (this->price_norm_type == price_norm_types::directional)
    {
        price_norm = this->directional_price_norm(lag_mults, structures);
    }
    else if (this->price_norm_type == price_norm_types::dmax)
    {
        price_norm = this->dmax_price_norm(lag_mults, structures);
    }

    if (price_norm > 0.0)
    {
        this->latest_price = latest_price / price_norm;
    }

    //this->current_aperture.normalize_dose();

    std::cout << this->energy << " MV: " << this->latest_price << std::endl;
}

double ControlPoint::dmax_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    // Dmax rule
    double dmax = 0.0;
    for (auto &structure : structures)
    {
        for (size_t i = 0; i < structure->masked_dose.size(); i++)
        {
            size_t vox_num = structure->masked_dose[i].first;
            if (this->current_aperture.dose[vox_num] > dmax)
                dmax = this->current_aperture.dose[vox_num];
        }
    }

    return dmax;
}

double ControlPoint::average_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    return this->current_aperture.average_dose();
}

double ControlPoint::directional_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures)
{
    double norm;

    // Directional norm rule
    double inner_dose = 0.0;
    double inner_lag = 0.0;
    for (auto &structure : structures)
    {
        for (size_t i = 0; i < structure->masked_dose.size(); i++)
        {
            size_t vox_num = structure->masked_dose[i].first;
            if (std::abs(lag_mults[vox_num]) > 0)
            {
                inner_dose += this->current_aperture.dose[vox_num] * this->current_aperture.dose[vox_num];
                inner_lag += lag_mults[vox_num] * lag_mults[vox_num];
            }
        }
    }

    norm = std::sqrt(inner_dose * inner_lag);
    return norm;
}

Aperture *ControlPoint::save_current_aperture()
{
    //this->calculate_aperture_dose(this->current_aperture);
    this->calculate_aperture_size();
    Aperture *new_aperture = new Aperture();
    std::copy(this->current_aperture.dose.begin(),
              this->current_aperture.dose.end(),
              std::back_inserter(new_aperture->dose));
    std::copy(this->current_aperture.rows.begin(),
              this->current_aperture.rows.end(),
              std::back_inserter(new_aperture->rows));
    new_aperture->weight = 0.0;
    new_aperture->active = true;
    new_aperture->size = this->current_aperture.size;
    new_aperture->num_beamlet_columns = this->current_aperture.num_beamlet_columns;
    this->apertures.push_back(new_aperture);
    return new_aperture;
}

nlohmann::json ControlPoint::write_statistics()
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
            aperture_sizes.push_back(aperture->size);
        }
    }

    stats["num_apertures"] = num_apertures;
    stats["weight_sum"] = weight_sum;
    stats["aperture_sizes"] = aperture_sizes;

    return stats;
}

void ControlPoint::output_statistics()
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

json ControlPoint::to_json()
{
    json cpt_json;
    cpt_json["energy"] = this->energy;
    cpt_json["FFF"] = this->FFF;
    cpt_json["beamlet_rows"] = this->num_beamlet_rows;
    cpt_json["beamlet_columns"] = this->num_beamlet_columns;
    cpt_json["iso_row_size"] = this->row_size;
    cpt_json["iso_col_size"] = this->col_size;

    cpt_json["gantry_angle"] = this->gantry_angle;
    cpt_json["couch_angle"] = this->couch_angle;
    cpt_json["col_angle"] = this->collimator_angle;
    cpt_json["arclength_scaling"] = this->arclength_scaling;

    cpt_json["iso"] = this->iso;
    cpt_json["sad"] = this->sad;
    cpt_json["is_arc"] = this->is_arc;

    //double relative_weight = max_global_weight / this->global_weight;
    double x_field_size = this->num_beamlet_columns * this->col_size;

    json cpt_apertures;

    // Epsilon to determine whether leaves are closed.
    double epsilon = 1e-3;
    for (size_t ap_i = 0; ap_i < this->apertures.size(); ap_i++)
    {
        if (this->apertures[ap_i]->active)
        {
            std::vector<std::vector<double>> ap_vector;
            for (size_t i = 0; i < this->apertures[ap_i]->rows.size(); i++)
            {
                double x_neg = this->apertures[ap_i]->rows[i].l_bound * this->col_size - 0.5 * x_field_size;
                double x_pos = this->apertures[ap_i]->rows[i].r_bound * this->col_size - 0.5 * x_field_size;
                std::vector<double> ap_row = {x_neg, x_pos};
                if (std::abs(x_pos - x_neg) < epsilon)
                {
                    ap_row[0] = 0.0;
                    ap_row[1] = 0.0;
                }
                ap_vector.push_back(ap_row);
            }

            json ap_json;
            ap_json["rows"] = ap_vector;
            ap_json["weight"] = this->apertures[ap_i]->weight;
            cpt_apertures.push_back(ap_json);
        }
    }

    cpt_json["apertures"] = cpt_apertures;

    if (this->particle == particle_type::electron)
    {
        cpt_json["particle"] = "electron";
    }
    else if (this->particle == particle_type::photon)
    {
        cpt_json["particle"] = "photon";
    }

    return cpt_json;
}
