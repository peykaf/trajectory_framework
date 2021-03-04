// Copyright 2018 Marc-Andre Renaud
#include <iostream>
#include <vector>

#include "json.hh"
#include "SparseDose.hh"
#include "RobustControlPoint.hh"
#include "ControlPoint.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t>> sparse_vector;

RobustControlPoint::RobustControlPoint(json &cpt_json) : ControlPoint(cpt_json)
{
    // Robust beamlet filenames are the only RobustControlPoint-specific feature.
    for (auto filenames : cpt_json["robust_beamlets"])
    {
        std::vector<std::string> parsed_filenames = filenames.get<std::vector<std::string>>();
        this->robust_filenames.push_back(parsed_filenames);
    }

    this->robust_apertures.resize(this->get_num_scenarios() - 1);
}

RobustControlPoint::RobustControlPoint(json &cpt_json, particle_type ptype) : RobustControlPoint(cpt_json)
{
    this->particle = ptype;

    // MLC transmission is only taken into account for photons
    if (this->particle == particle_type::electron)
    {
        this->do_mlc_transmission = false;
    }
}

void RobustControlPoint::load_beamlets()
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
        if (beamlet_filename != "Inactive")
        {
            this->beamlets.push_back(SparseDose(beamlet_filename));
        }
        else
        {
            this->beamlets.push_back(SparseDose());
        }
    }
    double norm = this->normalize_beamlets(this->beamlets);

    for (auto beamlet_set : this->robust_filenames)
    {
        std::vector<SparseDose> robust_set;
        for (auto beamlet_filename : beamlet_set)
        {
            if (beamlet_filename != "Inactive")
            {
                robust_set.push_back(SparseDose(beamlet_filename));
            }
            else
            {
                robust_set.push_back(SparseDose());
            }
        }
        this->normalize_beamlets(robust_set, norm);
        this->robust_beamlets.push_back(robust_set);
    }

    this->beamlets_loaded = true;
}

size_t RobustControlPoint::get_num_scenarios()
{
    return this->robust_filenames.size() + 1;
}

double RobustControlPoint::calculate_col_price(size_t row, size_t col, std::vector<double> &lag_mults, size_t scenario)
{
    double price = 0.0;
    size_t beamlet_offset = row * this->num_beamlet_columns + col;
    if (!this->beamlets[beamlet_offset].loaded)
        return price;

    if (scenario == 0)
    {
        for (sparse_vector::const_iterator it = this->beamlets[beamlet_offset].grid.begin();
             it != this->beamlets[beamlet_offset].grid.end();
             ++it)
        {
            price += lag_mults[it.index()] * (*it);
        }
    }
    else
    {
        for (sparse_vector::const_iterator it = this->robust_beamlets[scenario - 1][beamlet_offset].grid.begin();
             it != this->robust_beamlets[scenario - 1][beamlet_offset].grid.end();
             ++it)
        {
            price += lag_mults[it.index()] * (*it);
        }
    }

    return price;
}

void RobustControlPoint::price_aperture(std::vector<std::vector<double>> &lag_mults, std::vector<double> &kkt_mults, std::vector<std::vector<Structure *>> structure_sets)
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
            double col_val = 0.0;
            for (size_t sc_i = 0; sc_i < lag_mults.size(); sc_i++)
            {
                if (kkt_mults[sc_i] > 0.0)
                    col_val += kkt_mults[sc_i] * this->calculate_col_price(row, c_right - 1, lag_mults[sc_i], sc_i);
            }

            current_val += col_val;
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
        if (optimal_cl > this->max_l_column || optimal_cr < this->min_r_column) {
            optimal_val = 0.0;
            optimal_cl = 0;
            optimal_cr = 0;
        }

        this->latest_price += optimal_val;
        this->current_aperture.rows[row].l_bound = optimal_cl;
        this->current_aperture.rows[row].r_bound = optimal_cr;
    }

    // At this stage, aperture dose only needs to be calculated for nominal scenario.
    if (this->do_mlc_transmission)
    {
        this->calculate_transmission_aperture_dose(&this->current_aperture, this->beamlets);
    }
    else
    {
        this->calculate_aperture_dose(&this->current_aperture, this->beamlets);
    }

    std::vector<Structure *> structures = structure_sets[0];
    std::vector<double> lags = lag_mults[0];

    double price_norm = 0.0;
    if (this->price_norm_type == price_norm_types::directional)
    {
        price_norm = this->directional_price_norm(lags, structures);
    }
    else if (this->price_norm_type == price_norm_types::dmax)
    {
        price_norm = this->dmax_price_norm(lags, structures);
    }

    if (price_norm > 0.0)
    {
        this->latest_price = latest_price / price_norm;
    }

    std::cout << this->energy << " MV: " << this->latest_price << std::endl;
    //this->current_aperture.normalize_dose();
}

void RobustControlPoint::price_constrained_aperture(std::vector<std::vector<double>> &lag_mults, std::vector<double> &kkt_mults, std::vector<std::vector<Structure *>> structure_sets,
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
        }
        else
        {
            prev_l = 0;
            prev_r = this->num_beamlet_columns;
            next_l = 0;
            next_r = this->num_beamlet_columns;
        }

        // fuck me
        // Find the allowable range of left and right MLC leaf positions.
        int min_prev_l = std::max(0, prev_l - max_unit_prev);
        int min_next_l = std::max(0, next_l - max_unit_next);
        int min_cl = std::max(min_prev_l, min_next_l);

        int max_prev_cl = std::min(prev_l + max_unit_prev, int(this->num_beamlet_columns) - 1);
        int max_next_cl = std::min(next_l + max_unit_next, int(this->num_beamlet_columns) - 1);
        int max_cl = std::min(max_prev_cl, max_next_cl);

        int min_prev_r = std::max(1, prev_r - max_unit_prev);
        int min_next_r = std::max(1, next_r - max_unit_next);
        int min_cr = std::max(min_prev_r, min_next_r);

        int max_prev_cr = std::min(prev_r + max_unit_prev, int(this->num_beamlet_columns));
        int max_next_cr = std::min(next_r + max_unit_next, int(this->num_beamlet_columns));
        int max_cr = std::min(max_prev_cr, max_next_cr);

        // Calculate the prices for the beamlets between min_cl and min_cr to figure
        // out what the initial value of cl is going to be.
        std::vector<double> initial_prices;
        for (int cl = min_cl; cl < min_cr; cl++)
        {
            double col_val = 0.0;
            for (size_t sc_i = 0; sc_i < lag_mults.size(); sc_i++)
            {
                col_val += kkt_mults[sc_i] * this->calculate_col_price(row, cl, lag_mults[sc_i], sc_i);
            }

            initial_prices.push_back(col_val);
        }

        // Sometimes max_cl will be less than min_cr, we want to catch this.
        // This will happen if leaf constraints force us to have a sub-optimal open row.
        int opt_max_cl = std::min(min_cr, max_cl);

        int opt_init_cl = min_cl;
        // We may be forced to start with a negative price due to leaf constraints.
        // Don't want to initialize this to 0.
        double opt_init_price = -1e10;

        // <= instead of < because we can have completely closed initial state, ie cl = min_cr
        for (int i = 0; i <= opt_max_cl; i++)
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
        optimal_val = opt_init_price;
        current_val = opt_init_price;
        if (current_val <= min_val)
        {
            min_val = current_val;
        }

        for (int c_right = min_cr + 1; c_right <= max_cr; c_right++)
        {
            double col_val = 0.0;
            for (size_t sc_i = 0; sc_i < lag_mults.size(); sc_i++)
            {
                col_val += kkt_mults[sc_i] * this->calculate_col_price(row, c_right - 1, lag_mults[sc_i], sc_i);
            }

            // As long as we're allowed to move our right leaf
            current_val += col_val;

            // If the price was lowered by this beamlet, then it's a bad beamlet, skip over it.
            if (current_val <= min_val)
            {
                min_val = current_val;
                // We would like to move our left leaf forward, but we may not be able to.
                if (abs(cl - prev_l) < max_unit_prev && abs(cl - next_l) < max_unit_next)
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
        if (optimal_cl > this->max_l_column || optimal_cr < this->min_r_column) {
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
        this->calculate_transmission_aperture_dose(&this->current_aperture, this->beamlets);
    }
    else
    {
        this->calculate_aperture_dose(&this->current_aperture, this->beamlets);
    }

    std::vector<Structure *> structures = structure_sets[0];
    std::vector<double> lags = lag_mults[0];

    double price_norm = 0.0;
    if (this->price_norm_type == price_norm_types::directional)
    {
        price_norm = this->directional_price_norm(lags, structures);
    }
    else if (this->price_norm_type == price_norm_types::dmax)
    {
        price_norm = this->dmax_price_norm(lags, structures);
    }

    if (price_norm > 0.0)
    {
        this->latest_price = latest_price / price_norm;
    }

    //this->current_aperture.normalize_dose();

    std::cout << this->energy << " MV: " << this->latest_price << std::endl;
}

double RobustControlPoint::normalize_beamlets(std::vector<SparseDose> &beamlets)
{
    // All beamlet doses within this control point are normalized by the
    // maximum dose of all beamlets.
    double max_dose = 0.0;
    for (auto &beamlet : beamlets)
    {
        if (beamlet.loaded && beamlet.max_dose > max_dose)
        {
            max_dose = beamlet.max_dose;
        }
    }

    double norm = 1.0 / max_dose;
    for (auto &beamlet : beamlets)
    {
        if (beamlet.loaded)
            beamlet.normalize_dose(norm);
    }

    return norm;
}

void RobustControlPoint::normalize_beamlets(std::vector<SparseDose> &beamlets, double norm)
{
    for (auto &beamlet : beamlets)
    {
        if (beamlet.loaded)
            beamlet.normalize_dose(norm);
    }
}

std::vector<Aperture *> RobustControlPoint::save_robust_aperture()
{
    std::vector<Aperture *> new_apertures;

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
    new_apertures.push_back(new_aperture);

    for (size_t i = 0; i < this->robust_beamlets.size(); i++)
    {
        Aperture *new_aperture = new Aperture();
        std::copy(this->current_aperture.rows.begin(),
                  this->current_aperture.rows.end(),
                  std::back_inserter(new_aperture->rows));
        new_aperture->num_beamlet_columns = this->current_aperture.num_beamlet_columns;
        new_aperture->weight = 0.0;
        new_aperture->size = this->current_aperture.size;
        this->calculate_aperture_dose(new_aperture, this->robust_beamlets[i]);
        new_aperture->active = true;
        //new_aperture->normalize_dose();

        this->robust_apertures[i].push_back(new_aperture);
        new_apertures.push_back(new_aperture);
    }

    return new_apertures;
}

void RobustControlPoint::remove_inactive_apertures()
{
    this->apertures.erase(std::remove_if(this->apertures.begin(), this->apertures.end(),
                                         [](Aperture *ap) { return !(ap->active); }),
                          this->apertures.end());

    if (this->apertures.size() == 0)
        this->included = false;

    for (size_t set_i = 0; set_i < this->robust_apertures.size(); set_i++)
    {
        this->robust_apertures[set_i].erase(std::remove_if(this->robust_apertures[set_i].begin(), this->robust_apertures[set_i].end(),
                                                           [](Aperture *ap) { return !(ap->active); }),
                                            this->robust_apertures[set_i].end());
    }
}

double RobustControlPoint::calculate_aperture_dose(Aperture *aperture, const std::vector<SparseDose> &beamlet_set)
{
    if (aperture->dose.size() != this->total_voxels)
    {
        aperture->dose.resize(this->total_voxels);
    }

    std::fill(aperture->dose.begin(), aperture->dose.end(), 0.0);
    for (size_t row_i = 0; row_i < aperture->rows.size(); row_i++)
    {
        // Linearized beamlet offset
        size_t beamlet_offset = row_i * this->num_beamlet_columns;

        for (size_t beamlet_i = aperture->rows[row_i].l_bound;
             beamlet_i < aperture->rows[row_i].r_bound;
             beamlet_i++)
        {
            size_t index = beamlet_i + beamlet_offset;

            for (sparse_vector::const_iterator it = beamlet_set[index].grid.begin();
                 it != beamlet_set[index].grid.end();
                 ++it)
            {
                aperture->dose[it.index()] += *it;
            }
        }
    }

    double max_dose = 0.0;
    return max_dose;
}

double RobustControlPoint::calculate_transmission_aperture_dose(Aperture *aperture, const std::vector<SparseDose> &beamlet_set)
{
    if (aperture->dose.size() != this->total_voxels)
    {
        aperture->dose.resize(this->total_voxels);
    }

    std::fill(aperture->dose.begin(), aperture->dose.end(), 0.0);
    for (size_t row_i = 0; row_i < aperture->rows.size(); row_i++)
    {
        size_t beamlet_offset = row_i * this->num_beamlet_columns;

        for (size_t beamlet_i = 0; beamlet_i < this->num_beamlet_columns; beamlet_i++)
        {
            size_t index = beamlet_i + beamlet_offset;

            // MLC transmission factor for closed beamlets
            double scaling = this->transmission_factor;
            if (beamlet_i >= aperture->rows[row_i].l_bound && beamlet_i < aperture->rows[row_i].r_bound)
            {
                scaling = 1.0;
            }

            for (sparse_vector::const_iterator it = beamlet_set[index].grid.begin();
                 it != beamlet_set[index].grid.end();
                 ++it)
            {
                aperture->dose[it.index()] += *it * scaling;
            }
        }
    }

    double max_dose = 0.0;
    return max_dose;
}
