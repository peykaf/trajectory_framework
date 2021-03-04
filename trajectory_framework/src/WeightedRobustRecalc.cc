// Copyright 2018 Marc-Andre Renaud
#include <vector>
#include "json.hh"
#include <chrono>
#include <random>
#include <map>
#include "RobustRecalcControlPoint.hh"
#include "WeightedRobustRecalc.hh"
#include "RecalcAperture.hh"

using json = nlohmann::json;

WeightedRobustRecalc::WeightedRobustRecalc(json &json_input)
{
    std::cout << "Weighted robust recalc optimisation!" << std::endl;
    this->input_filename = json_input["name"];
    this->output_modalities = false;
    this->latest_cost = 0.0;
    this->num_weights = 0;

    this->from_json(json_input);
    this->num_weights = this->aperture_sets[0].size();
    this->num_hess_ele = this->num_weights * (this->num_weights + 1) / 2;
}

void WeightedRobustRecalc::from_json(json &json_input)
{
    this->total_voxels = json_input["total_voxels"];
    this->num_scenarios = json_input["num_scenarios"];

    // Output dose distributions per modality
    this->output_modalities = json_input.value("output_modalities", false);

    this->aperture_sets.resize(this->num_scenarios);
    this->structure_sets.resize(this->num_scenarios);

    int num_cpts = json_input["electron_cpts"].size();
    std::cout << "Number of Electron control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["electron_cpts"])
    {
        RobustRecalcControlPoint *cpt = new RobustRecalcControlPoint(cpt_json, RobustRecalcControlPoint::electron);
        this->electron_cpts.push_back(cpt);
    }

    num_cpts = json_input["photon_cpts"].size();
    std::cout << "Number of Photon control points: " << num_cpts << std::endl;

    for (auto &cpt_json : json_input["photon_cpts"])
    {
        RobustRecalcControlPoint *cpt = new RobustRecalcControlPoint(cpt_json, RobustRecalcControlPoint::photon);
        this->photon_cpts.push_back(cpt);
    }

    std::vector<double> sc_weights(this->num_scenarios, 1.0);
    if (json_input.find("scenario_weights") != json_input.end())
    {
        sc_weights = json_input["scenario_weights"].get<std::vector<double>>();
    }

    double sc_sum = 0.0;
    for (size_t i = 0; i < sc_weights.size(); i++)
    {
        sc_sum += sc_weights[i];
    }
    for (size_t i = 0; i < sc_weights.size(); i++)
    {
        sc_weights[i] /= sc_sum;
    }

    this->scenario_weights = sc_weights;

    std::cout << "Number of structures: " << json_input["structures"].size() << std::endl;
    for (size_t i = 0; i < this->num_scenarios; i++)
    {
        for (auto &struct_json : json_input["structures"])
        {
            bool include_robust = struct_json.value("robust", true);
            // if include_robust is false, only optimise the nominal scenario
            if (include_robust || i == 0)
            {
                Structure *new_struct = new Structure(struct_json);
                this->structure_sets[i].push_back(new_struct);
            }
        }
    }

    std::cout << std::endl
              << "SCENARIO WEIGHTS: " << std::endl;
    for (size_t i = 0; i < this->scenario_weights.size(); i++)
    {
        double wt = this->scenario_weights[i];
        std::cout << wt;
        if (i != this->scenario_weights.size() - 1)
            std::cout << ", ";
    }
    std::cout << std::endl
              << std::endl;

    // Output constraints at the end of the input file reading.
    for (auto &structure : this->structure_sets[0])
    {
        structure->output_constraints();
    }
    std::cout << std::endl;

    this->init_apertures();

    this->start_time = std::chrono::steady_clock::now();
}

void WeightedRobustRecalc::init_apertures() {
    for (auto &cpt : this->photon_cpts) {
        for (auto &aperture : cpt->apertures) {
            this->aperture_sets[0].push_back(aperture);
        }

        for (size_t sc_i = 0; sc_i < cpt->robust_apertures.size(); sc_i++) {
            for (auto &aperture : cpt->robust_apertures[sc_i]) {
                this->aperture_sets[sc_i+1].push_back(aperture);
            }
        }
    }

    for (auto &cpt : this->electron_cpts) {
        for (auto &aperture : cpt->apertures) {
            this->aperture_sets[0].push_back(aperture);
        }

        for (size_t sc_i = 0; sc_i < cpt->robust_apertures.size(); sc_i++) {
            for (auto &aperture : cpt->robust_apertures[sc_i]) {
                this->aperture_sets[sc_i+1].push_back(aperture);
            }
        }
    }
}

void WeightedRobustRecalc::get_starting_point(double *weights)
{
    for (size_t i = 0; i < this->aperture_sets[0].size(); i++)
    {
        weights[i] = this->aperture_sets[0][i]->weight;
    }
}

void WeightedRobustRecalc::update_weights(const double *weights)
{
    for (auto &apertures : this->aperture_sets)
    {
        for (size_t i = 0; i < apertures.size(); i++)
        {
            apertures[i]->set_weight(weights[i]);
        }
    }
}

void WeightedRobustRecalc::update_weights(const double *weights, std::vector<RecalcAperture *> apertures)
{
    for (size_t i = 0; i < apertures.size(); i++)
    {
        apertures[i]->set_weight(weights[i]);
    }
}

void WeightedRobustRecalc::calculate_structure_dose(Structure *structure, std::vector<RecalcAperture *> apertures)
{
    for (auto &vox : structure->masked_dose)
    {
        vox.second = 0.0;
        size_t vox_num = vox.first;
        for (auto &aperture : apertures)
        {
            vox.second += aperture->dose[vox_num] * aperture->weight;
        }
    }

    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) { return lhs.second < rhs.second; });
}

void WeightedRobustRecalc::update_doses()
{
    for (size_t ap_index = 0; ap_index < this->aperture_sets.size(); ap_index++)
    {
        for (auto &structure : this->structure_sets[ap_index])
        {
            this->calculate_structure_dose(structure, this->aperture_sets[ap_index]);
        }
    }
}

void WeightedRobustRecalc::update_doses(const double *weights)
{
    // Every aperture set has an associated structure set where doses are calculated.
    this->update_weights(weights);
    this->update_doses();
}

double WeightedRobustRecalc::calculate_cost_function(const double *weights,
                                                     bool new_weights)
{
    double cost = 0.0;

    if (new_weights)
    {
        this->update_doses(weights);
    }

    for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
    {
        if (this->scenario_weights[sc_i] > 0.0)
        {
            double sc_cost = 0.0;
            for (auto &structure : this->structure_sets[sc_i])
            {
                sc_cost += structure->calculate_cost();
            }

            cost += sc_cost * this->scenario_weights[sc_i];
        }
    }

    this->latest_cost = cost;
    return cost;
}

void WeightedRobustRecalc::calculate_gradient(double *values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        this->update_doses(weights);
    }

    for (size_t wt_i = 0; wt_i < this->num_weights; wt_i++)
    {
        double grad_value = 0.0;
        for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
        {
            if (this->scenario_weights[sc_i] > 0.0)
            {
                double scenario_grad = 0.0;
                for (auto &structure : this->structure_sets[sc_i])
                {
                    scenario_grad += structure->calculate_gradient(this->aperture_sets[sc_i][wt_i]->dose);
                }
                grad_value += scenario_grad * this->scenario_weights[sc_i];
            }
        }
        values[wt_i] = grad_value;
    }
}

void WeightedRobustRecalc::calculate_hessian(double *hess_values, const double *weights, bool new_weights)
{
    if (new_weights)
    {
        this->update_doses(weights);
    }

    int idx = 0;
    for (size_t row = 0; row < this->num_weights; row++)
    {
        for (size_t col = 0; col <= row; col++)
        {
            double hess = 0.0;
            for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
            {
                if (this->scenario_weights[sc_i] > 0.0)
                {
                    double sc_hess = 0.0;
                    for (auto &structure : this->structure_sets[sc_i])
                    {
                        sc_hess += structure->calculate_hessian(this->aperture_sets[sc_i][row]->dose,
                                                                this->aperture_sets[sc_i][col]->dose);
                    }

                    hess += this->scenario_weights[sc_i] * sc_hess;
                }
            }

            hess_values[idx] = hess;
            idx++;
        }
    }
}

std::vector<double> WeightedRobustRecalc::calculate_total_dose()
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : this->aperture_sets[0])
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

std::vector<double> WeightedRobustRecalc::calculate_scenario_dose(size_t sc_i)
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : this->aperture_sets[sc_i])
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

std::vector<double> WeightedRobustRecalc::calculate_modality_dose(std::vector<RecalcAperture *> apertures)
{
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : apertures)
    {
        for (size_t i = 0; i < total_dose.size(); i++)
        {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

int WeightedRobustRecalc::get_num_apertures()
{
    return this->num_weights;
}

void WeightedRobustRecalc::update_data()
{
    json iter_data;
    iter_data["active_apertures"] = this->aperture_sets[0].size();
    iter_data["cost_function_value"] = this->latest_cost;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - this->start_time).count();

    iter_data["running_time"] = elapsed;

    json e_cpt_stats;
    for (auto &cpt : this->electron_cpts)
    {
        e_cpt_stats.push_back(cpt->write_statistics());
    }
    iter_data["electron_cpt_stats"] = e_cpt_stats;

    json p_cpt_stats;
    for (auto &cpt : this->photon_cpts)
    {
        p_cpt_stats.push_back(cpt->write_statistics());
    }
    iter_data["photon_cpt_stats"] = p_cpt_stats;

    this->plan_data = iter_data;
}

void WeightedRobustRecalc::output_plan_statistics()
{
    std::cout << this->electron_cpts.size() << " electron control points" << std::endl;
    for (auto &cpt : this->electron_cpts)
    {
        cpt->output_statistics();
    }

    std::cout << "---------------------------------------------------------" << std::endl;

    std::cout << this->photon_cpts.size() << " photon control points" << std::endl;
    for (auto &cpt : this->photon_cpts)
    {
        cpt->output_statistics();
    }
}

Coordinates WeightedRobustRecalc::get_coordinates()
{
    if (this->electron_cpts.size())
    {
        return this->electron_cpts[0]->coords;
    }
    else
    {
        if (this->photon_cpts.size() > 1)
        {
            return this->photon_cpts[1]->coords;
        }
        else
        {
            return this->photon_cpts[0]->coords;
        }
    }
}

void WeightedRobustRecalc::output_cost()
{
    for (auto &structure : this->structure_sets[0])
    {
        structure->output_cost();
    }
}

void WeightedRobustRecalc::export_modality_doses()
{
    std::vector<int> electron_energies;
    std::vector<int> photon_energies;
    std::map<int, std::vector<RecalcAperture *>> electron_apertures;
    std::map<int, std::vector<RecalcAperture *>> photon_apertures;

    for (auto &cpt : this->electron_cpts)
    {
        if (electron_apertures.count(cpt->energy) == 0)
        {
            electron_apertures[cpt->energy] = std::vector<RecalcAperture *>();
            electron_energies.push_back(cpt->energy);
        }

        for (auto &ap : cpt->apertures)
        {
            if (ap->active)
            {
                electron_apertures[cpt->energy].push_back(ap);
            }
        }
    }

    for (auto &cpt : this->photon_cpts)
    {
        if (photon_apertures.count(cpt->energy) == 0)
        {
            photon_apertures[cpt->energy] = std::vector<RecalcAperture *>();
            photon_energies.push_back(cpt->energy);
        }

        for (auto &ap : cpt->apertures)
        {
            if (ap->active)
            {
                photon_apertures[cpt->energy].push_back(ap);
            }
        }
    }

    Coordinates coords = this->get_coordinates();

    for (auto energy : electron_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MeV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(electron_apertures[energy]);
        this->write_dose(dose_filename, coords, modality_dose);
    }

    for (auto energy : photon_energies)
    {
        std::string dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(energy) + "MV.3ddose";
        std::cout << "Exporting modality dose to " << dose_filename << "\n";
        std::vector<double> modality_dose = this->calculate_modality_dose(photon_apertures[energy]);
        this->write_dose(dose_filename, coords, modality_dose);
    }
}

void WeightedRobustRecalc::export_scenario_doses(bool interm)
{
    std::string buffer;
    std::string dose_filename;

    Coordinates coords = this->get_coordinates();

    for (size_t sc_i = 0; sc_i < this->num_scenarios; sc_i++)
    {
        if (sc_i == 0)
        {
            dose_filename = "combined_doses/" + remove_extension(input_filename);
        }
        else
        {
            dose_filename = "combined_doses/" + remove_extension(input_filename) + "_" + std::to_string(sc_i);
        }

        if (interm)
        {
            dose_filename += "_interm.3ddose";
        }
        else
        {
            dose_filename += ".3ddose";
        }

        std::cout << "Exporting scenario dose to " << dose_filename << "\n";

        std::vector<double> total_dose = this->calculate_scenario_dose(sc_i);
        this->write_dose(dose_filename, coords, total_dose);
    }
}

void WeightedRobustRecalc::write_dose(const std::string filename, const Coordinates coords, const std::vector<double> dose)
{
    std::ofstream dose_file(filename);

    dose_file << "  " << coords.num_voxels[0] << "  " << coords.num_voxels[1]
              << "  " << coords.num_voxels[2] << "\n";

    for (int i = 0; i < coords.num_voxels[0] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[0] + coords.topleft[0];
    }
    dose_file << "\n";

    for (int i = 0; i < coords.num_voxels[1] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[1] + coords.topleft[1];
    }
    dose_file << "\n";

    for (int i = 0; i < coords.num_voxels[2] + 1; i++)
    {
        if (i != 0)
            dose_file << " ";
        dose_file << i * coords.voxel_size[2] + coords.topleft[2];
    }
    dose_file << "\n";

    for (size_t i = 0; i < dose.size(); i++)
    {
        if (i != 0)
            dose_file << " ";
        if (dose[i] > 0)
        {
            dose_file << dose[i];
        }
        else
        {
            dose_file << "0";
        }
    }
    dose_file << "\n";

    for (size_t i = 0; i < dose.size(); i++)
    {
        dose_file << "0 ";
    }
    dose_file << "\n";

    dose_file.close();
}

void WeightedRobustRecalc::export_final_dose(bool interm)
{
    this->export_scenario_doses(interm);
    if (this->output_modalities)
        this->export_modality_doses();
}

void WeightedRobustRecalc::export_final_weights(bool interm)
{
    std::string buffer;
    std::string plan_filename;

    if (interm)
    {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + "_interm.weights";
    }
    else
    {
        plan_filename = "combined_doses/" + remove_extension(this->input_filename) + ".weights";
    }

    std::cout << "Exporting final weights to " << plan_filename << "\n";
    std::ofstream plan_file(plan_filename);

    json plan_json;
    json structs_json;
    json electron_cpts_json;
    json photon_cpts_json;

    for (auto &roi : this->structure_sets[0])
    {
        structs_json.push_back(roi->to_json());
    }

    for (auto &cpt : this->electron_cpts)
    {
        electron_cpts_json.push_back(cpt->to_json());
    }

    for (auto &cpt : this->photon_cpts)
    {
        photon_cpts_json.push_back(cpt->to_json());
    }

    plan_json["structures"] = structs_json;
    plan_json["electron_cpts"] = electron_cpts_json;
    plan_json["photon_cpts"] = photon_cpts_json;
    plan_json["cost_function_value"] = this->latest_cost;
    plan_file << plan_json.dump();

    plan_file.close();
}

void WeightedRobustRecalc::export_data()
{
    this->update_data();
    std::string filename = "combined_doses/" + remove_extension(this->input_filename) + ".data";
    std::ofstream myfile(filename);
    myfile << this->plan_data.dump();
    myfile.close();
}
