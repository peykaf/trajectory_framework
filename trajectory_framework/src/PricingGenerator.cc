#include <string>
#include <fstream>
#include "json.hh"

#include "PricingGenerator.hh"

PricingGenerator::PricingGenerator() {
    this->max_apertures = 0;
    this->latest_cost = 0.0;
    this->previous_cost = 0.0;
    this->iteration_number = 0;
}

void PricingGenerator::calculate_structure_dose(Structure *structure) {
    for (auto &vox : structure->masked_dose) {
        vox.second = 0.0;
        size_t vox_num = vox.first;
        for (auto &aperture : this->active_weights) {
            vox.second += aperture->dose[vox_num] * aperture->weight;
        }
    }

    std::sort(structure->masked_dose.begin(),
              structure->masked_dose.end(),
              [](const std::pair<int, double>& lhs,const std::pair<int, double>& rhs){ return lhs.second < rhs.second; });
}

void PricingGenerator::output_cost() {
    for (auto &structure : this->structures) {
        structure->output_cost();
    }
}

void PricingGenerator::refresh_doses() {
    for (auto &structure : this->structures) {
        this->calculate_structure_dose(structure);
    }
}

void PricingGenerator::update_weights(const double *weights) {
    for (size_t i = 0; i < this->active_weights.size(); i++) {
        this->active_weights[i]->set_weight(weights[i]);
    }
}

double PricingGenerator::calculate_cost_function(const double *weights,
                                                bool new_weights) {
    double cost = 0.0;

    if (new_weights) {
        this->update_weights(weights);
        this->refresh_doses();
    }

    for (auto &structure: this->structures) {
        cost += structure->calculate_cost();
    }

    this->latest_cost = cost;
    return cost;
}

void PricingGenerator::calculate_gradient(double *grad_values,
                                         const double *weights,
                                         bool new_weights) {
    if (new_weights) {
        this->update_weights(weights);
        this->refresh_doses();
    }

    size_t idx = 0;
    for (auto &aperture : this->active_weights) {
        double grad_value = 0.0;
        for (auto &structure : this->structures) {
            grad_value += structure->calculate_gradient(aperture->dose);
        }

        grad_values[idx] = grad_value;
        idx += 1;
    }
}

void PricingGenerator::calculate_hessian(double* hess_values,
                                        const double *weights,
                                        bool new_weights) {
    if (new_weights) {
        this->update_weights(weights);
        this->refresh_doses();
    }

    int idx = 0;
    for (size_t row = 0; row < this->active_weights.size(); row++) {
        for (size_t col = 0; col <= row; col++) {
            double hess = 0.0;
            for (auto &structure : this->structures) {
                hess += structure->calculate_hessian(this->active_weights[row]->dose,
                                                     this->active_weights[col]->dose);
            }

            hess_values[idx] = hess;
            idx++;
        }
    }
}

void PricingGenerator::calculate_lag_mults() {
    std::fill(this->lag_mults.begin(), this->lag_mults.end(), 0.0);
    for (auto &structure : this->structures) {
        structure->calculate_lag_mults(this->lag_mults);
    }
}

std::vector<double> PricingGenerator::calculate_total_dose() {
    std::vector<double> total_dose(this->total_voxels, 0.0);
    for (auto &aperture : this->active_weights) {
        for(size_t i = 0; i < total_dose.size(); i++) {
            total_dose[i] += aperture->dose[i] * aperture->weight;
        }
    }

    return total_dose;
}

void PricingGenerator::export_data() {
    std::string filename = "combined_doses/" + remove_extension(this->input_filename) + ".data";
    std::ofstream myfile(filename);
    myfile << this->plan_data.dump();
    myfile.close();

    std::string filename3 = "combined_doses/" + remove_extension(this->input_filename) + ".cost";
    std::ofstream myfile3(filename3);
    myfile3 << nlohmann::json(this->cost_data).dump();
    myfile3.close();

    std::string filename4 = "combined_doses/" + remove_extension(this->input_filename) + ".dvh";
    std::ofstream myfile4(filename4);
    myfile4 << nlohmann::json(this->dvh_data).dump();
    myfile4.close();

    this->export_progress();
}

void PricingGenerator::export_progress() {
    std::string filename2 = "combined_doses/" + remove_extension(this->input_filename) + ".progress";
    std::ofstream myfile2(filename2);
    myfile2 << this->active_weights.size() << "/" << this->max_apertures << std::endl;
    myfile2.close();

}