#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "json.hh"

#include "string_utils.hh"
#include "FMStructure.hh"
#include "FMHardConstraint.hh"
#include "FMDVConstraint.hh"

using json = nlohmann::json;

FMStructure::FMStructure() {
    this->voxel_volume = 1.0;
}

FMStructure::FMStructure(json &json_struct) {
    this->voxel_volume = 1.0;
    this->name = json_struct["name"];

    this->mask = json_struct["mask"].get<std::vector<int> >();
    this->masked_dose.resize(this->mask.size());
    std::fill(this->masked_dose.begin(), this->masked_dose.end(), 0.0);
    this->total_voxels = this->mask.size();

    if (json_struct.find("hard_constraints") != json_struct.end()) {
        for (auto& hc_json : json_struct["hard_constraints"]) {
            FMHardConstraint hard_constraint(hc_json);
            this->hard_constraints.push_back(hard_constraint);
        }
    }

    if (json_struct.find("dv_constraints") != json_struct.end()) {
        for (auto& dv_json : json_struct["dv_constraints"]) {
            FMDVConstraint dv_constraint(dv_json);
            dv_constraint.set_volume(this->total_voxels);
            this->dv_constraints.push_back(dv_constraint);
        }
    }

    if (json_struct.find("voxel_volume") != json_struct.end()) {
        this->voxel_volume = json_struct["voxel_volume"];
    }

}

void FMStructure::output_constraints() {
    std::cout << this->name << std::endl;
    if (this->hard_constraints.size() > 0) {
        std::cout << "Hard constraints" << std::endl;
        for (size_t i = 0; i < this->hard_constraints.size(); i++) {
            if (this->hard_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->hard_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Weight: " << this->hard_constraints[i].weight << std::endl;
        }
    }

    if (this->dv_constraints.size() > 0) {
        std::cout << "DV constraints" << std::endl;
        for (size_t i = 0; i < this->dv_constraints.size(); i++) {
            if (this->dv_constraints[i].type == UPPER_LIMIT) {
                std:: cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->dv_constraints[i].threshold << "Gy @ "
                << this->dv_constraints[i].percent_volume << "%%" << std::endl;
            std::cout << "Weight: " << this->dv_constraints[i].weight << std::endl;
        }
    }
}

json FMStructure::to_json() {
    json struct_json;

    struct_json["name"] = this->name;

    if (this->hard_constraints.size() > 0) {
        json hc_constraints;
        for (auto &hc : this->hard_constraints) {
            json hc_json;
            hc_json["threshold"] = hc.threshold;
            hc_json["weight"] = hc.weight;
            hc_json["constraint_type"] = hc.type;
            hc_constraints.push_back(hc_json);
        }
        struct_json["hard_constraints"] = hc_constraints;
    }

    if (this->dv_constraints.size() > 0) {
        json dv_json_constraints;
        for (auto &dv : this->dv_constraints) {
            json dv_json;
            dv_json["threshold"] = dv.threshold;
            dv_json["weight"] = dv.weight;
            dv_json["constraint_type"] = dv.type;
            dv_json["volume"] = dv.percent_volume;

            dv_json_constraints.push_back(dv_json);
        }
        struct_json["dv_constraints"] = dv_json_constraints;
    }

    return struct_json;

}

void FMStructure::write_constraints(std::ofstream &outfile) {
    outfile << "FMStructure" << std::endl;
    outfile << "Name = " << this->name << std::endl;
    outfile << "Num hard constraints = " << this->hard_constraints.size() << std::endl;
    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        outfile << "constraint type = ";
        if (this->hard_constraints[i].type == UPPER_LIMIT) {
            outfile << "Upper limit";
        } else {
            outfile << "Lower limit";
        }
        outfile << std::endl;
        outfile << "Threshold = " << this->hard_constraints[i].threshold << std::endl;
        outfile << "Weight = " << this->hard_constraints[i].weight << std::endl;
    }
    outfile << "Num dv constraints = " << this->dv_constraints.size() << std::endl;
    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        outfile << "constraint type = ";
        if (this->dv_constraints[i].type == UPPER_LIMIT) {
            outfile << "Upper limit";
        } else {
            outfile << "Lower limit";
        }
        outfile << std::endl;
        outfile << "Threshold = " << this->dv_constraints[i].threshold << std::endl;
        outfile << "Percent volume = " << this->dv_constraints[i].percent_volume << std::endl;
        outfile << "Weight = " << this->dv_constraints[i].weight << std::endl;
    }
}

void FMStructure::output_cost() {
    std::cout << this->name << std::endl;
    // Dose inside structures is sorted.
    if (this->hard_constraints.size() > 0) {
        std::cout << "Hard constraints" << std::endl;
        for (size_t i = 0; i < this->hard_constraints.size(); i++) {
            if (this->hard_constraints[i].type == UPPER_LIMIT) {
                std::cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->hard_constraints[i].threshold << "Gy" << std::endl;
            std::cout << "Cost: " << this->hard_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }

    if (this->dv_constraints.size() > 0) {
        std::cout << "DV constraints" << std::endl;
        for (size_t i = 0; i < this->dv_constraints.size(); i++) {
            if (this->dv_constraints[i].type == UPPER_LIMIT) {
                std:: cout << "Upper Limit: ";
            } else {
                std::cout << "Lower Limit: ";
            }
            std::cout << this->dv_constraints[i].threshold << "Gy @ "
                << this->dv_constraints[i].percent_volume << "%" << std::endl;
            std::cout << "Cost: " << this->dv_constraints[i].latest_cost << std::endl;
        }
        std::cout << std::endl;
    }
}

double FMStructure::calculate_cost() {
    double struct_cost = 0.0;

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        struct_cost += this->hard_constraints[i].calculate_cost(this->masked_dose);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        struct_cost += this->dv_constraints[i].calculate_cost(this->masked_dose);
    }

    return struct_cost;
}

double FMStructure::calculate_gradient(std::vector<double> &dose) {
    double gradient = 0.0;

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
        gradient += this->hard_constraints[i].calculate_gradient(this->masked_dose, dose);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        gradient += this->dv_constraints[i].calculate_gradient(this->masked_dose, dose);
    }

    return gradient;
}

double FMStructure::calculate_hessian(std::vector<double> &cpt_one, std::vector<double> &cpt_two) {
    double hessian = 0.0;

    for (size_t i = 0; i < this->hard_constraints.size(); i++) {
       hessian += this->hard_constraints[i].calculate_hessian(this->masked_dose,
                                                               cpt_one,
                                                               cpt_two);
    }

    for (size_t i = 0; i < this->dv_constraints.size(); i++) {
        hessian += this->dv_constraints[i].calculate_hessian(this->masked_dose,
                                                             cpt_one,
                                                             cpt_two);
    }

    return hessian;
}


void FMStructure::set_voxel_volume(double vol) {
    this->voxel_volume = vol;
}
