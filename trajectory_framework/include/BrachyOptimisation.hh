// Copyright 2015 Marc-André Renaud
#ifndef BrachyOptimisation_h
#define BrachyOptimisation_h 1

#include <fstream>
#include <string>
#include <vector>
#include "json.hh"
#include "PricingGenerator.hh"
#include "Phantom.hh"
#include "Structure.hh"
#include "DwellPosition.hh"

class BrachyOptimisation : public PricingGenerator {
    public:
        BrachyOptimisation();
        BrachyOptimisation(nlohmann::json &input_json);

        void from_json(nlohmann::json &input_json);

        bool add_new_weight();

        void export_final_dose(bool interm=false);
        void export_final_weights(bool interm=false);

    private:
        void price_unused_dwells();
        void update_data();

        int num_control_points;

        nlohmann::json plan_data;
        std::vector<nlohmann::json> cost_data;
        std::vector<nlohmann::json> dvh_data;

        std::vector<DwellPosition *> dwells;
};

#endif
