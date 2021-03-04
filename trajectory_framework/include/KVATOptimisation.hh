// Copyright 2017 Marc-André Renaud
#ifndef KVATOptimisation_H
#define KVATOptimisation_H 1

#include <vector>
#include "json.hh"

#include "PricingGenerator.hh"
#include "Phantom.hh"
#include "Structure.hh"
#include "KVATControlPoint.hh"

class KVATControlPoint;

class KVATOptimisation : public PricingGenerator {
    public:
        KVATOptimisation();
        KVATOptimisation(nlohmann::json &json_input);

        void from_json(nlohmann::json &json_input);

        bool add_new_weight();

        void prune_cpts();

        void export_final_dose(bool interm=false);
        void export_final_weights(bool interm=false);

    private:
        void price_unused_cpts();
        void output_plan_statistics();

        void add_beamlet(int);

        double avg_aperture_size();
        void update_data();

        bool check_convergence();

        int num_control_points;

        std::vector<KVATControlPoint *> cpts;

        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
