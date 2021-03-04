// Copyright 2018 Marc-Andre Renaud
#ifndef WeightedRobustArcOptimisation_h
#define WeightedRobustArcOptimisation_h
#include <vector>
#include <chrono>
#include "json.hh"

#include "TrajectoryGenerator.hh"
#include "WeightClass.hh"
#include "Structure.hh"
#include "RobustControlPoint.hh"
#include "WeightedRobustOptimisation.hh"

class WeightedRobustArcOptimisation : public WeightedRobustOptimisation  {
    public:
        WeightedRobustArcOptimisation() {};
        WeightedRobustArcOptimisation(nlohmann::json &json_input);
        void export_progress();
        bool add_new_weight();
    protected:
        size_t count_electron_apertures();
        void price_unused_photon_cpts();
};

#endif