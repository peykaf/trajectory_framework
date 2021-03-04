// Copyright 2018 Marc-Andre Renaud
#ifndef RobustArcOptimisation_h
#define RobustArcOptimisation_h
#include <vector>
#include <chrono>
#include "json.hh"

#include "TrajectoryGenerator.hh"
#include "WeightClass.hh"
#include "Structure.hh"
#include "RobustControlPoint.hh"
#include "RobustOptimisation.hh"

class RobustArcOptimisation : public RobustOptimisation  {
    public:
        RobustArcOptimisation() {};
        RobustArcOptimisation(nlohmann::json &json_input);
        void export_progress();
    protected:
        size_t count_electron_apertures();
        void price_unused_photon_cpts();
        virtual bool add_new_aperture();
};

#endif
