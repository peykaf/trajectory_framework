// Copyright 2017 Marc-André Renaud
#ifndef MixedArcOptimisation_H
#define MixedArcOptimisation_H 1

#include "json.hh"

#include "MixedOptimisation.hh"

class MixedArcOptimisation : public MixedOptimisation {
    public:
        MixedArcOptimisation() {};
        MixedArcOptimisation(nlohmann::json &json_input);
        virtual void export_progress();
        bool add_new_weight();

    protected:
        size_t count_electron_apertures();
        virtual void price_unused_photon_cpts();
};

#endif
