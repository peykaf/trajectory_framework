// Copyright 2018 Marc-Andre Renaud
#include <iostream>
#include <vector>

#include "json.hh"
#include "SparseDose.hh"
#include "RobustRecalcControlPoint.hh"
#include "RecalcControlPoint.hh"

using json = nlohmann::json;
typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t> > sparse_vector;

RobustRecalcControlPoint::RobustRecalcControlPoint(json &cpt_json) : RecalcControlPoint(cpt_json)
{
    // Robust beamlet filenames are the only RobustRecalcControlPoint-specific feature.
    this->robust_apertures.resize(cpt_json["robust_apertures"].size());

    for (size_t i = 0; i < this->robust_apertures.size(); i++) {
        for (auto &ap : cpt_json["robust_apertures"][i]) {
            this->robust_apertures[i].push_back(new RecalcAperture(ap));
        }
    }
}

RobustRecalcControlPoint::RobustRecalcControlPoint(json &cpt_json, particle_type ptype) : RobustRecalcControlPoint(cpt_json)
{
    this->particle = ptype;
}

size_t RobustRecalcControlPoint::get_num_scenarios()
{
    return this->robust_apertures.size() + 1;
}