#ifndef RobustRecalcControlPoint_H
#define RobustRecalcControlPoint_H 1

#include "json.hh"

#include "RecalcControlPoint.hh"
#include "RecalcAperture.hh"

using json = nlohmann::json;

class RobustRecalcControlPoint : public RecalcControlPoint
{
  public:
    RobustRecalcControlPoint() {}
    RobustRecalcControlPoint(json &cpt_json);
    RobustRecalcControlPoint(json &cpt_json, particle_type ptype);

    size_t get_num_scenarios();
    std::vector<std::vector<RecalcAperture *>> robust_apertures;
};

#endif