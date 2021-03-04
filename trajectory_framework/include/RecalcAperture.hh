#ifndef RECALC_APERTURE_H
#define RECALC_APERTURE_H 1
#include <vector>
#include "json.hh"
#include "WeightClass.hh"

class RecalcApertureRow
{
  public:
    int row_index;
    double l_bound;
    double r_bound;
};

class RecalcAperture : public WeightClass
{
  public:
    RecalcAperture(nlohmann::json ap_json);

    std::vector<RecalcApertureRow> rows;
    nlohmann::json to_json();
};
#endif
