#ifndef RecalcControlPoint_H
#define RecalcControlPoint_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"

#include "RecalcAperture.hh"
#include "Structure.hh"
#include "Coordinates.hh"

class RecalcControlPoint
{
  public:
    enum particle_type
    {
        photon,
        electron
    };

    RecalcControlPoint();
    RecalcControlPoint(nlohmann::json &cpt_json);
    RecalcControlPoint(nlohmann::json &cpt_json, particle_type ptype);

    void output_statistics();

    nlohmann::json to_json();
    nlohmann::json write_statistics();

    int arc_index;

    particle_type particle;
    double energy;
    bool FFF;

    double gantry_angle;
    double couch_angle;
    double collimator_angle;
    std::vector<double> iso;
    double sad;
    double arclength_scaling;

    std::vector<RecalcAperture *> apertures;

    Coordinates coords;
  protected:
    Coordinates read_header(std::string);
    Coordinates read_3ddose_header(std::string);
    Coordinates read_minidos_header(std::string);
    Coordinates read_bindos_header(std::string);
};

#endif
