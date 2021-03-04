#ifndef DwellPosition_H
#define DwellPosition_H 1
#include <vector>
#include <string>
#include "json.hh"
#include "SparseDose.hh"
#include "WeightClass.hh"

class DwellPosition : public WeightClass {
  public:
    DwellPosition();
    DwellPosition(nlohmann::json &dwell_json);

    bool included;
    double latest_price;

    int catheter_index;
    int angle_index;
    std::vector<double> position;
    double shield_angle;

    std::string dose_filename;

    void load_beamlets(std::vector<std::string> beamlet_filenames);
    void load_beamlets();
    void free_beamlets();
    void read_dose(std::string);
    void read_3ddose_header(std::string);
    void read_binary_header(std::string);

    void price_dwell(std::vector<double> &lag_mults);

    nlohmann::json to_json();
    nlohmann::json write_statistics();
    void output_statistics();
};

#endif
