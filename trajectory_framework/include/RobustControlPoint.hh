#ifndef RobustControlPoint_H
#define RobustControlPoint_H 1

#include "json.hh"

#include "ControlPoint.hh"
#include "Aperture.hh"
#include "SparseDose.hh"

using json = nlohmann::json;

class RobustControlPoint : public ControlPoint
{
  public:
    RobustControlPoint() {}
    RobustControlPoint(json &cpt_json);
    RobustControlPoint(json &cpt_json, particle_type ptype);
    virtual void remove_inactive_apertures();

    size_t get_num_scenarios();
    double calculate_col_price(size_t row, size_t col, std::vector<double> &lag_mults, size_t scenario);
    void price_aperture(std::vector<std::vector<double>> &lag_mults, std::vector<double> &kkt_mults, std::vector<std::vector<Structure *>> structure_sets);
    void price_constrained_aperture(std::vector<std::vector<double>> &lag_mults, std::vector<double> &kkt_mults, std::vector<std::vector<Structure *>> structure_sets, ControlPoint *prev_cpt, ControlPoint *next_cpt);
    std::vector<Aperture *> save_robust_aperture();

  protected:
    std::vector<std::vector<std::string>> robust_filenames;
    std::vector<std::vector<SparseDose>> robust_beamlets;
    std::vector<std::vector<Aperture *>> robust_apertures;

    virtual void load_beamlets();

    double calculate_aperture_dose(Aperture *aperture, const std::vector<SparseDose> &beamlet_set);
    double calculate_transmission_aperture_dose(Aperture *aperture, const std::vector<SparseDose> &beamlet_set);
    double normalize_beamlets(std::vector<SparseDose> &beamlets);
    void normalize_beamlets(std::vector<SparseDose> &beamlets, double norm);
};

#endif