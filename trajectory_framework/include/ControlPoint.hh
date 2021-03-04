#ifndef ControlPoint_H
#define ControlPoint_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"

#include "SparseDose.hh"
#include "WeightClass.hh"
#include "Aperture.hh"
#include "Structure.hh"
#include "enums.hh"

class ControlPoint
{
  public:


    ControlPoint();
    ControlPoint(nlohmann::json &cpt_json);
    ControlPoint(nlohmann::json &cpt_json, particle_type ptype);
    Aperture *save_current_aperture();
    void update_aperture_dose(Aperture *, size_t, size_t, size_t);

    void output_statistics();

    void normalize_beamlets();

    void set_price_norm(price_norm_types);

    virtual void remove_inactive_apertures();

    void price_aperture(std::vector<double> &lag_mults, std::vector<Structure *> structures);

    void price_constrained_aperture(std::vector<double> &lag_mults,
                                    std::vector<Structure *> structures,
                                    ControlPoint *prev_cpt,
                                    ControlPoint *next_cpt);

    bool is_row_deliverable(ControlPoint *, ControlPoint *, size_t, size_t, size_t);
    nlohmann::json to_json();
    nlohmann::json write_statistics();

    size_t num_beamlet_rows;
    size_t num_beamlet_columns;

    int arc_index;

    bool included;
    bool beamlets_loaded;
    bool dummy;
    double latest_price;

    double col_size;
    double row_size;
    particle_type particle;

    double gantry_angle;
    double couch_angle;
    double collimator_angle;
    std::vector<double> iso;
    double sad;

    double max_fluence_rate;
    double max_leaf_speed;
    double max_overtravel;
    double arclength_scaling;
    double gantry_speed;
    double transmission_factor;

    int max_l_column;
    int min_r_column;

    double energy;
    bool FFF;
    bool is_arc;

    double global_weight;
    double price_norm;

    std::vector<std::string> beamlet_filenames;
    Aperture current_aperture;
    std::vector<Aperture *> apertures;

    // Phantom dimensions
    size_t total_voxels;
    int num_voxels[3];
    double topleft[3];
    double voxel_size[3];

  protected:
    void read_header(std::string);
    void read_3ddose_header(std::string);
    void read_minidos_header(std::string);
    void read_bindos_header(std::string);
    void load_beamlets(std::vector<std::string> beamlet_filenames);
    virtual void load_beamlets();

    void free_beamlets();
    void init_aperture();

    double calculate_aperture_dose(Aperture &aperture);
    double calculate_transmission_aperture_dose(Aperture &aperture);
    void calculate_aperture_size();

    double calculate_col_price(size_t row, size_t col, std::vector<double> &lag_mults);
    double directional_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures);
    double dmax_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures);
    double average_price_norm(std::vector<double> &lag_mults, std::vector<Structure *> structures);

    bool do_mlc_transmission;
    std::vector<SparseDose> beamlets;
    price_norm_types price_norm_type;

};

#endif
