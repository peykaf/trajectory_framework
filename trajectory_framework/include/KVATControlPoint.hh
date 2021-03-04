#ifndef KVATControlPoint_H
#define KVATControlPoint_H 1
#include <vector>
#include <string>
#include <fstream>
#include "json.hh"

#include "Aperture.hh"
#include "SparseDose.hh"
#include "WeightClass.hh"
#include "KVATBeamlet.hh"
#include "Structure.hh"

class KVATControlPoint
{
  public:
    KVATControlPoint();
    KVATControlPoint(nlohmann::json &cpt_json);
    void normalize_beamlets();

    void load_beamlets(std::vector<std::string> beamlet_filenames);
    void load_beamlets();
    void free_beamlets();
    double calculate_beamlet_price(size_t row, size_t beamlet_id, std::vector<double> &lag_mults);
    void price_beamlets(std::vector<double> &lag_mults, std::vector<Structure *> structures);
    double directional_price_norm(size_t row, size_t beamlet_id, std::vector<double> &lag_mults, std::vector<Structure *> structures);
    KVATBeamlet *add_best_beamlet();

    nlohmann::json write_statistics();
    void output_statistics();
    nlohmann::json to_json();

    void write_weight(std::ofstream &ofile, double max_global_weight){};

    size_t num_beamlet_rows;
    size_t num_beamlet_columns;

    int arc_index;

    bool included;
    bool beamlets_loaded;
    double latest_price;

    double gantry_angle;
    double couch_angle;
    double collimator_angle;
    std::vector<double> iso;
    double sad;

    double max_fluence_rate;
    double arclength_scaling;
    double gantry_speed;

    double energy;

    double global_weight;
    double price_norm;

    int best_beamlet;

    // Phantom dimensions
    size_t total_voxels;
    int num_voxels[3];
    double topleft[3];
    double voxel_size[3];

    std::vector<SparseDose> beamlets;
    std::vector<std::string> beamlet_filenames;
    std::vector<KVATBeamlet *> active_beamlets;

    protected:
        void read_header(std::string);
        void read_3ddose_header(std::string);
        void read_minidos_header(std::string);
        void read_bindos_header(std::string);
};

#endif
