#ifndef WeightClass_H
#define WeightClass_H 1

#include <vector>
#include <fstream>
#include "json.hh"

class WeightClass {
  public:
    WeightClass();
      double weight;
      double norm;
  	std::vector<double> dose;

    // Phantom dimensions
    size_t total_voxels;
    int num_voxels[3];
    double topleft[3];
    double voxel_size[3];
    int weight_type;
    bool active;

	void set_weight(double wt) {
		this->weight = wt;
    };

    double get_weight() {
        return this->weight;
    }

    double normalize_dose();
    double average_dose();

    virtual void deactivate() {
        this->dose.clear();
        this->dose.shrink_to_fit();
        this->active = false;
    };

    virtual nlohmann::json to_json() = 0;

};

#endif
