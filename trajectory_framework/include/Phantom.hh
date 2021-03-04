#ifndef Phantom_h
#define Phantom_h 1

#include <string>
#include <vector>

class Phantom {
  public:
    Phantom();
    Phantom(std::string phant_filename);
    ~Phantom();

    float get_x_min();
    float get_y_min();
    float get_z_min();

    float get_x_max();
    float get_y_max();
    float get_z_max();

    float get_x_size();
    float get_y_size();
    float get_z_size();

    float get_num_x();
    float get_num_y();
    float get_num_z();

    std::vector<std::string> get_mat_names();
    int get_mat_id(int vox_number);
    float get_density(int vox_number);

    std::string p_filename;
    int num_mats;

    void load_from_EGSphant(std::string);

    int num_voxels[3];
    float voxel_size[3];
    float topleft[3];

  private:
    std::vector<int> material_ids;
    std::vector<float> densities;
    std::vector<std::string> material_names;
};

#endif
