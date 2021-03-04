#include <locale>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include "Phantom.hh"
#include "string_utils.hh"

Phantom::Phantom(std::string phant_filename) {
    load_from_EGSphant(phant_filename);
}

Phantom::Phantom() {}

Phantom::~Phantom() {}

void Phantom::load_from_EGSphant(std::string phant_filename) {
    float dummy;
    std::vector<std::string> mats;
    std::string buffer, buffer2;
    std::stringstream convertor;

    std::ifstream phantom_file(phant_filename);
    if (!phantom_file) {
        std::cout << "Could not find " << phant_filename << std::endl;
        exit(-1);
    }
    p_filename = phant_filename;
    std::getline(phantom_file, buffer);
    convertor.str(buffer);

    convertor >> num_mats;
    for (int i = 0; i < num_mats; i++) {
        std::getline(phantom_file, buffer);
        material_names.push_back(trim(buffer));
    }
    std::getline(phantom_file, buffer);  // dummy

    std::getline(phantom_file, buffer);  // number of voxels
    sscanf(buffer.c_str(), "%i %i %i", &num_voxels[0], &num_voxels[1], &num_voxels[2]);
    for (int i = 0; i < 3; i++) {
        std::getline(phantom_file, buffer);
        sscanf(buffer.c_str(), "%f %f %*s", &topleft[i], &dummy);
        voxel_size[i] = dummy - topleft[i];
    }

    int big_counter = 0;
    for (int k = 0; k < num_voxels[2]; k++) {
        for (int j = 0; j < num_voxels[1]; j++) {
            std::getline(phantom_file, buffer);  // Material values
            const char *buf_dup = buffer.c_str();
            for (int i = 0; i < num_voxels[0]; i++) {
                material_ids.push_back(buf_dup[i] - '0' - 1);
                big_counter++;
            }
        }
        std::getline(phantom_file, buffer);
    }

    big_counter = 0;
    for (int k = 0; k < num_voxels[2]; k++) {
        for (int j = 0; j < num_voxels[1]; j++) {
            std::getline(phantom_file, buffer);  // Density values
            std::vector<std::string> densityValues = split(buffer, ' ');
            for (int i = 0; i < num_voxels[0]; i++) {
                densities.push_back(std::atof(densityValues[i].c_str()));
                big_counter++;
            }
        }
        std::getline(phantom_file, buffer);
    }

    phantom_file.close();
}

float Phantom::get_x_min() {
    return topleft[0];
}

float Phantom::get_y_min() {
    return topleft[1];
}

float Phantom::get_z_min() {
    return topleft[2];
}

float Phantom::get_x_max() {
    return (topleft[0] + num_voxels[0] * voxel_size[0]);
}

float Phantom::get_y_max() {
    return (topleft[1] + num_voxels[1] * voxel_size[1]);
}

float Phantom::get_z_max() {
    return (topleft[2] + num_voxels[2] * voxel_size[2]);
}

float Phantom::get_x_size() {
    return voxel_size[0];
}

float Phantom::get_y_size() {
    return voxel_size[1];
}

float Phantom::get_z_size() {
    return voxel_size[2];
}

float Phantom::get_num_x() {
    return num_voxels[0];
}

float Phantom::get_num_y() {
    return num_voxels[1];
}

float Phantom::get_num_z() {
    return num_voxels[2];
}

std::vector<std::string> Phantom::get_mat_names() {
    return material_names;
}

int Phantom::get_mat_id(int vox_number) {
    return material_ids[vox_number];
}

float Phantom::get_density(int vox_number) {
    return densities[vox_number];
}
