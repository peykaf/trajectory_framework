#include <vector>
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "SparseDose.hh"

using boost::numeric::ublas::compressed_vector;
using boost::numeric::ublas::unbounded_array;

SparseDose::SparseDose() {
    this->loaded = 0;
}

SparseDose::SparseDose(std::string filename) : SparseDose() {
    this->from_file(filename);
}

void SparseDose::from_3ddose(std::string filename, bool norm) {
    int num_voxels[3];
    std::string buffer;
    // Save the input filename to know how to name the final dose.
    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer);  // number of voxels
    std::stringstream convertor(buffer);
    convertor >> num_voxels[0] >> num_voxels[1] >> num_voxels[2];
    size_t total_vox = num_voxels[0] * num_voxels[1] * num_voxels[2];

    this->grid = compressed_vector<float, 0, unbounded_array<int32_t> >(total_vox);

    std::getline(beamlet_file, buffer);  // x voxel coordinates
    std::getline(beamlet_file, buffer);  // y voxel coordinates
    std::getline(beamlet_file, buffer);  // z voxel coordinates
    std::getline(beamlet_file, buffer);  // Dose values

    char *buf_dup = strdup(buffer.c_str());
    grid.resize(total_vox);

    this->max_dose = 0.0;

    char *dumstr = strtok(buf_dup, " ");
    float value;
    size_t counter = 0;
    while (dumstr != NULL) {
        value = atof(dumstr);
        if (value > 0) {
            grid.push_back(counter, value);
            if (value > this->max_dose) {
                this->max_dose = value;
            }
        }
        counter++;
        dumstr = strtok(NULL, " ");
    }

    free(buf_dup);
    this->loaded = 1;
    beamlet_file.close();
}

float SparseDose::normalize_dmax() {
    float max_dose = *std::max_element(this->grid.begin(), this->grid.end());
    this->grid /= max_dose;
    return max_dose;

}

void SparseDose::normalize_dose(float norm) {
    this->grid *= norm;
}

void SparseDose::from_minidos(std::string filename, bool norm){
    int num_nonzero;

    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);

    if (!beamlet_file.is_open()) {
        std::cout << "Could not find file: " << filename << std::endl;
        throw;
    }
    beamlet_file.read((char*)num_voxels, 3 * sizeof(int));
    beamlet_file.read((char*)voxel_size, 3 * sizeof(float));
    beamlet_file.read((char*)topleft, 3 * sizeof(float));
    beamlet_file.read((char*)&num_nonzero, sizeof(int));

    this->max_dose = 0.0;
    size_t total_vox = num_voxels[0] * num_voxels[1] * num_voxels[2];
    this->grid = compressed_vector<float, 0, unbounded_array<int32_t> >(total_vox);
    for (int i = 0; i < num_nonzero; i++) {
        std::pair<int, float> dose;
        beamlet_file.read((char*)&dose.first, sizeof(int));
        beamlet_file.read((char*)&dose.second, sizeof(float));
        //grid(dose.first) = dose.second;
        grid.push_back(dose.first, dose.second);
        if (dose.second > this->max_dose) this->max_dose = dose.second;
    }

    this->loaded = 1;
    beamlet_file.close();
}

void SparseDose::from_bindos(std::string filename){
    int num_nonzero;

    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);

    if (!beamlet_file.is_open()) {
        std::cout << "Could not find file: " << filename << std::endl;
        throw;
    }
    beamlet_file.read((char*)num_voxels, 3 * sizeof(int));

    float *x_voxels = new float[num_voxels[0] + 1];
    float *y_voxels = new float[num_voxels[1] + 1];
    float *z_voxels = new float[num_voxels[2] + 1];

    beamlet_file.read((char*)x_voxels, (num_voxels[0] + 1) * sizeof(float));
    beamlet_file.read((char*)y_voxels, (num_voxels[1] + 1) * sizeof(float));
    beamlet_file.read((char*)z_voxels, (num_voxels[2] + 1) * sizeof(float));

    beamlet_file.read((char*)&num_nonzero, sizeof(int));

    this->max_dose = 0.0;
    if (num_nonzero > 0) {
        int *voxels = new int[num_nonzero];
        float *doses = new float[num_nonzero];

        beamlet_file.read((char*)voxels, num_nonzero * sizeof(int));
        beamlet_file.read((char*)doses, num_nonzero * sizeof(float));


        size_t total_vox = num_voxels[0] * num_voxels[1] * num_voxels[2];
        this->grid = compressed_vector<float, 0, unbounded_array<int32_t> >(total_vox);

        for (size_t i = 0; i < static_cast<size_t>(num_nonzero); i++) {
            grid.push_back(voxels[i], doses[i]);
            if (doses[i] > this->max_dose) this->max_dose = doses[i];
        }

        this->loaded = 1;

        delete[] voxels;
        delete[] doses;
    }

    voxel_size[0] = x_voxels[1] - x_voxels[0];
    voxel_size[1] = y_voxels[1] - y_voxels[0];
    voxel_size[2] = z_voxels[1] - z_voxels[0];

    topleft[0] = x_voxels[0];
    topleft[1] = y_voxels[0];
    topleft[2] = z_voxels[0];

    beamlet_file.close();

    delete[] x_voxels;
    delete[] y_voxels;
    delete[] z_voxels;

}

void SparseDose::from_file(std::string filename, bool norm) {
  std::string extension;
  std::vector<std::string> buffer;

    extension = find_extension(filename);
    if (extension == "3ddose") {
            this->from_3ddose(filename, norm);
    } else if (extension == "minidos") {
            this->from_minidos(filename, norm);
    } else if (extension == "bindos") {
            this->from_bindos(filename);
    } else {
            throw "Dose file extension not recognised";
    }

    if (norm) {
        // Ideally I really don't want to be doing this...
        float norm = 1.0 / this->max_dose;
        this->grid *= norm;
    }
}

void SparseDose::from_file(std::string filename, size_t total_vox, bool norm) {
    this->from_file(filename, norm);
}

