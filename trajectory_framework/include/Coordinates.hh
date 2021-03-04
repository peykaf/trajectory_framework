// Copyright 2018 Marc-Andre Renaud
#ifndef Coordinates_H
#define Coordinates_H

class Coordinates {
    public:
    Coordinates() {}
    Coordinates(int *num_vox, double *topl, double *vox_size) {
        this->num_voxels[0] = num_vox[0];
        this->num_voxels[1] = num_vox[1];
        this->num_voxels[2] = num_vox[2];

        this->topleft[0] = topl[0];
        this->topleft[1] = topl[1];
        this->topleft[2] = topl[2];

        this->voxel_size[0] = vox_size[0];
        this->voxel_size[1] = vox_size[1];
        this->voxel_size[2] = vox_size[2];
    }

    int num_voxels[3];
    double topleft[3];
    double voxel_size[3];
};

#endif