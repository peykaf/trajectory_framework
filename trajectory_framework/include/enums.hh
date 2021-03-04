#ifndef ENUM_H
#define ENUM_H

enum class particle_type
{
    photon,
    electron
};

enum class price_norm_types
{
    dmax,
    directional,
    none
};

enum class mixing_scheme_choices
{
    best_each_particle,
    alternating_particle,
    best_directional,
    best_per_modality,
    all_apertures,
    random
};


#endif