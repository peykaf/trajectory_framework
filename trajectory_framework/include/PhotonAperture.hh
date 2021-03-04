#ifndef PhotonAperture_H
#define PhotonAperture_H 1

#include "Aperture.hh"


class PhotonAperture : public Aperture {
public:
	PhotonAperture(size_t num_ap) : Aperture(num_ap) {}
	PhotonAperture() {}
};

#endif
