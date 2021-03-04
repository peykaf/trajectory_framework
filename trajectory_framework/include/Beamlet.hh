#ifndef BEAMLET_H
#define BEAMLET_H 1
#include <vector>

class Beamlet {
    public:
        double x_size;
        double y_size;
        std::vector<double> dose;
};

#endif
