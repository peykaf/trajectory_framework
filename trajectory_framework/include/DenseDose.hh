#ifndef DOSE_H
#define DOSE_H 1
#include <vector>
#include "string_utils.hh"
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/tokenizer.hpp>

class DenseDose
{
  public:
    void from_file(std::string filename, bool norm = true);
    void from_file(std::string filename, size_t total_vox, bool norm = true);
    void from_3ddose(std::string filename, bool norm = true);
    void from_3ddose(std::string filename, size_t total_vox, bool norm = true);
    void from_minidos(std::string filename, bool norm = false);
    void from_bindos(std::string filename, bool norm = false);

    double &operator[](std::size_t idx)
    {
        return grid[idx];
    };

    const double &operator[](std::size_t idx) const
    {
        return grid[idx];
    };

    typename std::vector<double>::iterator begin()
    {
        return grid.begin();
    }

    typename std::vector<double>::iterator end()
    {
        return grid.end();
    }
    std::vector<double> grid;

  private:
    int total_voxels;
};

#endif
