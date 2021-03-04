#ifndef SparseDose_H
#define SparseDose_H 1
#include <boost/numeric/ublas/vector_sparse.hpp>
#include "string_utils.hh"

typedef boost::numeric::ublas::compressed_vector<float, 0, boost::numeric::ublas::unbounded_array<int32_t> > sparse_vector;

class SparseDose {
    public:
        SparseDose();
        SparseDose(std::string filename);
        void from_file(std::string filename, bool norm = false);
        void from_file(std::string filename, size_t total_vox, bool norm = false);
        void from_minidos(std::string filename, bool norm);
        void from_3ddose(std::string filename, bool norm);
        void from_bindos(std::string filename);
        void normalize_dose(float norm);
        float normalize_dmax();

        typename sparse_vector::iterator begin() {
            return grid.begin();
        }

        typename sparse_vector::iterator end() {
            return grid.end();
        }

        sparse_vector grid;

        float max_dose;
        int num_voxels[3];
        float voxel_size[3];
        float topleft[3];
        int loaded;
};

#endif
