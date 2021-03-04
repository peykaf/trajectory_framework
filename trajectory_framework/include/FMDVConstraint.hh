#ifndef FMDVConstraint_H
#define FMDVConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class FMDVConstraint {
    public:
        FMDVConstraint();
        FMDVConstraint(nlohmann::json &dv_json);
        double calculate_cost(std::vector<double> &struct_dose);

        double calculate_gradient(std::vector<double> &struct_dose,
                                  std::vector<double> &cpt_dose);

        double calculate_hessian(std::vector<double> &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        void find_violator_doses(std::vector<double> &struct_dose);

        void set_volume(double vol, int total_voxels);
        void set_volume(int total_voxels);

        double threshold;
        double percent_volume;
        double weight;
        double dv_threshold;
        int voxels_considered;
        int type;
        int allowed_violators;

        double latest_cost;
        std::vector<std::pair<int, double> > violator_doses;

    private:
        std::size_t start;

        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
