#ifndef FMHardConstraint_H
#define FMHardConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class FMHardConstraint {
    public:
        FMHardConstraint();
        FMHardConstraint(nlohmann::json &hc_json);
        double calculate_cost(std::vector<double> &struct_dose);

        double calculate_gradient(std::vector<double> &struct_dose,
                                                  std::vector<double> &cpt_dose);

        double calculate_hessian(std::vector<double> &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);


        double threshold;
        double weight;
        int type;
        double latest_cost;

    private:
        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
