#ifndef PTVConstraint_H
#define PTVConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class PTVConstraint {
    public:
        PTVConstraint();
        PTVConstraint(nlohmann::json &hc_json);
        double calculate_cost(std::vector<std::pair<int, double> > &struct_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  std::vector<double> &cpt_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  SparseDose &cpt_dose);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two);

        void calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults);

        nlohmann::json to_json();

        double threshold;
        double weight;
        double latest_cost;
};

#endif
