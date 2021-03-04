#ifndef HardConstraint_H
#define HardConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class HardConstraint {
    public:
        HardConstraint();
        HardConstraint(nlohmann::json &hc_json);
        double calculate_cost(std::vector<std::pair<int, double> > &struct_dose);
        double calculate_unsorted_cost(std::vector<std::pair<int, double> > &struct_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  std::vector<double> &cpt_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  SparseDose &cpt_dose);

        double calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  std::vector<double> &cpt_dose);

        double calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                                  SparseDose &cpt_dose);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two);

        double calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        double calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two);

        void calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults);

        void calculate_unsorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults);


        nlohmann::json to_json();

        double threshold;
        double weight;
        int type;
        double latest_cost;

    private:
        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
