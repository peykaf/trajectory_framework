#ifndef DVConstraint_H
#define DVConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class DVConstraint {
    public:
        DVConstraint();
        DVConstraint(nlohmann::json &dv_json);
        double calculate_unsorted_cost(std::vector<std::pair<int, double> > &struct_dose);

        double calculate_sorted_cost(std::vector<std::pair<int, double> > &struct_dose);

        double calculate_cost(std::vector<std::pair<int, double> > &struct_dose,
                                  bool sorted);

        double calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  SparseDose &cpt_dose);

        double calculate_sorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  SparseDose &cpt_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  SparseDose &cpt_dose,
                                  bool sorted);

        double calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two);

        double calculate_sorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 SparseDose &cpt_one,
                                 SparseDose &cpt_two,
                                 bool sorted);

        void calculate_unsorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults);

        void calculate_sorted_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults);

        void calculate_lag_mults(std::vector<std::pair<int, double> > &struct_dose,
                               std::vector<double> &lag_mults,
                               bool sorted);

        void set_volume(double vol, int total_voxels);
        void set_volume(int total_voxels);

        double calculate_unsorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        double calculate_sorted_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two);

        double calculate_hessian(std::vector<std::pair<int, double> > &struct_dose,
                                 std::vector<double> &cpt_one,
                                 std::vector<double> &cpt_two,
                                 bool sorted);

        double calculate_unsorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  std::vector<double> &cpt_dose);

        double calculate_sorted_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  std::vector<double> &cpt_dose);

        double calculate_gradient(std::vector<std::pair<int, double> > &struct_dose,
                                  std::vector<double> &cpt_dose,
                                  bool sorted);

        nlohmann::json to_json();

        double threshold;
        double percent_volume;
        double weight;
        int voxels_considered;
        int type;

        double latest_cost;

    private:
        std::size_t start;

        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
