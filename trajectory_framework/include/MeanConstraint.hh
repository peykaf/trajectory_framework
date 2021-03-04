#ifndef MeanConstraint_H
#define MeanConstraint_H 1
#include <vector>
#include <utility>
#include "json.hh"
#include "SparseDose.hh"

class MeanConstraint {
    public:
        MeanConstraint();
        MeanConstraint(nlohmann::json &mean_json);
        void set_volume(double vol, int total_voxels);
        void set_volume(int total_voxels);

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


        double sum_struct_aperture_dose(std::vector<std::pair<int,double> > &struct_dose,
                                        SparseDose &cpt_dose);

        double sum_struct_aperture_dose(std::vector<std::pair<int,double> > &struct_dose,
                                        std::vector<double> &cpt_dose);

        double calculate_mean_dose(std::vector<std::pair<int, double> > &struct_dose);

        nlohmann::json to_json();

        double threshold;
        double weight;
        int type;
        double percent_volume;
        int voxels_considered;
        double fraction;
        
        double latest_cost;


    private:
        int start;
        int real_start;
        int end;
        const int LOWER_LIMIT = 0;
        const int UPPER_LIMIT = 1;
};

#endif
