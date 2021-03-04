#ifndef PricingGenerator_H
#define PricingGenerator_H 1

#include <vector>
#include <string>
#include <list>
#include "json.hh"

#include "TrajectoryGenerator.hh"
#include "Structure.hh"
#include "WeightClass.hh"

class PricingGenerator : public TrajectoryGenerator {
  public:
    PricingGenerator();
    virtual void calculate_structure_dose(Structure *structure);
    virtual void update_weights(const double *);
    void export_data();
    virtual void export_progress();

    virtual bool add_new_weight() = 0;
    virtual void export_final_dose(bool interm=false) = 0;
    virtual void export_final_weights(bool interm=false) = 0;

    double calculate_cost_function(const double *weights,
                                           bool new_weights = true);
    void calculate_gradient(double *grad_values,
                            const double *weights,
                            bool new_weights = true);

    void calculate_hessian(double* hess_values,
                           const double *weights,
                           bool new_weights = true);

    void calculate_jacobian() {}

    virtual std::vector<double> calculate_total_dose();
    void output_cost();

    std::vector<Structure *> get_structures() {return this->structures;}

    std::vector<WeightClass *> active_weights;
    std::list<int> last_apertures;
    std::list<double> last_cost;

  protected:
    int iteration_number;
    double latest_cost;
    double previous_cost;
    size_t total_voxels;
    size_t max_apertures;

    std::vector<double> lag_mults;
    std::vector<Structure *> structures;

    nlohmann::json plan_data;
    std::vector<nlohmann::json> cost_data;
    std::vector<nlohmann::json> dvh_data;

    void calculate_lag_mults();
    void refresh_doses();

    const int LOWER_LIMIT = 0;
    const int UPPER_LIMIT = 1;
};

#endif
