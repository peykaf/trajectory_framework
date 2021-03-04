#ifndef TrajectoryGenerator_H
#define TrajectoryGenerator_H 1

#include <vector>
#include <string>

class TrajectoryGenerator {
  public:
    virtual double calculate_cost_function(const double *weights, bool new_weights = true) {return 0.0;};
    virtual void calculate_gradient(double *grad_values, const double *weights, bool new_weights = true) {};
    virtual void calculate_hessian(double *hess_values, const double *weights, bool new_weights = true) {};
    virtual void calculate_jacobian(double *values, const double *weights, bool new_weights = true) {};
    virtual void calculate_constraints(double *, const double *, bool) {};
    virtual void export_final_dose(const double *final_weights) {};

    size_t num_weights;
    size_t num_hess_ele;

    std::string input_filename;
    std::vector<double> starting_weights;
};

#endif
