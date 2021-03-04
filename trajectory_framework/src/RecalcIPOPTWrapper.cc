#include <assert.h>
#include "string_utils.hh"
#include "RecalcIPOPTWrapper.hh"
#include "PricingGenerator.hh"
#include "WeightedRobustRecalc.hh"
#include "IpIpoptApplication.hpp"

using namespace Ipopt;

RecalcIPOPTWrapper::RecalcIPOPTWrapper(WeightedRobustRecalc *robustObj) {
    this->robust_obj = robustObj;
}

void RecalcIPOPTWrapper::run() {
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    this->setIpoptSettings(app, this->robust_obj->input_filename);

    ApplicationReturnStatus status;
    status = app->OptimizeTNLP(this);
}

void RecalcIPOPTWrapper::setIpoptSettings(SmartPtr<IpoptApplication> &app, std::string input_filename) {
    app->RethrowNonIpoptException(true);
    app->Options()->SetNumericValue("tol", 1e-8);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma97");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");

    std::string stripped_input = remove_extension(input_filename);
    app->Options()->SetStringValue("output_file", stripped_input + ".out");

}


RecalcIPOPTWrapper::~RecalcIPOPTWrapper() {}

// returns the size of the problem
bool RecalcIPOPTWrapper::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style) {
    n = this->robust_obj->num_weights;

    m = 0;

    // constraint jacobian is empty since there are no constraints.
    nnz_jac_g = 0;

    // hessian is symmetric so we only need lower triangular part.
    nnz_h_lag = n * (n + 1) / 2;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool RecalcIPOPTWrapper::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u) {
    for (Index i = 0; i < n; i++) {
        x_l[i] = 0.0;
    }

    for (Index i = 0; i < n; i++) {
        x_u[i] = 2e19;
    }

    return true;
}

bool RecalcIPOPTWrapper::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda) {
    // initialize to the given starting point
    this->robust_obj->get_starting_point(x);
    return true;
}

// returns the value of the objective function
bool RecalcIPOPTWrapper::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    obj_value = this->robust_obj->calculate_cost_function(x, new_x);

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool RecalcIPOPTWrapper::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
    this->robust_obj->calculate_gradient(grad_f, x, new_x);
    return true;
}

// return the value of the constraints: g(x)
bool RecalcIPOPTWrapper::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
  return true;
}

// return the structure or values of the jacobian
bool RecalcIPOPTWrapper::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values) {
  return true;
}

// return the structure or values of the hessian
bool RecalcIPOPTWrapper::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values) {
    if (values == NULL) {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.

        Index idx = 0;
        for (Index row = 0; row < n; row++) {
            for (Index col = 0; col <= row; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }
    } else {
        this->robust_obj->calculate_hessian(values, x, new_x);
    }

    return true;
}

void RecalcIPOPTWrapper::finalize_solution(SolverReturn status,
      Index n, const Number* x, const Number* z_L, const Number* z_U,
      Index m, const Number* g, const Number* lambda,
      Number obj_value, const IpoptData* ip_data,
      IpoptCalculatedQuantities* ip_cq) {
    std::cout << std::endl << std::endl
    << "Solution of the weight variables, x" << std::endl;

    for (Index i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;

    this->robust_obj->update_weights(x);
    this->robust_obj->export_final_dose();
    this->robust_obj->export_final_weights();
    this->robust_obj->export_data();
}
