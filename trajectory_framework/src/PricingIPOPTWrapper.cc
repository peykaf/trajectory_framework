#include <assert.h>
#include "string_utils.hh"
#include "PricingIPOPTWrapper.hh"
#include "PricingGenerator.hh"
#include "IpIpoptApplication.hpp"

using namespace Ipopt;

PricingIPOPTWrapper::PricingIPOPTWrapper(PricingGenerator *pricingObj) {
    this->pricing_obj = pricingObj;
}

void PricingIPOPTWrapper::run() {
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    this->setIpoptSettings(app, this->pricing_obj->input_filename);

    ApplicationReturnStatus status;
    while (this->pricing_obj->add_new_weight()) {
        status = app->OptimizeTNLP(this);
    }

    if (this->pricing_obj->active_weights.size() > 0) {
        this->pricing_obj->export_final_dose();
        this->pricing_obj->export_final_weights();
        this->pricing_obj->export_data();
    }
}

void PricingIPOPTWrapper::setIpoptSettings(SmartPtr<IpoptApplication> &app, std::string input_filename) {
    app->RethrowNonIpoptException(true);
    app->Options()->SetNumericValue("tol", 1e-8);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma97");
    //app->Options()->SetStringValue("print_info_string", "yes");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    //app->Options()->SetStringValue("ma97_automatic_scaling", "none");
    //app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    //app->Options()->SetStringValue("accept_every_trial_step", "yes");
    //app->Options()->SetIntegerValue("max_soc", 0);
    //app->Options()->SetIntegerValue("watchdog_trial_iter_max", 10);
    //app->Options()->SetNumericValue("nlp_scaling_min_value", 0.0001);
    //app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    //app->Options()->SetStringValue("derivative_test", "second-order");
    //app->Options()->SetStringValue("derivative_test_print_all", "yes");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 0.000001);
    // app->Options()->SetNumericValue("point_perturbation_radius", 0);

    std::string stripped_input = remove_extension(input_filename);
    app->Options()->SetStringValue("output_file", stripped_input + ".out");

}


PricingIPOPTWrapper::~PricingIPOPTWrapper() {}

// returns the size of the problem
bool PricingIPOPTWrapper::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style) {
    n = this->pricing_obj->num_weights;

    m = 0;

    // constraint jacobian is empty since there are no constraints.
    nnz_jac_g = 0;

    // hessian is symmetric so we only need lower triangular part.
    //nnz_h_lag = n * (n + 1) / 2;
    nnz_h_lag = this->pricing_obj->num_hess_ele;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool PricingIPOPTWrapper::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u) {
    for (Index i = 0; i < n; i++) {
        x_l[i] = 0.0;
    }

    for (Index i = 0; i < n; i++) {
        x_u[i] = 2e19;
    }

    return true;
}

bool PricingIPOPTWrapper::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda) {
    // initialize to the given starting point
    for (Index i = 0; i < n; i++) {
        x[i] = this->pricing_obj->active_weights[i]->weight;
    }

    return true;
}

// returns the value of the objective function
bool PricingIPOPTWrapper::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    obj_value = this->pricing_obj->calculate_cost_function(x, new_x);

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool PricingIPOPTWrapper::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
    this->pricing_obj->calculate_gradient(grad_f, x, new_x);
    return true;
}

// return the value of the constraints: g(x)
bool PricingIPOPTWrapper::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
  return true;
}

// return the structure or values of the jacobian
bool PricingIPOPTWrapper::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values) {
  return true;
}

// return the structure or values of the hessian
bool PricingIPOPTWrapper::eval_h(Index n, const Number* x, bool new_x,
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
        this->pricing_obj->calculate_hessian(values, x, new_x);
    }

    return true;
}

void PricingIPOPTWrapper::finalize_solution(SolverReturn status,
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
    this->pricing_obj->update_weights(x);

}
