#include <assert.h>
#include "string_utils.hh"
#include "IpoptOptimiser.hh"
#include "TrajectoryGenerator.hh"
#include "IpIpoptApplication.hpp"

using namespace Ipopt;

IpoptOptimiser::IpoptOptimiser(TrajectoryGenerator *traj) {
    this->traj_obj = traj;
}

void IpoptOptimiser::run() {
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);
    std::string stripped_input = remove_extension(this->traj_obj->input_filename);
    app->Options()->SetNumericValue("tol", 1e-8);
    app->Options()->SetIntegerValue("print_level", 5);
    // app->Options()->SetIntegerValue("max_iter", 10000);
    app->Options()->SetIntegerValue("max_iter", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetStringValue("linear_solver", "ma97");
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("derivative_test", "second-order");
    // app->Options()->SetStringValue("derivative_test_print_all", "no");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 0.000001);
    // app->Options()->SetNumericValue("point_perturbation_radius", 0);
    app->Options()->SetStringValue("output_file", stripped_input + ".out");
    // The following overwrites the default name (ipopt.opt) of the
    // options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << std::endl << "*** Error during initialization!" << std::endl;
        return;
    }

    status = app->OptimizeTNLP(this);
}

IpoptOptimiser::~IpoptOptimiser() {}

// returns the size of the problem
bool IpoptOptimiser::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style) {
    n = traj_obj->num_weights;

    m = 0;

    // constraint jacobian is empty since there are no constraints.
    nnz_jac_g = 0;

    // hessian is symmetric so we only need lower triangular part.
    //nnz_h_lag = n * (n + 1) / 2;
    nnz_h_lag = traj_obj->num_hess_ele;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

bool IpoptOptimiser::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u) {
    for (Index i = 0; i < n; i++) {
            x_l[i] = 0.0;
    }

    for (Index i = 0; i < n; i++) {
        x_u[i] = 2e19;
    }

    return true;
}

bool IpoptOptimiser::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda) {
    /*
    if (traj_obj->starting_weights.size() == n) {
        std::cout << "Starting from old weights!" << std::endl;
        for (Index i = 0; i < n; i++) {
            x[i] = traj_obj->starting_weights[i];
        }
    } else {
    */
        for (Index i = 0; i < n; i++) {
            x[i] = 0.0;
        }

    return true;
}

// returns the value of the objective function
bool IpoptOptimiser::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    obj_value = traj_obj->calculate_cost_function(x, new_x);

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptOptimiser::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
    traj_obj->calculate_gradient(grad_f, x, new_x);
    return true;
}

// return the value of the constraints: g(x)
bool IpoptOptimiser::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
  return true;
}

// return the structure or values of the jacobian
bool IpoptOptimiser::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values) {
  return true;
}

// return the structure or values of the hessian
bool IpoptOptimiser::eval_h(Index n, const Number* x, bool new_x,
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
        assert(idx == num_hess_ele);
    } else {
        traj_obj->calculate_hessian(values, x, new_x);
    }

    return true;
}

void IpoptOptimiser::finalize_solution(SolverReturn status,
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

    traj_obj->export_final_dose(x);
}
