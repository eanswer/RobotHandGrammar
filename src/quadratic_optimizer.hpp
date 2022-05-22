#pragma once
#include <igl/min_quad_with_fixed.h>
#include <igl/active_set.h>
#include <igl/parallel_for.h>
#include <Eigen/Sparse>
#include <mutex>

bool quadratic_optimization(
	const Eigen::SparseMatrix<double>& Q, // Quadratic Optimization Coefficients,
	const Eigen::MatrixXi& b, // Indices of fixed boundary conditions
	const Eigen::MatrixXd& bc, // Values of fixed boundary conditions
	const Eigen::VectorXd& ux, // Weight upper bounds - pass an empty vector for no bounds
	const Eigen::VectorXd& lx, // Weight lower bounds - pass an empty vector for no bounds
	Eigen::MatrixXd& W, // Output: Optimized Function
	int num_iterations = 8 // Optional Input: Max iterations for solver
	)
{
	int n = Q.rows();
	int m = bc.cols();

	// Our Optimization has none of these terms or constraints:
	Eigen::VectorXd c = Eigen::VectorXd::Zero(n);
	Eigen::SparseMatrix<double> A(0, n), Aeq(0, n), Aieq(0, n);
	Eigen::VectorXd Beq(0, 1), Bieq(0, 1);

	// Set Solver Options: Number of iterations to run
	auto solver_params = igl::active_set_params();
	solver_params.max_iter = num_iterations;

	// Solve without upper/lower bounds constraints to get an initial guess
	igl::min_quad_with_fixed_data<double> mqwf;
	igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, mqwf);
	igl::min_quad_with_fixed_solve(mqwf, c, bc, Beq, W);

	// This is essentially the "first iteration"
	solver_params.max_iter--;

	std::cerr << "finished preconditioning and initial guess" << std::endl;

	// Compute Weights for each handle separately
	bool error = false;
	std::mutex critical;
	const auto& optimize_weight = [&](const int i)
	{
		// Quicker exit for parallel_for
		if (error)
		{
			return;
		}
        {
            std::lock_guard<std::mutex> lock(critical);
            std::cout << "BBW: Computing weight for handle " << i + 1 << " out of " << m <<
                    "." << std::endl;
        }
		Eigen::VectorXd bci = bc.col(i);
		Eigen::VectorXd Wi;
		// use initial guess
		Wi = W.col(i);
		igl::SolverStatus ret = igl::active_set(
			Q, c, b, bci, Aeq, Beq, Aieq, Bieq, lx, ux, solver_params, Wi);
		switch (ret)
		{
		case igl::SOLVER_STATUS_CONVERGED:
			break;
		case igl::SOLVER_STATUS_MAX_ITER:
			std::cerr << "active_set: max iter without convergence." << std::endl;
			break;
		case igl::SOLVER_STATUS_ERROR:
		default:
			std::cerr << "active_set error." << std::endl;
			error = true;
		}
		W.col(i) = Wi;
	};
    #ifdef WIN32
        for (int i = 0; i < m; ++i)
            optimize_weight(i);
    #else
        igl::parallel_for(m,optimize_weight, 5);
    #endif

	return !error;
}