#pragma once

#include <igl/min_quad_with_fixed.h>
#include <igl/invert_diag.h>
#include <Eigen/Sparse>
#include <igl/active_set.h>
#include <igl/normalize_row_sums.h>
#include <igl/massmatrix.h>

#include "cotangent_laplacian.hpp"
#include "quadratic_optimizer.hpp"

// Complete this function by constructing the Q matrix which is the quadratic coefficients for the laplacian energy.
// V is the vertex matrix of shape (n, 3), each row is the position of a vertex in the mesh
// F is the element index matrix of shape (m, 4), each row is the vertex indices of a tetrahedron
// b and bc is the boundary condition comes from your handles, do not touch them
// W is the output weights of shape (n, c)
void bounded_biharmonic_weights(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F,
	const Eigen::VectorXi& b,		
	const Eigen::MatrixXd& bc,
	Eigen::MatrixXd& W
) {

	// Compute the Cotangent Laplacian and construct the quadratic term of the optimization
	Eigen::SparseMatrix<double> L, M, Mi;

	// call the cotangent_laplacian function you wrote in Part 3.2 to get matrix L
	cotangent_laplacian(V, F, L);

	// we help you get the mass matrix M and its inverse Mi here
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Mi);
	
	Eigen::SparseMatrix<double> Q; // you only need to construct matrix Q to complete the pipeline
	Q = (L * Mi * L).eval();
	
	// Weights must be between 0 and 1. Make arrays of these
	// to act as constraint arrays for the solver
	int n = V.rows();
	Eigen::VectorXd ux = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd lx = Eigen::VectorXd::Zero(n);

	// call our solver to solve the quadratic programming
	quadratic_optimization(Q, b, bc, ux, lx, W, 100);
	
	// std::cerr << W << std::endl;

	igl::normalize_row_sums(W, W);
}