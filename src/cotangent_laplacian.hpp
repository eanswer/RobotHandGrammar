#pragma once

#include <Eigen/Sparse>
#include <igl/edge_lengths.h>
#include <igl/face_areas.h>
#include <igl/dihedral_angles.h>
#include <igl/volume.h>
#include "dihedral_sine.hpp"

// TODO: HW3
// Assignment 3, Part 3.2.
/* Implement your code here. */
// Implement the function to compute the cotangent laplacian matrix L 
// V is the vertex matrix of shape (n, 3), each row is the position of a vertex in the mesh
// F is the element index matrix of shape (m, 4), each row is the vertex indices of a tetrahedron
// L is the output cotangent laplacian matrix of shape (n, n), and it's a sparse matrix.
void cotangent_laplacian(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F, 
	Eigen::SparseMatrix<double>& L) 
{
	L.resize(V.rows(), V.rows());
	L.reserve(17 * V.rows()); // Reserve about as much memory as we expect the Laplacian to use

	// Hints:
	// 1. For each tetrahedron, loop over each of its edge,
	//    consider which part of the L matrix this edge in this tetrahedron contributes to
	// 2. compute the cos and sin of the dihedral angle by the law of diehedral angles http://mathworld.wolfram.com/Tetrahedron.html
	//	  specifically, compute the sin and cos of dihedral angles from the edge lengths, face areas and tet volume
	// 3. build the triplets <row, col, value> in IJV

	// Useful Array for iterating over the edges of the tetrahedrons
	// The vertex indices of the e-th edge of the i-th tetrahedron are
	// F(i, edges(e, 0)) and F(i, edges(e, 1))
	Eigen::MatrixXd edges;
	edges.resize(6, 2);
	edges <<
		1, 2,
		2, 0,
		0, 1,
		3, 0,
		3, 1,
		3, 2;

	// Compute Properties of the Tetrahedrons that you will need to compute the Laplacian
	// Face properties are indexed like the vertices opposite them
	// Edge properties are indexed in the same order as the edges array above
	int m = F.rows(); // Number of facets (tetrahedrons)
	Eigen::Matrix<double, Eigen::Dynamic, 6> l; // Edge lengths
	igl::edge_lengths(V, F, l);	
	Eigen::Matrix<double, Eigen::Dynamic, 4> s; // Face areas
	igl::face_areas(l, s);
	Eigen::Matrix<double, Eigen::Dynamic, 6> cos_theta, theta; // Dihedral Angles of each edge and its cosine
	igl::dihedral_angles_intrinsic(l, s, theta, cos_theta);
	Eigen::Matrix<double, Eigen::Dynamic, 1> vol; // Volume of each tetrahedron
	igl::volume(l, vol);
	Eigen::Matrix<double, Eigen::Dynamic, 6> sin_theta(m, 6); // Sine of each edge's dihedral angle
	dihedral_sine(vol, s, l, sin_theta);


	// Compute weighted Cotangents
	Eigen::MatrixXd C = (1. / 6.) * l.array() * cos_theta.array() / sin_theta.array();

	// Add Cotangents together in Laplacian Matrix
	std::vector<Eigen::Triplet<double> > IJV;
	IJV.reserve(F.rows() * edges.rows() * 4);
	// Loop over triangles
	for (int i = 0; i < F.rows(); i++)
	{
		// loop over edges of element
		for (int e = 0; e < edges.rows(); e++)
		{
			int source = F(i, edges(e, 0));
			int dest = F(i, edges(e, 1));
			IJV.push_back(Eigen::Triplet<double>(source, dest, C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(dest, source, C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(source, source, -C(i, e)));
			IJV.push_back(Eigen::Triplet<double>(dest, dest, -C(i, e)));
		}
	}
	// Set From Triplets Sums all Triplets with the same indices
	L.setFromTriplets(IJV.begin(), IJV.end());
}