#pragma once
#include <Eigen/Core>

void ComputeNearestNeighborWeights(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXd& C, 
	Eigen::MatrixXd& W)
{
	W.resize(V.rows(), C.rows());
	for (int i = 0; i < V.rows(); ++i) {
		int nearest = 0;
		for (int j = 0; j < C.rows(); ++j) {
			W(i, j) = 0;
			if ((V.row(i) - C.row(j)).norm() < (V.row(i) - C.row(nearest)).norm()) {
				nearest = j;
			}
		}
		W(i, nearest) = 1;
	}
}