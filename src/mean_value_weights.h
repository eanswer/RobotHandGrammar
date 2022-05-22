#pragma once

#include "Common/Common.h"

// implements by following the paper "Mean Value Coordinates for Closed Triangular Meshes"
// Ju et al. SIGGRAPH 2005, https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf

void ComputeMeanValueWeights(
    const Eigen::MatrixXd& V, // vertices to compute weights
    const Eigen::MatrixXd& P, // handle points on the cages           
    const Eigen::MatrixXi& Tri, // triangle elements of the cages.
    Eigen::MatrixXd& W
);