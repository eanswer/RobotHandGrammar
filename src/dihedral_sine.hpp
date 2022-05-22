#pragma once

#include <Eigen/Core>

void dihedral_sine(
	const Eigen::VectorXd& vol,
	const Eigen::Matrix<double, Eigen::Dynamic, 4> & s,
	const Eigen::Matrix<double, Eigen::Dynamic, 6> & l,
	Eigen::Matrix<double, Eigen::Dynamic, 6>& sin_theta)
{
	// Use the law of sines to compute the sine of each dihedral angle
	// http://mathworld.wolfram.com/Tetrahedron.html
	sin_theta.col(0) = vol.array() / ((2. / (3. * l.col(0).array())).array() * s.col(1).array() * s.col(2).array());
	sin_theta.col(1) = vol.array() / ((2. / (3. * l.col(1).array())).array() * s.col(2).array() * s.col(0).array());
	sin_theta.col(2) = vol.array() / ((2. / (3. * l.col(2).array())).array() * s.col(0).array() * s.col(1).array());
	sin_theta.col(3) = vol.array() / ((2. / (3. * l.col(3).array())).array() * s.col(3).array() * s.col(0).array());
	sin_theta.col(4) = vol.array() / ((2. / (3. * l.col(4).array())).array() * s.col(3).array() * s.col(1).array());
	sin_theta.col(5) = vol.array() / ((2. / (3. * l.col(5).array())).array() * s.col(3).array() * s.col(2).array());
}