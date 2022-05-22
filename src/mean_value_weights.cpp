#include "mean_value_weights.h"

// implements by following the paper "Mean Value Coordinates for Closed Triangular Meshes"
// Ju et al. SIGGRAPH 2005, https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf

void ComputeMeanValueWeights(
    const Eigen::MatrixXd& V, // vertices to compute weights
    const Eigen::MatrixXd& P, // handle points on the cages           
    const Eigen::MatrixXi& Tri, // triangle elements of the cages.
    Eigen::MatrixXd& W
) {
    W.resize(V.rows(), P.rows());
    W.setZero();

    int n = V.rows(), m = P.rows();
    double EPS = 1e-7;

    Eigen::VectorXd d(m);
    Eigen::MatrixXd u(m, 3);
    for (int i = 0;i < n;i++) {
        const Eigen::Vector3d& v = V.row(i);
        
        bool flag = false; // flag = true if boundary cases are triggered

        // d(j) = |P(j) - v|
        // u(j) = (P(j) - v) / d(j)
        for (int j = 0;j < m;j++) {
            const Eigen::Vector3d &p = P.row(j);
            Eigen::Vector3d vec = p - v;
            d(j) = vec.norm();
            if (d(j) < EPS) {
                W.row(i).setZero();
                W(i, j) = 1.;
                flag = true;
                break;
            }
            u.row(j) = vec / d(j);
        }

        if (flag) {
            continue;
        }

        // compute the contribution from each triangle.

        for (int j = 0;j < Tri.rows();j++) {
            flag = false;

            // compute l and theta
            // l(k) = |u(k + 1) - u(k - 1)|
            // theta(k) = 2arcsin(l(k) / 2)
            std::vector<double> theta;
            double h = 0.;
            for (int k = 0;k < 3;k++) {
                // int p1 = Tri(j, (k + 1) % 3), p2 = Tri(j, (k + 2) % 3);
                // std::cerr << "p1 = " << p1 << ", p2 = " << p2 << std::endl;
                const Eigen::Vector3d &u1 = u.row(Tri(j, (k + 1) % 3)), &u2 = u.row(Tri(j, (k + 2) % 3));
                double l = (u1 - u2).norm();
                theta.push_back(2. * asin(l / 2.));
                h += theta[k] / 2.;
            }

            if (pi - h < EPS) { // v(i) lines on tri(j), use 2D barycentric coordinates
                W.row(i).setZero();
                for (int k = 0;k < 3;k++) {
                    W(i, Tri(j, k)) = sin(theta[k]) * d(Tri(j, (k + 2) % 3)) * d(Tri(j, (k + 1) % 3));
                    // std::cerr << sin(theta[k]) * d(Tri(j, (k + 2) % 3) * d(Tri(j, (k + 1) % 3))) << ", sin = " << sin(theta[k]) << ", d(-1) = " << d(Tri(j, (k + 2) % 3)) << ", d(+1) = " << d(Tri(j, (k + 1) % 3)) << std::endl;
                }
                flag = true;
                break;
            }

            // sign = sign(det[u1, u2, u3])
            double sign = -1.;
            const Eigen::Vector3d &u0 = u.row(Tri(j, 0)), &u1 = u.row(Tri(j, 1)), &u2 = u.row(Tri(j, 2));
            if (u0.dot(u1.cross(u2)) > EPS) {
                sign = 1.;
            }

            // c(k) = (2 * sin(h) * sin(h - theta(k))) / (sin(theta(k + 1)) * sin(theta(k - 1))) - 1
            // s(k) = sign * sqrt(1 - c(k)^2)
            std::vector<double> c, s;
            for (int k = 0;k < 3;k++) {
                c.push_back((2. * sin(h) * sin(h - theta[k])) / (sin(theta[(k + 1) % 3]) * sin(theta[(k + 2) % 3])) - 1.);
                s.push_back(sign * sqrt(max(0., 1. - c[k] * c[k])));
                if (fabs(s[k]) < EPS) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                continue;
            }
            
            // w(k) = (theta(k) - c(k + 1) * theta(k - 1) - c(k - 1) * theta(k + 1)) / (d(k) * sin(theta(k + 1)) * s(k - 1))
            for (int k = 0;k < 3;k++) {
                double w = (theta[k] - c[(k + 1) % 3] * theta[(k + 2) % 3] - c[(k + 2) % 3] * theta[(k + 1) % 3])
                            / (d(Tri(j, k)) * sin(theta[(k + 1) % 3]) * s[(k + 2) % 3]); 
                W(i, Tri(j, k)) += w;
            }
        }
    }
}