#pragma once

#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace ibm {
    class LinearSolver {
    public:
        static void GMRES(VectorXd* out, void (*linearSystem)(VectorXd* out, VectorXd v), const VectorXd rhs, const VectorXd initial, const int dimVector, const int maxDimSubspace, const double maxTolerance);
    };
}