
#include <cmath>

#include "LinearSolver.h"


using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace ibm
{
    void LinearSolver::GMRES(VectorXd* out, void (*linearSystem)(VectorXd* out, VectorXd v), const VectorXd rhs, const VectorXd initial, const int dimVector, const int maxDimSubspace, const double maxTolerance)
    {
        int dimSubspace = maxDimSubspace;

        VectorXd* resBasis = new VectorXd[maxDimSubspace + 1];
        for (int i = 0; i < maxDimSubspace + 1; ++i)
        {
            resBasis[i].resize(dimVector);
        }
        MatrixXd hessianMatrix(maxDimSubspace + 1, maxDimSubspace);

        VectorXd sines(maxDimSubspace);
        VectorXd cosines(maxDimSubspace);
        VectorXd residual(maxDimSubspace + 1);

        VectorXd tmpVec(maxDimSubspace);
        linearSystem(&tmpVec, initial);


        resBasis[0] = rhs - tmpVec;

        residual(0) = resBasis[0].norm();

        resBasis[0] /= resBasis[0].norm();

        for (int j = 0; j < maxDimSubspace; ++j)
        {
            linearSystem(&resBasis[j + 1], resBasis[j]);
            for (int i = 0; i <= j; ++i)
            {
                hessianMatrix(i, j) = resBasis[i].dot(resBasis[j + 1]);
                resBasis[j + 1] -= hessianMatrix(i, j) * resBasis[i];
            }

            hessianMatrix(j + 1, j) = resBasis[j + 1].norm();
            resBasis[j + 1] /= hessianMatrix(j + 1, j);




            for (int i = 0; i < j; ++i)
            {
                double tmp = cosines(i) * hessianMatrix(i, j) + sines(i) * hessianMatrix(i + 1, j);
                hessianMatrix(i + 1, j) = -sines(i) * hessianMatrix(i, j) + cosines(i) * hessianMatrix(i + 1, j);
                hessianMatrix(i, j) = tmp;
            }

            cosines(j) = hessianMatrix(j, j) / std::sqrt(hessianMatrix(j, j) * hessianMatrix(j, j) + hessianMatrix(j + 1, j) * hessianMatrix(j + 1, j));
            sines(j) = hessianMatrix(j + 1, j) / std::sqrt(hessianMatrix(j, j) * hessianMatrix(j, j) + hessianMatrix(j + 1, j) * hessianMatrix(j + 1, j));

            hessianMatrix(j, j) = cosines(j) * hessianMatrix(j, j) + sines(j) * hessianMatrix(j + 1, j);
            hessianMatrix(j + 1, j) = 0.0;



            residual(j + 1) = -sines(j) * residual(j);
            residual(j) = cosines(j) * residual(j);

            double error = std::abs(residual(j + 1));


            if (error < maxTolerance)
            {
                dimSubspace = j + 1;
                break;
            }
        }

        //solve system and get y
        VectorXd y(dimSubspace);



        y = hessianMatrix.block(0, 0, dimSubspace, dimSubspace).triangularView<Eigen::Upper>().solve(residual.segment(0, dimSubspace));

        *out = initial;
        for (int j = 0; j < dimSubspace; ++j)
        {
            *out += resBasis[j] * y(j);
        }

        delete[] resBasis;
    }
}
