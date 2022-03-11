#define _USE_MATH_DEFINES

#include <cmath>
#include "Kernel1D.h"


namespace ibm
{
    double Kernel1D::Cos(double x)
    {
        if (x < -2)
        {
            return 0.0;
        }

        if (x > 2)
        {
            return 0.0;
        }

        return (std::cos(x * M_PI_2) + 1) * 0.25;
    }
}