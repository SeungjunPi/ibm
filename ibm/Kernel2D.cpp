#include "Kernel2D.h"
#include "Kernel1D.h"

namespace ibm
{
    double Kernel2D::HexaLinear(const double x, const double y)
    {
        if (x <= 0 && x >= -1)
        {
            if (y >= 0 && y <= x + 1)
            {
                return x - y + 1;
            }

            if (y < 0 && y >= x)
            {
                return x + 1;
            }

            if (y >= -1 && y < x)
            {
                return y + 1;
            }
        }

        if (x > 0 && x <= 1)
        {
            if (y <= 1 && y >= x)
            {
                return 1 - y;
            }

            if (y >= 0 && y < x)
            {
                return 1 - x;
            }

            if (y < 0 && y >= x)
            {
                return y - x + 1;
            }
        }

        return 0.0;
    }

    double Kernel2D::Cos(double x, double y)
    {
        return Kernel1D::Cos(x) * Kernel1D::Cos(y);
    }
}
