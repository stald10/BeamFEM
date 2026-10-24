using System;
using System.Collections.Generic;
using System.Text;

namespace BeamFEM.Models.Core.Util;

public static class MathUtils
{
    public static int ToleratedSign(double value, double tolerance)
    {
        if (value < -tolerance) return -1;
        if (value > tolerance) return 1;
        return 0;
    }
}
