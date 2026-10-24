using BeamFEM.Models.Core.Matrix;
using System;
using System.Collections.Generic;
using System.Text;

namespace BeamFEM.Models.Core.Solver;


public interface ISolver : IDisposable
{
    public IGlobalMatrix Matrix { set; }

    public record struct Params(
        double Eps = 1e-10,
        int MaxIterations = 1000
    );

    (double residual, int iterations) Solve(ReadOnlySpan<double> rhsVector, Span<double> solution, Params paramz = new());
}
