using BeamFEM.Models.FEM.Mesh;
using BeamFEM.Models.FEM.Elements;
using BeamFEM.Models.Core.CoordinateSystem;
using System;
using System.Collections.Generic;
using System.Text;
using Telma.Extensions;
using Telma;

namespace BeamFEM.Models.FEM.Problem;

public sealed class StationarySolution<TSpace, TBoundary>(
    IMeshWithBoundaries<TSpace, TBoundary> mesh,
    double[] coefficients
)
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
{
    public IMeshWithBoundaries<TSpace, TBoundary> Mesh { get; } = mesh;
    private readonly double[] _coefficients = coefficients;
    public ReadOnlySpan<double> Coefficients => _coefficients;

    private SearchTree<TSpace>? _searchTree = null;

    public Vector2D TensorEvalute(TSpace point, CoordinateSpace coordSpace = CoordinateSpace.Mesh) 
    {
        point = coordSpace switch
        {
            CoordinateSpace.Physical => Mesh.CoordinateSystem.Transform(point),
            CoordinateSpace.Mesh => point,
            _ => throw new NotSupportedException()
        };
        _searchTree ??= SearchTree<TSpace>.BuildFor(Mesh);
        var element = _searchTree.FindElementAt(point) ??
            throw new ArgumentException("Failed to find element at the point", nameof(point));
        return element.TensorEvaluate(_coefficients,point);
    }

    public double Evaluate(TSpace point, CoordinateSpace coordSpace = CoordinateSpace.Mesh)
    {
        point = coordSpace switch
        {
            CoordinateSpace.Physical => Mesh.CoordinateSystem.Transform(point),
            CoordinateSpace.Mesh => point,
            _ => throw new NotSupportedException()
        };
        _searchTree ??= SearchTree<TSpace>.BuildFor(Mesh);
        var element = _searchTree.FindElementAt(point) ??
            throw new ArgumentException("Failed to find element at the point", nameof(point));
        return element.Evaluate(_coefficients, point);
    }

    /// <summary>
    /// Numerically calculates L2 norm difference with <paramref name="analyticalSolution"/>
    /// </summary>
    /// <param name="analyticalSolution">Function defined in mesh-space coordinates</param>
    public double Difference(Func<TSpace, double> analyticalSolution)
    {
        double totalDifference = 0;
        foreach (var element in Mesh.FiniteElements)
        {
            var masterCs = element.MasterElementCoordinateSystem;

            double elementDifference = 0;
            foreach (var q in element.Quadratures)
            {
                var mp = masterCs.InverseTransform(q.Point); // mesh-space point
                var analyticalValue = analyticalSolution(mp);
                var femValue = element.EvaluateAtElementSpace(_coefficients, q.Point);
                var difference = analyticalValue - femValue;
                elementDifference += q.Weight * difference * difference * Math.Abs(masterCs.Jacobian(q.Point));
            }
            totalDifference += elementDifference;
        }
        return totalDifference;
    }

    public double DifferenceTensor(Func<TSpace, Vector2D> analyticSolution) 
    {
        double totalDifference = 0;
        foreach (var element in Mesh.FiniteElements)
        {
            var masterCs = element.MasterElementCoordinateSystem;

            double elementDifference = 0;
            foreach (var q in element.Quadratures)
            {
                var mp = masterCs.InverseTransform(q.Point); // mesh-space point
                var analyticalValue = analyticSolution(mp);
                var femValue = element.EvaluateAtElementTensorSpace(_coefficients, q.Point);
                var difference = analyticalValue - femValue;
                elementDifference += q.Weight * difference * difference * Math.Abs(masterCs.Jacobian(q.Point));
            }
            totalDifference += elementDifference;
        }
        return totalDifference;
    } 
}

