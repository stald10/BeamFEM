using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.FEM.Basis;
using BeamFEM.Models.FEM.Mesh;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using Telma;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Elements;

public static class FiniteElementExtensions 
{
    extension<TSpace, TBoundary>(IFiniteElementBase<TSpace, TBoundary> self)
        where TSpace : IVectorBase<TSpace>
        where TBoundary : IVectorBase<TBoundary>
    {
        /// <summary>
        /// Evaluates mesh function value at given point inside element
        /// </summary>
        /// <param name="point">Point defined in mesh-space coordinates</param>
        /// <param name="coefficients">Basis function coefficients indexed by global dof indices (q array)</param>
        public double Evaluate(ReadOnlySpan<double> coefficients, TSpace point)
        {
            Debug.Assert(self.ContainsPoint(point, epsilon: 1e-12));
            var ep = self.MasterElementCoordinateSystem.Transform(point); // element-space point
            return self.EvaluateAtElementSpace(coefficients, ep);
        }

        public Vector2D TensorEvaluate(ReadOnlySpan<double> coefficients, TSpace point) 
        {
            var ep = self.MasterElementCoordinateSystem.Transform(point); // element-space point
            return self.EvaluateAtElementTensorSpace(coefficients, ep);
        }

        /// <summary>
        /// Evaluates mesh function value at given point inside element
        /// </summary>
        /// <param name="point">Point defined in element-space (master element) coordinates</param>
        /// <param name="coefficients">Basis function coefficients indexed by global dof indices (q array)</param>
        public double EvaluateAtElementSpace(ReadOnlySpan<double> coefficients, TBoundary point)
        {
            double value = 0.0;
            for (int i = 0; i < self.DOF.Count; ++i)
            {
                double q = coefficients[self.DOF.Dof[i]];
                value += q * self.Basis[i].Value(point);
            }

            return value;
        }

        public Vector2D EvaluateAtElementTensorSpace(ReadOnlySpan<double> coefficients, TBoundary point) 
        {
            double u = 0.0;
            double v = 0.0;
            for (int i = 0; i < self.Basis.Length; ++i)
            {
                int indexU = self.DOF.Dof[2 * i];
                int indexV = self.DOF.Dof[2 * i + 1];

                double q1 = coefficients[indexU];
                double q2 = coefficients[indexV];

                double basis = self.Basis[i].Value(point);

                u += q1 * basis;
                v += q2 * basis;
            }

            return new Vector2D(u, v);
        }
    }
}

/// <summary>
/// Utility methods to inline finite element member calls
/// </summary>
public static class FiniteElementComposition
{
    extension<TSpace, TBoundary>(IFiniteElementBase<TSpace, TBoundary> self)
        where TSpace : IVectorBase<TSpace>
        where TBoundary : IVectorBase<TBoundary>
    {
        public IMesh<TSpace> Mesh => self.Geometry.Mesh;
        public ReadOnlySpan<int> Vertices => self.Geometry.Vertices;
        public IEnumerable<Edge> Edges => self.Geometry.Edges;
        public int EdgeCount => self.Geometry.EdgeCount;
        public bool ContainsPoint(TSpace point, double epsilon = 1e-12) => self.Geometry.ContainsPoint(point, epsilon);
        public ICoordinateTransform<TSpace, TBoundary> MasterElementCoordinateSystem => self.Geometry.MasterElementCoordinateSystem;

        public IEnumerable<Quadratures.Node<TBoundary>> Quadratures => self.BasisSet.Quadratures;
        public ReadOnlySpan<IBasisFunction<TBoundary>> Basis => self.BasisSet.Basis;

        public IEnumerable<TSpace> VertexCoords
        {
            get
            {
                for (int i = 0; i < self.Vertices.Length; ++i)
                    yield return self.Mesh[self.Vertices[i]];
            }
        }
    }

    extension<TSpace, TBoundary>(IBoundaryElement<TSpace, TBoundary> self)
        where TSpace : IVectorBase<TSpace>
        where TBoundary : IVectorBase<TBoundary>
    {
        public IMesh<TSpace> Mesh => self.Geometry.Mesh;
        public ReadOnlySpan<int> Vertices => self.Geometry.Vertices;
        public IEnumerable<Edge> Edges => self.Geometry.Edges;
        public int EdgeCount => self.Geometry.EdgeCount;
        public ICoordinateTransform<TSpace, TBoundary> MasterElementCoordinateSystem => self.Geometry.MasterElementCoordinateSystem;

        public IEnumerable<Quadratures.Node<TBoundary>> Quadratures => self.BasisSet.Quadratures;
        public ReadOnlySpan<IBasisFunction<TBoundary>> Basis => self.BasisSet.Basis;
    }

    extension(IDofManager self)
    {
        public int Count => self.Dof.Length;
    }

    extension<TSpace>(IElementGeometryBase<TSpace> self)
        where TSpace : IVectorBase<TSpace>
    {
        public int VertexCount => self.Vertices.Length;
    }
}