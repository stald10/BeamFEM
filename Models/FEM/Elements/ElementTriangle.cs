using BeamFEM.Models.Core.IBasis;
using BeamFEM.Models.Core.IElement;
using BeamFEM.Models.FEM.Basis;
using BeamFEM.Models.FEM.DofManager;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using Telma;

namespace BeamFEM.Models.FEM.Elements
{
    public class ElementTriangle : IFiniteElement
    {
        public Vector2D[] points { get => _p; }
        public double[] globalDof { get => _globalD; }

        public ElementTriangle(Vector2D[] _points, double[] _globalDof) 
        {
            _p = _points;
            _globalD = _globalDof;
        }

        private Vector2D[] _p;
        private double[] _globalD;

        public IBasisSet2D basis() => new BasisSet<Vector2D>
            (
                    Quadratures.TriangleOrder18,
                    TriangleBasis.Hierarchical.L3,
                    TriangleBasis.Hierarchical.L1,
                    TriangleBasis.Hierarchical.L2,
                    TriangleBasis.Hierarchical.L3L1,
                    TriangleBasis.Hierarchical.L1L2,
                    TriangleBasis.Hierarchical.L2L3,
                    new OrientedBasisFunction2D(TriangleBasis.Hierarchical.L3L1L3SubL1),
                    new OrientedBasisFunction2D(TriangleBasis.Hierarchical.L1L2L1SubL2),
                    new OrientedBasisFunction2D(TriangleBasis.Hierarchical.L2L3L2SubL3),
                    TriangleBasis.Hierarchical.L1L2L3
            );

        public sealed class Dof(IBasisSet2D basis) : DofManagerSys(dofCount: 10)
        {
            private readonly IBasisSet2D _basis = basis;

            public override int NumberOfDofOnVertex => 1;
            public override int NumberOfDofOnEdge => 2;
            public override int NumberOfDofOnElement => 1;

            public override void SetVertexDof(int localVertexIndex, int n, int dofIndex)
            {
                AssertIsValidVertexDofNumber(n);
                Debug.Assert(0 <= localVertexIndex && localVertexIndex < 3);
                _dof[localVertexIndex] = dofIndex;
            }

            public override void SetEdgeDof(int localEdgeIndex, bool isOrientationFlipped, int n, int dofIndex)
            {
                AssertIsValidEdgeDofNumber(n);
                Debug.Assert(0 <= localEdgeIndex && localEdgeIndex < 3);

                int basisIndex = 3 + 3 * n + localEdgeIndex;
                _dof[basisIndex] = dofIndex;
                if (n == 1)
                    ((OrientedBasisFunction2D)_basis.Basis[basisIndex]).IsOrientationFlipped = isOrientationFlipped;
            }

            public override void SetElementDof(int n, int dofIndex)
            {
                AssertIsValidElementDofNumber(n);
                _dof[9] = dofIndex;
            }
        }

    }
}
