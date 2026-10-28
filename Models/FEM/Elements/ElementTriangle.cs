using BeamFEM.Models.Core.IBasis;
using BeamFEM.Models.Core.IElement;
using BeamFEM.Models.FEM.Basis;
using System.Diagnostics;
using System.Text;
using Telma;

namespace BeamFEM.Models.FEM.Elements
{

    public class HierarchicalCubicTriangleTensorFactory : IFiniteElementFactory2D
    {
        public IFiniteElement2D CreateElement(IMesh2D mesh, int[] vertices, int materialIndex)
        {
            var basis = DefaultBasis();
            return new FiniteElement2D(
                Geometry: new TriangleGeometry(vertices) { Mesh = mesh },
                DOF: new Dof(basis),
                BasisSet: basis,
                MaterialIndex: materialIndex
            );
        }

        public static IBasisSet2D DefaultBasis() => new BasisSet2D(
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

        public sealed class Dof(IBasisSet2D basis) : DofManager(dofCount: 20)
        {
            private readonly IBasisSet2D _basis = basis;

            public override int NumberOfDofOnVertex => 2;
            public override int NumberOfDofOnEdge => 4;
            public override int NumberOfDofOnElement => 2;

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


    public class HierarchicalCubicTriangleFactory : IFiniteElementFactory2D
    {

        public IFiniteElement2D CreateElement(IMesh2D mesh, int[] vertices, int materialIndex)
        {
            var basis = DefaultBasis();
            return new FiniteElement2D(
                Geometry: new TriangleGeometry(vertices) { Mesh = mesh },
                DOF: new Dof(basis),
                BasisSet: basis,
                MaterialIndex: materialIndex
            );
        }

        public static IBasisSet2D DefaultBasis() => new BasisSet2D(
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

        public sealed class Dof(IBasisSet2D basis) : DofManager(dofCount: 10)
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
