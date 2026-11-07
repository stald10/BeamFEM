using Telma.Extensions;
using BeamFEM.Models.FEM.Elements;
using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.Core.Matrix;
using System;
using Telma;
using System.Reflection.Metadata;
using System.Reflection.Metadata.Ecma335;

namespace BeamFEM.Models.FEM.Integrator;

public class NumericIntegratorTensor<TSpace, TBoundary, TOps> : IIntegrator<TSpace, TBoundary, TOps>
        where TSpace : IVectorBase<TSpace>
        where TBoundary : IVectorBase<TBoundary>
        where TOps : IMatrixOperations<TSpace, TSpace, TOps>
{
    private NumericIntegratorTensor() { }
    public static readonly NumericIntegratorTensor<TSpace, TBoundary, TOps> Instance = new();


    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    ) => throw new NotSupportedException();

    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    ) => throw new NotSupportedException();

    private static void CalculateTensorLocalLoad<TB>(
        IFiniteElementBase<TSpace, TB> element,
        Func<TSpace, TSpace> source,
        Span<double> outLoad,
        Func<TSpace, double> thinkness) where TB : IVectorBase<TB>
    {
        int n = element.DOF.Count;
        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        outLoad.Clear();

        foreach (var q in element.Quadratures)
        {
            var ep = q.Point;
            var mp = masterCs.InverseTransform(ep);
            var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp));


            double Weight = thinkness(mp) * jacobian * q.Weight;

            if (source(mp) is Vector2D B) 
            {
                for (int i = 0; i < n / 2; ++i)
                {
                    var shape = element.Basis[i].Value(ep);
                    outLoad[2 * i] += shape * Weight * B.X;
                    outLoad[2 * i + 1] += shape * Weight* B.Y;
                }
            }
        }
    }

    private static LocalMatrix CalculateLocalMass<TB>(
     IFiniteElementBase<TSpace, TB> element,
     Func<TSpace, double> gamma
)
    where TB : IVectorBase<TB>
    {
        int n = element.DOF.Count;

        

        var mass = new LocalMatrix(n);

        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;
        foreach (var q in element.Quadratures)
        {
            var ep = q.Point;
            var mp = masterCs.InverseTransform(ep);

            var jac1 = masterCs.Jacobian(ep);
            var jac2 = meshCs.Jacobian(mp);

            var jacobian = jac1 * jac2;
            //var jacobian = Math.Abs(masterCs.Jacobian(ep));

            var matrixN = new TensorMatrix(2, n);

            for (int i = 0; i < element.Basis.Length; ++i) 
            {
                int dofU = 2 * i;
                int dofV = 2 * i + 1;

                matrixN[0, dofU] = element.Basis[i].Value(ep);
                matrixN[1, dofV] = element.Basis[i].Value(ep);

            }

            var matrixNTran = matrixN.Tranpose();

            var NTN = matrixNTran * matrixN;

            var Weight = jacobian * q.Weight;
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    var val = NTN[i, j] * Weight;
                    mass[i, j] += val;
                    if (i != j)
                    {
                        mass[j, i] += val;
                    }
                }
            }
        }
        return mass;
    }


    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, double> gamma)
                => CalculateLocalMass<TBoundary>(element, gamma);
    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> gamma)
                => throw new NotSupportedException();
    public LocalMatrix CalculateLocalStiffness(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> lambda)
                => throw new NotSupportedException();

    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
        => CalculateTensorLocalLoad<TSpace>(element, source, outLoad, thinkness);


    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
            => CalculateTensorLocalLoad<TBoundary>(element, source, outLoad, thinkness);

    public LocalMatrix CalculateLocalTensorStiffness(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> thinkness, double E, double poisson)
    {
        int n = element.DOF.Count;

        //double E = 1e+5;
        //double poisson = 0.3;

        double _E = E / (1 - poisson * poisson);

        var matrixD = new TensorMatrix(3, 3);


        var stiffness = new LocalMatrix(n);

        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        var masterInvJ = masterCs.InverseJacoby();
        var meshInvJ = meshCs.InverseJacoby();

        matrixD[0, 0] = 1;
        matrixD[0, 1] = poisson;
        matrixD[1, 0] = matrixD[0, 1];
        matrixD[1, 1] = matrixD[0, 0];
        matrixD[2, 2] = 0;

        matrixD *= _E;

        matrixD[2, 2] = E / (2 * (1 + poisson));

        foreach (var q in element.Quadratures)
        {
            var ep = q.Point;
            var mp = masterCs.InverseTransform(ep);
            var invJ = masterInvJ.MulAt<TSpace, TSpace, TSpace, TOps>(ep, meshInvJ, mp);

            var jac1 = masterCs.Jacobian(ep);
            var jac2 = meshCs.Jacobian(mp);

            var jacobian = jac1 * jac2;
            //var jacobian = Math.Abs(masterCs.Jacobian(ep));

            var matrixB = new TensorMatrix(3, n);


            for (int i = 0; i < 3; ++i)
            {

                if (element.Basis[i].Derivatives(ep) * invJ is Vector2D deriv) 
                {
                    int dofU = 2 * i;
                    int dofV = 2 * i + 1;

                    matrixB[0, dofU] = deriv.X;
                    matrixB[1, dofV] = deriv.Y;
                    matrixB[2, dofU] = deriv.Y;
                    matrixB[2, dofV] = deriv.X;
                }


            }
            for (int e = 0; e < 3; ++e)
            {
                int quadBasisIdx = 3 + e;
                int cubicBasisIdx = 6 + e;
                int edgeStartPos = 6 + 4 * e;


                if (element.Basis[quadBasisIdx].Derivatives(ep) * invJ is Vector2D derivQuad) 
                {
                    matrixB[0, edgeStartPos + 0] = derivQuad.X;
                    matrixB[1, edgeStartPos + 1] = derivQuad.Y;
                    matrixB[2, edgeStartPos + 0] = derivQuad.Y;
                    matrixB[2, edgeStartPos + 1] = derivQuad.X;
                }
                if ( element.Basis[cubicBasisIdx].Derivatives(ep) * invJ is Vector2D derivCubic) 
                {
                    matrixB[0, edgeStartPos + 2] = derivCubic.X;
                    matrixB[1, edgeStartPos + 3] = derivCubic.Y;
                    matrixB[2, edgeStartPos + 2] = derivCubic.Y;
                    matrixB[2, edgeStartPos + 3] = derivCubic.X;
                }

            }

            {

                int centerStartPos = 18;

                if (element.Basis[9].Derivatives(ep)*invJ is Vector2D derivCenter)
                {
                    matrixB[0, centerStartPos + 0] = derivCenter.X;
                    matrixB[1, centerStartPos + 1] = derivCenter.Y;
                    matrixB[2, centerStartPos + 0] = derivCenter.Y;
                    matrixB[2, centerStartPos + 1] = derivCenter.X;
                }
            }

            var matrixBT = matrixB.Tranpose();
            var BTDB = matrixBT * matrixD * matrixB;

            double Weight = thinkness(mp) * jacobian * q.Weight;

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    var val = BTDB[i, j] * Weight;
                    stiffness[i, j] += val;
                    if (i != j)
                    {
                        stiffness[j, i] += val;
                    }
                }
            }
        }

        return stiffness;
    }
}



public class NumericIntegrator<TSpace, TBoundary, TOps> : IIntegrator<TSpace, TBoundary, TOps>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
    where TOps : IMatrixOperations<TSpace, TSpace, TOps>
{
    public static readonly NumericIntegrator<TSpace, TBoundary, TOps> Instance = new();
    private NumericIntegrator() { }

    public LocalMatrix CalculateLocalStiffness(
       IFiniteElementBase<TSpace, TSpace> element,
       Func<TSpace, double> lambda
    )
    {
        int n = element.DOF.Count;
        var stiffness = new LocalMatrix(n);

        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        var masterInvJ = masterCs.InverseJacoby();
        var meshInvJ = meshCs.InverseJacoby();

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                var value = 0.0;
                foreach (var q in element.Quadratures)
                {
                    var ep = q.Point; // master element-space point
                    var mp = masterCs.InverseTransform(ep); // mesh-space point
                    var invJ = masterInvJ.MulAt<TSpace, TSpace, TSpace, TOps>(ep, meshInvJ, mp);
                    var gradPhiI = element.Basis[i].Derivatives(ep);
                    var gradPhiJ = element.Basis[j].Derivatives(ep);

    

                    // product = grad(phi_i)^T * J^(-T) * J(-1) * grad(phi_j)
                    var product = (gradPhiI * invJ) * (gradPhiJ * invJ);

                    //var jacobian = Math.Abs(1.0 / invJ.Det());
                    var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp));
                    value += lambda(mp) * product * jacobian * q.Weight;
                }
                stiffness[i, j] = value;
                stiffness[j, i] = value;
            }
        }
        return stiffness;
    }

    public LocalMatrix CalculateLocalMass(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, double> gamma
    ) => CalculateLocalMass<TBoundary>(element, gamma);

    public LocalMatrix CalculateLocalMass(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> gamma
    ) => CalculateLocalMass<TSpace>(element, gamma);

    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    ) => CalculateLocalLoad<TBoundary>(element, source, outLoad);

    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    ) => CalculateLocalLoad<TSpace>(element, source, outLoad);

    private static LocalMatrix CalculateLocalMass<TB>(
         IFiniteElementBase<TSpace, TB> element,
         Func<TSpace, double> gamma
    )
        where TB : IVectorBase<TB>
    {
        int n = element.DOF.Count;
        var mass = new LocalMatrix(n);

        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                var value = 0.0;
                foreach (var q in element.Quadratures)
                {
                    var ep = q.Point; // master element-space point
                    var mp = masterCs.InverseTransform(ep); // mesh-space point
                    var phiI = element.Basis[i].Value(ep);
                    var phiJ = element.Basis[j].Value(ep);
                    var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp)); // FIXME: may be optimized (by IsConstant)
                    value += gamma(mp) * phiI * phiJ * jacobian * q.Weight;
                }
                mass[i, j] = value;
                mass[j, i] = value;
            }
        }
        return mass;
    }

    private static void CalculateLocalLoad<TB>(
        IFiniteElementBase<TSpace, TB> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    )
        where TB : IVectorBase<TB>
    {
        int n = element.DOF.Count;
        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        for (int j = 0; j < n; j++)
        {
            var value = 0.0;
            foreach (var q in element.Quadratures)
            {
                var ep = q.Point; // master element-space point
                var mp = masterCs.InverseTransform(q.Point); // mesh-space point
                var psiJ = element.Basis[j].Value(ep);
                var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp));  // FIXME: may be optimized (by IsConstant)
                value += source(mp) * psiJ * q.Weight * jacobian;
            }
            outLoad[j] = value;
        }
    }

    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
           => throw new NotSupportedException();

    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
           => throw new NotSupportedException();

    public LocalMatrix CalculateLocalTensorStiffness(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> thinkness, double E, double poisson)
        => throw new NotSupportedException();
}

