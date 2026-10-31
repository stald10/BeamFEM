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
                var mp = masterCs.InverseTransform(ep); // mesh-space point
                var psiJ = element.Basis[j].Value(ep);
                var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp));  // FIXME: may be optimized (by IsConstant)
                value += source(mp) * psiJ * q.Weight * jacobian;
            }
            outLoad[j] = value;
        }
    }

    private static void CalculateTensorLocalLoad<TB>(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, TSpace> source,
        Span<double> outLoad,
        Func<TSpace, double> thinkness) where TB : IVectorBase<TB>
    {
        int n = element.DOF.Count;
        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

    }

    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, double> gamma)
                => throw new NotSupportedException();
    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> gamma)
                => throw new NotSupportedException();
    public LocalMatrix CalculateLocalStiffness(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> lambda)
                => throw new NotSupportedException();

    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
        => CalculateTensorLocalLoad<TSpace>(element, source, outLoad, thinkness);


    public void CalculateTensorLocalLoad(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, TSpace> source, Span<double> outLoad, Func<TSpace, double> thinkness)
    {
        throw new NotImplementedException();
    }

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
            var jacobian = Math.Abs(masterCs.Jacobian(ep) * meshCs.Jacobian(mp));

            var matrixB = new TensorMatrix(3, n);

            for (int i = 0; i < n / 2; ++i)
            {
                var shape = element.Basis[i].Derivatives(ep);
                var shapeNatural = invJ * shape;

                if (shapeNatural is Vector2D deriv)
                {

                    matrixB[0, 2 * i] = deriv.X;
                    matrixB[1, 2 * i + 1] = deriv.Y;
                    matrixB[2, 2 * i + 1] = deriv.X;
                    matrixB[2, 2 * i] = deriv.Y;

                }
            }

            var matrixBT = matrixB.Tranpose();
            var BTDB = matrixBT * matrixD * matrixB;


            //лямбда выступает в роли толщины, чтобы не переделывать параметры метода
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

                    var str = gradPhiI.GetType();

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

