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


    public void CalculateLocalLoad(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, double> source, Span<double> outLoad)
    {
        throw new NotImplementedException();
    }

    public void CalculateLocalLoad(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> source, Span<double> outLoad)
    {
        throw new NotImplementedException();
    }

    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TBoundary> element, Func<TSpace, double> gamma)
    {
        throw new NotImplementedException();
    }

    public LocalMatrix CalculateLocalMass(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> gamma)
    {
        throw new NotImplementedException();
    }

    public LocalMatrix CalculateLocalStiffness(IFiniteElementBase<TSpace, TSpace> element, Func<TSpace, double> lambda)
    {
        int n = element.DOF.Count;
        var stiffness = new LocalMatrix(n);

        var masterCs = element.MasterElementCoordinateSystem;
        var meshCs = element.Mesh.CoordinateSystem;

        var masterInvJ = masterCs.InverseJacoby();
        var meshInvJ = meshCs.InverseJacoby();


        double[][] C = new double[3][];
        for (int i = 0; i < 3; ++i)
        {
            C[i] = new double[3];
        }

        double E = 1e+6;
        double poisson = 1e+3;

        double _E = E / (1 - poisson * poisson);

        C[0][0] = _E;
        C[0][1] = poisson * _E;

        C[1][0] = C[0][1];
        C[1][1] = C[0][0];

        C[2][2] = E / (2 * (1 + poisson));

        double[][] B = new double[3][];

        for (int i = 0; i < 3; ++i)
        {
            B[i] = new double[n];
        }

        for (int i = 0; i < n/2; ++i)
        {
            double valX = 0;
            double valY = 0;
            foreach (var q in element.Quadratures)
            {
                var ep = q.Point;
                var shape = element.Basis[i].Derivatives(ep);

                if (shape is Vector2D deriv) {valX += deriv.X; valY += deriv.Y;}

            }

            B[0][2 * i] = valX;
            B[1][2 * i + 1] = valY;

            B[2][2 * i + 1] = valX;
            B[2][2 * i] = valY;
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


}

