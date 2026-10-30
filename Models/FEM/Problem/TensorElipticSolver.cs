using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.Core.Matrix;
using BeamFEM.Models.Core.Solver;
using BeamFEM.Models.FEM.Assembly;
using BeamFEM.Models.FEM.Integrator;
using BeamFEM.Models.FEM.Mesh;
using System;
using System.Collections.Generic;
using System.Text;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Problem
{
    /// <summary>
    ///     -div(σ) = b in Ω
    ///     σ*n = g on Г_n
    ///     u = u-  on Г_d
    ///     
    ///     σ = λtr(ε)I + 2με
    ///     
    ///     ε = 1/2*(∇u + (∇u)^T)
    /// </summary>
    /// <typeparam name="TSpace">Domain function space type</typeparam>
    /// <typeparam name="TBoundary">Boundary function space type</typeparam>
    /// <typeparam name="TOps">Matrix operations type for the function space</typeparam>
    public sealed class TensorElipticSolver<TSpace, TBoundary, TOps>(
           IMatrixFactory matrixFactory,
           IIntegrator<TSpace, TBoundary, TOps> integrator,
           ISolver algebraicSolver
    )
         where TSpace : IVectorBase<TSpace>
         where TBoundary : IVectorBase<TBoundary>
         where TOps : IMatrixOperations<TSpace, TSpace, TOps>
    {
        private readonly IMatrixFactory _matrixFactory = matrixFactory;
        private readonly IIntegrator<TSpace, TBoundary, TOps> _integrator = integrator;
        private readonly ISolver _algebraicSolver = algebraicSolver;

        public StationarySolution<TSpace, TBoundary> Solve(
               EllipticProblem<TSpace> problem,
               ISolver.Params solverParams = new())
        {
            const double time = 0.0;
            var mesh = (IMeshWithBoundaries<TSpace, TBoundary>)problem.Mesh;
            var dofManager = DofManager.NumerateDof(mesh, problem.BoundaryConditions);
            var assembler = new Assembler<TSpace, TBoundary, TOps>(mesh, dofManager, _matrixFactory, _integrator);

            var globalMatrix = assembler.CreateGlobalMatrix();
            var rhsVector = assembler.CreateRhsVector();
            var solution = new double[dofManager.TotalDofCount];

            void AddLocalMatrixWithFixedLoadContribution(LocalMatrix local, ReadOnlySpan<int> indices)
            {
                globalMatrix.AddLocalMatrix(local, indices);
                assembler.AddFixedLoadContribution(local, indices, dofManager.AsFixedSpan(solution), rhsVector);
            }

            assembler.CalculateStiffness(
                id => problem.Materials[id].Lambda,
                AddLocalMatrixWithFixedLoadContribution
            );

            _algebraicSolver.Matrix = globalMatrix;
            _algebraicSolver.Solve(
                rhsVector,
                dofManager.AsFreeSpan(solution),
                solverParams
            );

            return new(
                mesh: mesh,
                coefficients: solution
            );
        }
    }
}
