using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.Core.Matrix;
using BeamFEM.Models.FEM.Elements;
using System;
using System.Collections.Generic;
using System.Text;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Integrator;

/// <summary>
/// Performs numerical integration for finite element assembly over domain and boundary elements.
/// </summary>
/// <typeparam name="TSpace">Domain function space type</typeparam>
/// <typeparam name="TBoundary">Boundary function space type</typeparam>
/// <typeparam name="TOps">Matrix operations type</typeparam>
public interface IIntegrator<TSpace, TBoundary, TOps>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
    where TOps : IMatrixOperations<TSpace, TSpace, TOps>
{
    /// <summary>
    /// Computes local stiffness matrix: G_ij = ∫_Ω λ (∇φ_i)·(∇φ_j) dΩ
    /// </summary>
    /// <param name="element">Domain element</param>
    /// <param name="lambda">Material property λ(x) defined in mesh-space coordinates</param>
    /// <returns>Local stiffness matrix of size n × n, where n = <paramref name="element"/>.DOF.Count</returns>
    LocalMatrix CalculateLocalStiffness(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> lambda
    );

    /// <summary>
    /// Computes local mass matrix: M_ij = ∫_S γ φ_i φ_j dS
    /// </summary>
    /// <param name="element">Boundary element</param>
    /// <param name="gamma">Coefficient function γ(x) defined in mesh-space coordinates</param>
    /// <returns>Local mass matrix of size n × n, where n = <paramref name="element"/>.DOF.Count</returns>
    LocalMatrix CalculateLocalMass(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, double> gamma
    );

    /// <summary>
    /// Computes local mass matrix: M_ij = ∫_Ω γ φ_i φ_j dΩ
    /// </summary>
    /// <param name="element">Domain element</param>
    /// <param name="gamma">Coefficient function γ(x) defined in mesh-space coordinates</param>
    /// <returns>Local mass matrix of size n × n, where n = <paramref name="element"/>.DOF.Count</returns>
    LocalMatrix CalculateLocalMass(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> gamma
    );

    /// <summary>
    /// Computes local load vector: b_i = ∫_S s φ_i dS
    /// </summary>
    /// <param name="element">Boundary element</param>
    /// <param name="source">Source function s(x) defined in mesh-space coordinates</param>
    /// <param name="outLoad">Output vector of size n, where n = <paramref name="element"/>.DOF.Count</param>
    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    );


    /// <summary>
    /// Computes local load vector: b_i = ∫_Ω s φ_i dΩ
    /// </summary>
    /// <param name="element">Domain element</param>
    /// <param name="source">Source function s(x) defined in mesh-space coordinates</param>
    /// <param name="outLoad">Output vector of size n, where n = length of <see cref="IFiniteElementBase{TSpace, TSpace}.DOF"/></param>
    public void CalculateLocalLoad(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, double> source,
        Span<double> outLoad
    );


    /// <summary>
    /// Вычисление локального вектора f_1^e = ∫_Ω^e h*N^T*b dΩ^e
    /// </summary>
    /// <param name="element">Элемент</param>
    /// <param name="source">вектор b</param>
    /// <param name="outLoad">Результат</param>
    public void CalculateTensorLocalLoad(
        IFiniteElementBase<TSpace, TSpace> element,
        Func<TSpace, TSpace> source,
        Span<double> outLoad,
        Func<TSpace, double> thinkness
    );

    /// <summary>
    /// Вычисление вкалада краевого: f_2^e = ∫_Г^e h*N^T*t dГ^e
    /// </summary>
    /// <param name="element">Граничный элемент</param>
    /// <param name="source">вектор нагрузок t</param>
    /// <param name="outLoad">Результат</param>
    /// <param name="thinkness">Толщина пластины</param>
    public void CalculateTensorLocalLoad(
        IFiniteElementBase<TSpace, TBoundary> element,
        Func<TSpace, TSpace> source,
        Span<double> outLoad,
        Func<TSpace, double> thinkness
    );


    LocalMatrix CalculateLocalTensorStiffness(
    IFiniteElementBase<TSpace, TSpace> element,
    Func<TSpace, double> thinkness,
    double E,
    double poisson
);




}
