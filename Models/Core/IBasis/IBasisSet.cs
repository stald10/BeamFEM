using BeamFEM.Models.FEM.Basis;
using System;
using System.Collections.Generic;
using System.Text;
using Telma;
using Telma.Extensions;

namespace BeamFEM.Models.Core.IBasis;

public interface IBasisSet<TSpace> where TSpace : IVectorBase<TSpace>
{
    IEnumerable<Quadratures.Node<TSpace>> Quadratures { get; }
    ReadOnlySpan<IBasisFunction<TSpace>> Basis { get; }
}


public readonly struct BasisSet<TSpace>(
    Func<IEnumerable<Quadratures.Node<TSpace>>> quadratures,
    params IBasisFunction<TSpace>[] basis
) : IBasisSet<TSpace> where TSpace : IVectorBase<TSpace>
{
    public IEnumerable<Quadratures.Node<TSpace>> Quadratures => quadratures();
    private readonly IBasisFunction<TSpace>[] _basis = basis;
    public ReadOnlySpan<IBasisFunction<TSpace>> Basis => _basis;
}

