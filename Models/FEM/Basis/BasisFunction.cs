using Splat.ModeDetection;
using System;
using System.Collections.Generic;
using System.Text;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Basis;



public interface IBasisFunction<TSpace> where TSpace : IVectorBase<TSpace>
{
    double Value(TSpace point);
    TSpace Derivatives(TSpace point);
}


public sealed class OrientedBasisFunction<TSpace>(IBasisFunction<TSpace> basis) :
    IBasisFunction<TSpace>
    where TSpace : IVectorBase<TSpace>
{
    private int _sign = 1;
    private readonly IBasisFunction<TSpace> _basis = basis;

    public bool IsOrientationFlipped { get => _sign == -1; set => _sign = value ? -1 : 1; }

    public double Value(TSpace point) => _basis.Value(point) * _sign;
    public TSpace Derivatives(TSpace point) => _basis.Derivatives(point) * _sign;
}
