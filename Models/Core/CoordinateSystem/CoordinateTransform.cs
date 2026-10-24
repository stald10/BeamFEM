using System;
using System.Collections.Generic;
using System.Text;
using Telma.Extensions;

namespace BeamFEM.Models.Core.CoordinateSystem;

public interface ICoordinateTransform<TSource, TTarget>
    where TSource : IVectorBase<TSource>
    where TTarget : IVectorBase<TTarget>
{
    static abstract bool IsLinear { get; }

    TTarget Transform(TSource sourcePoint);
    TSource InverseTransform(TTarget targetPoint);
    double Jacobian(TTarget targetPoint);
    IJacobyMatrix<TTarget, TSource> InverseJacoby();
}
