using BeamFEM.Models.FEM.Basis;
using BeamFEM.Models.Core.IBasis;
using System;
using System.Collections.Generic;
using System.Text;
using Telma;

namespace BeamFEM.Models.FEM.Elements;

public static class SegmentBasis
{
    /// <summary>N0(ξ) = 1 - ξ</summary>
    public static readonly IBasisFunction1D N0 = new SegmentN0();
    /// <summary>N1(ξ) = ξ</summary>
    public static readonly IBasisFunction1D N1 = new SegmentN1();
    /// <summary>N0N1(ξ) = ξ (1 - ξ)</summary>
    public static readonly IBasisFunction1D N0N1 = new SegmentN0N1();
    /// <summary>N0N1N0SubN1(ξ) = ξ (1 - ξ)(2ξ - 1)</summary>
    public static readonly IBasisFunction1D N0N1N0SubN1 = new SegmentN0N1N0SubN1();

    private readonly struct SegmentN0 : IBasisFunction1D
    {
        public double Value(Vector1D point) => 1 - point;
        public Vector1D Derivatives(Vector1D _) => -1;
    }

    private readonly struct SegmentN1 : IBasisFunction1D
    {
        public double Value(Vector1D point) => point;
        public Vector1D Derivatives(Vector1D _) => 1;
    }

    private readonly struct SegmentN0N1 : IBasisFunction1D
    {
        public double Value(Vector1D point) => point * (1 - point);
        public Vector1D Derivatives(Vector1D point) => 1 - 2 * point;
    }

    private readonly struct SegmentN0N1N0SubN1 : IBasisFunction1D
    {
        public double Value(Vector1D point) => point * (1 - point) * (2 * point - 1);
        public Vector1D Derivatives(Vector1D point) => -6 * point * (point - 1) - 1;
    }
}
