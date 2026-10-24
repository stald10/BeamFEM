using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.FEM.Mesh;
using System;
using System.Collections.Generic;
using System.Text;
using Telma;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Elements;

public class SegmentGeometry<TSpace> : ElementGeometry<TSpace, Vector1D>
    where TSpace : struct, IVectorBase<TSpace>
{
    public SegmentGeometry(int[] vertexIndices) : base(vertexIndices)
    {
        if (vertexIndices.Length != 2)
            throw new ArgumentException($"Expected 2 vertices for {nameof(SegmentGeometry<>)} but got {vertexIndices.Length}");
    }

    public override IEnumerable<Edge> Edges => [new(Vertices[0], Vertices[1])];
    public override int EdgeCount => 1;

    public override bool ContainsPoint(TSpace point, double epsilon = 1E-12)
    {
        if (TSpace.Dimensions != 1) return false;
        var xi = MasterElementCoordinateSystem.Transform(point);
        return -epsilon <= xi && xi <= 1 + epsilon;
    }

    public override ICoordinateTransform<TSpace, Vector1D> MasterElementCoordinateSystem =>
        new SegmentParametrization<TSpace>(Mesh[Vertices[0]], Mesh[Vertices[1]]);
}

public sealed class TriangleGeometry : VolumeElementGeometry2D
{
    public TriangleGeometry(int[] vertexIndices) : base(vertexIndices)
    {
        if (vertexIndices.Length != 3)
            throw new ArgumentException($"Expected 3 vertices for {nameof(TriangleGeometry)} but got {vertexIndices.Length}");
    }

    public override IEnumerable<Edge> Edges
    {
        get
        {
            yield return new(Vertices[0], Vertices[1]);
            yield return new(Vertices[1], Vertices[2]);
            yield return new(Vertices[2], Vertices[0]);
        }
    }

    public override int EdgeCount => 3;

    public override bool ContainsPoint(Vector2D point, double epsilon = 1E-12)
    {
        var (xi, eta) = MasterElementCoordinateSystem.Transform(point);
        return
            xi >= -epsilon &&
            eta >= -epsilon &&
            (xi + eta) <= 1 + epsilon;
    }

    public override ICoordinateTransform<Vector2D, Vector2D> MasterElementCoordinateSystem =>
        new BarycentricCoordinateSystem(
            Mesh[Vertices[0]], Mesh[Vertices[1]], Mesh[Vertices[2]]
        );
}


public readonly record struct Edge(int I, int J)
{
    public Edge Sorted() => Sorted(out _);
    public Edge Sorted(out bool isOrientationFlipped)
    {
        isOrientationFlipped = I > J;
        return isOrientationFlipped ? new(J, I) : this;
    }
}

public interface IElementGeometry<TSpace>
    where TSpace : IVectorBase<TSpace>
{
    IMesh<TSpace> Mesh { get; }

    ReadOnlySpan<int> Vertices { get; }
    IEnumerable<Edge> Edges { get; }
    int EdgeCount { get; }
}

public interface IElementGeometryBase<TSpace>
    where TSpace : IVectorBase<TSpace>
{
    IMesh<TSpace> Mesh { get; }
    ReadOnlySpan<int> Vertices { get; }
    IEnumerable<Edge> Edges { get; }
    int EdgeCount { get; }

    bool ContainsPoint(TSpace point, double epsilon = 1e-12);
}

public interface IElementGeometry<TSpace, TBoundary> : IElementGeometryBase<TSpace>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
{
    ICoordinateTransform<TSpace, TBoundary> MasterElementCoordinateSystem { get; }
}

public abstract class ElementGeometry<TSpace, TBoundary>(int[] vertexIndices) : IElementGeometry<TSpace, TBoundary>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
{
    public required IMesh<TSpace> Mesh { get; init; }

    protected int[] _vertexIndices = vertexIndices;
    public ReadOnlySpan<int> Vertices => _vertexIndices;
    public abstract IEnumerable<Edge> Edges { get; }
    public abstract int EdgeCount { get; }

    public abstract bool ContainsPoint(TSpace point, double epsilon = 1e-12);
    public abstract ICoordinateTransform<TSpace, TBoundary> MasterElementCoordinateSystem { get; }
}