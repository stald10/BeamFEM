using System;
using System.Collections.Generic;
using System.Text;

namespace BeamFEM.Models.FEM.Elements;

public static class FiniteElements
{
    public static class Segment
    {
        //public static readonly IBoundaryElementFactory2D Linear = new LinearSegmentFactory();
        //public static readonly IBoundaryElementFactory2D HierarchicalQuadratic = new HierarchicalQuadraticSegmentFactory();
        public static readonly IBoundaryElementFactory2D HierarchicalCubic = new HierarchicalCubicSegmentFactory();
        public static readonly IBoundaryElementFactory2D HierarchicalTensorCubic = new HierarchicalCubicSegmentTensorFactory();
        //public static readonly IBoundaryElementFactory2D LagrangeCubic = new LagrangeCubicSegmentFactory();
        //public static readonly IBoundaryElementFactory2D LagrangeQuadratic = new LagrangeQuadraticSegmentFactory();
        //public static readonly IBoundaryElementFactory2D Hermite = new HermiteSegmentFactory();
    }

    public static class Triangle
    {
        //public static readonly IFiniteElementFactory2D Linear = new LinearTriangleFactory();
        //public static readonly IFiniteElementFactory2D HierarchicalQuadratic = new HierarchicalQuadraticTriangleFactory();
        public static readonly IFiniteElementFactory2D HierarchicalCubic = new HierarchicalCubicTriangleFactory();
        public static readonly IFiniteElementFactory2D HierarchicalTensorCubic = new HierarchicalCubicTriangleTensorFactory();
        //public static readonly IFiniteElementFactory2D LagrangeCubic = new LagrangeCubicTriangleFactory();
    }


}
