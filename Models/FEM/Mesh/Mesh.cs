using BeamFEM.Models.Core.CoordinateSystem;
using BeamFEM.Models.FEM.Elements;
using BeamFEM.Topology;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata;
using Telma;
using Telma.Extensions;

namespace BeamFEM.Models.FEM.Mesh;

public interface IMesh<TSpace>
    where TSpace : IVectorBase<TSpace>
{
    TSpace this[int i] { get; }
    int VertexCount { get; }
    ICoordinateTransform<TSpace, TSpace> CoordinateSystem { get; }
    IEnumerable<IFiniteElement<TSpace>> FiniteElements { get; }
}


// FIXME: This is dirty workaround because TBoundary cannot be deduced from TSpace
public interface IMeshWithBoundaries<TSpace, TBoundary> : IMesh<TSpace>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
{
    IEnumerable<IBoundaryElement<TSpace, TBoundary>> BoundaryElements { get; }
}


public sealed class Mesh<TSpace, TBoundary>(ICoordinateTransform<TSpace, TSpace> coordinateSystem) :
    IMesh<TSpace>, IMeshWithBoundaries<TSpace, TBoundary>
    where TSpace : IVectorBase<TSpace>
    where TBoundary : IVectorBase<TBoundary>
{
    private readonly List<TSpace> _vertices = [];
    private readonly List<IFiniteElement<TSpace>> _finiteElements = [];
    private readonly List<IBoundaryElement<TSpace, TBoundary>> _boundaryElements = [];

    public TSpace this[int i] => _vertices[i];
    public int VertexCount => _vertices.Count;

    public IEnumerable<IFiniteElement<TSpace>> FiniteElements => _finiteElements;
    public IEnumerable<IBoundaryElement<TSpace, TBoundary>> BoundaryElements => _boundaryElements;

    public ICoordinateTransform<TSpace, TSpace> CoordinateSystem { get; } = coordinateSystem;

    public Mesh() : this(coordinateSystem: IdentityTransform<TSpace>.Instance) { }

    // Mesh construction:

    public int AddVertex(TSpace vertex)
    {
        _vertices.Add(vertex);
        return _vertices.Count - 1;
    }

    public IFiniteElement<TSpace> AddElement(IFiniteElementFactory<TSpace> factory, int[] vertices, int materialIndex)
    {
        var element = factory.CreateElement(this, vertices, materialIndex);
        _finiteElements.Add(element);
        return element;
    }

    public IBoundaryElement<TSpace, TBoundary> AddBoundary(IBoundaryElementFactory<TSpace, TBoundary> factory, int[] vertices, int boundaryIndex)
    {
        var element = factory.CreateBoundary(this, vertices, boundaryIndex);
        _boundaryElements.Add(element);
        return element;
    }
}

public sealed class Grid()
{
    public static (Vector2D[] _p, Vector3D[] _cl, List<Vector3D> _bound) GridGenerator(TopologyStructer topologyParametr) 
    {
        //var topologyParametr = TopologyReader.ReadJSON(path);

        int Nx = topologyParametr.Nx;
        int Ny = topologyParametr.Ny;

        int Lx = topologyParametr.Lx;
        int Ly = topologyParametr.Ly;

        double dx = Lx / (double)(Nx - 1);
        double dy = Ly / (double)(Ny - 1);

        Vector2D[] p = new Vector2D[Nx*Ny];
        Vector3D[] cl = new Vector3D[2 * (Nx - 1) * (Ny - 1)];
        List<Vector3D> Bound = new();

        int index = 0;
        for (int i = 0; i < Ny; ++i) 
        {
            double y = i * dy;
            for (int j = 0; j < Nx; ++j) 
            {
                double x = j * dx;
                p[index] = new Vector2D(x, y);
                index++;
            }
        }
        index = 0;
        for (int i = 0; i < Ny-1; ++i) 
        {
            for (int j = 0; j < Nx-1; ++j) 
            {
                double index1 = j + i * Nx;
                double index2 = index1 + 1;
                double index3 = index2 + Nx;
                double index4 = index1 + Nx;

                cl[index] = new Vector3D(index1, index2, index3);
                index++;
                cl[index] = new Vector3D(index1, index3, index4);
                index++;
            }
        }

        int count = 0;
        for (int i = 1; i < Nx; ++i) 
        {
            Bound.Add(new Vector3D(i - 1, i, 1));
            count++;
        }
        for (int i = Nx*Ny - Nx; i < Nx*Ny - 1; ++i) 
        {
            Bound.Add(new Vector3D(i, i + 1, 1));
            count++;
        }
        for (int i = 0; i < Ny-1; ++i) 
        {
            Bound.Add(new Vector3D(i * Nx, (i + 1) * Nx, 0));
            count++;
            
        }
        for (int i = 0; i < Ny-1; ++i) 
        {
            Bound.Add(new Vector3D((i+1) * Nx - 1, (i + 2) * Nx - 1, 1));
            count++;
        }


        return (p,cl, Bound);
    }
}
