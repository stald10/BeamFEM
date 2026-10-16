using BeamFEM.Models.FEM.Elements;
using BeamFEM.Topology;
using System.Xml.Linq;
using Telma;

namespace BeamFEM.Models.FEM.Mesh
{
    public sealed class Grid()
    {
        public static (Vector2D[] _p, Vector3D[] _cl) GridGenerator(string path) 
        {
            var topologyParametr = TopologyReader.ReadJSON(path);

            int Nx = topologyParametr.Nx;
            int Ny = topologyParametr.Ny;

            int Lx = topologyParametr.Lx;
            int Ly = topologyParametr.Ly;

            var dx = Lx / (Nx - 1);
            var dy = Ly / (Ny - 1);

            Vector2D[] p = new Vector2D[Nx*Ny];
            Vector3D[] cl = new Vector3D[2 * (Nx - 1) * (Ny - 1)];

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

                    cl[index] = new Vector3D(index1, index3, index4);
                    index++;
                    cl[index] = new Vector3D(index1, index2, index3);
                    index++;
                }
            }
            return (p,cl);
        }
    }
}
