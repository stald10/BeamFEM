using BeamFEM.Models.FEM.Elements;
using BeamFEM.Topology;
using System.Xml.Linq;
using Telma;

namespace BeamFEM.Models.FEM.Mesh
{
    public sealed class Mesh(List<ElementTriangle> _element)
    {
        public static void ReadMesh() 
        {
            var topologyParametr = TopologyReader.ReadJSON("C:\\Users\\Lossky\\Documents\\GitHub\\BeamFEM\\Topology\\topology.json");

            int Nx = topologyParametr.Nx;
            int Ny = topologyParametr.Ny;



            Vector2D[] p = new Vector2D[Nx*Ny];
        }
    }
}
