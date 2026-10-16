using System;
using System.Collections.Generic;
using System.Text;
using System.Text.Json;

namespace BeamFEM.Topology
{
    public class TopologyStructer
    {
       //Ширина-Высота
       public  int Lx;
       public  int Ly;

       //Количество узлов сетки вдоль X, вдоль Y
       public  int Nx;
       public  int Ny;
    }       

    public class TopologyReader 
    {
        static public TopologyStructer ReadJSON(string fileName) => JsonSerializer.Deserialize<TopologyStructer>(File.ReadAllText(fileName))!;            
    }
}           
             