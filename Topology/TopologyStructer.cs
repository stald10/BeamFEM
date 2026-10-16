using System;
using System.Collections.Generic;
using System.Text;
using System.Text.Json;

namespace BeamFEM.Topology
{
    public class TopologyStructer
    {
       //Ширина-Высота
       public int Lx { get; set; }
       public  int Ly { get; set; }

        //Количество узлов сетки вдоль X, вдоль Y
        public  int Nx { get; set; }
        public  int Ny { get; set; }
    }       

    public class TopologyReader 
    {
        static public TopologyStructer ReadJSON(string fileName) => JsonSerializer.Deserialize<TopologyStructer>(File.ReadAllText(fileName))!;            
    }
}           
             