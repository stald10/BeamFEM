using System;
using System.Collections.Generic;
using System.Text;
using Telma;
using BeamFEM.Models.Core.IBasis;

namespace BeamFEM.Models.Core.IElement
{
    public interface IFiniteElement
    {

        Vector2D[] points { get; }

        double[] globalDof { get; }

        IBasisSet2D basis();
    }
}
