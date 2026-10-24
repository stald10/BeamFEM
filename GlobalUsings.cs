global using IBasisFunction1D = BeamFEM.Models.FEM.Basis.IBasisFunction<Telma.Vector1D>;
global using IBasisFunction2D = BeamFEM.Models.FEM.Basis.IBasisFunction<Telma.Vector2D>;
global using IBasisFunction3D = BeamFEM.Models.FEM.Basis.IBasisFunction<Telma.Vector3D>;

global using OrientedBasisFunction1D = BeamFEM.Models.FEM.Basis.OrientedBasisFunction<Telma.Vector1D>;
global using OrientedBasisFunction2D = BeamFEM.Models.FEM.Basis.OrientedBasisFunction<Telma.Vector2D>;
global using OrientedBasisFunction3D = BeamFEM.Models.FEM.Basis.OrientedBasisFunction<Telma.Vector3D>;

global using IBasisSet1D = BeamFEM.Models.Core.IBasis.IBasisSet<Telma.Vector1D>;
global using IBasisSet2D = BeamFEM.Models.Core.IBasis.IBasisSet<Telma.Vector2D>;
global using IBasisSet3D = BeamFEM.Models.Core.IBasis.IBasisSet<Telma.Vector3D>;

global using BasisSet1D = BeamFEM.Models.Core.IBasis.BasisSet<Telma.Vector1D>;
global using BasisSet2D = BeamFEM.Models.Core.IBasis.BasisSet<Telma.Vector2D>;
global using BasisSet3D = BeamFEM.Models.Core.IBasis.BasisSet<Telma.Vector3D>;

global using IJacobyMatrix1X1 = BeamFEM.Models.Core.CoordinateSystem.IJacobyMatrix<Telma.Vector1D, Telma.Vector1D>;
global using IJacobyMatrix2X2 = BeamFEM.Models.Core.CoordinateSystem.IJacobyMatrix<Telma.Vector2D, Telma.Vector2D>;
global using IJacobyMatrix3X3 = BeamFEM.Models.Core.CoordinateSystem.IJacobyMatrix<Telma.Vector3D, Telma.Vector3D>;
global using IJacobyMatrix1X2 = BeamFEM.Models.Core.CoordinateSystem.IJacobyMatrix<Telma.Vector1D, Telma.Vector2D>;
global using IJacobyMatrix2X3 = BeamFEM.Models.Core.CoordinateSystem.IJacobyMatrix<Telma.Vector2D, Telma.Vector3D>;

global using ConstantJacobyMatrix1X1 = BeamFEM.Models.Core.CoordinateSystem.ConstantJacobyMatrix<
    Telma.Vector1D, Telma.Vector1D,    BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops1X1
>;                                     
global using ConstantJacobyMatrix2X2 = BeamFEM.Models.Core.CoordinateSystem.ConstantJacobyMatrix<
    Telma.Vector2D, Telma.Vector2D,    BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops2X2
>;                                     
global using ConstantJacobyMatrix3X3 = BeamFEM.Models.Core.CoordinateSystem.ConstantJacobyMatrix<
    Telma.Vector3D, Telma.Vector3D,    BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops3X3
>;                                    
global using ConstantJacobyMatrix1X2 = BeamFEM.Models.Core.CoordinateSystem.ConstantJacobyMatrix<
    Telma.Vector1D, Telma.Vector2D,    BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops1X2
>;                         
global using ConstantJacobyMatrix2X3 = BeamFEM.Models.Core.CoordinateSystem.ConstantJacobyMatrix<
    Telma.Vector2D, Telma.Vector3D,    BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops2X3
>;


global using IFiniteElement1D = BeamFEM.Models.FEM.Elements.IFiniteElement<Telma.Vector1D>;
global using IFiniteElement2D = BeamFEM.Models.FEM.Elements.IFiniteElement<Telma.Vector2D>;
global using IFiniteElement3D = BeamFEM.Models.FEM.Elements.IFiniteElement<Telma.Vector3D>;
global using FiniteElement1D =  BeamFEM.Models.FEM.Elements.VolumeElement<Telma.Vector1D>;
global using FiniteElement2D =  BeamFEM.Models.FEM.Elements.VolumeElement<Telma.Vector2D>;
global using FiniteElement3D =  BeamFEM.Models.FEM.Elements.VolumeElement<Telma.Vector3D>;

global using IFiniteElementFactory1D = BeamFEM.Models.FEM.Elements.IFiniteElementFactory<Telma.Vector1D>;
global using IFiniteElementFactory2D = BeamFEM.Models.FEM.Elements.IFiniteElementFactory<Telma.Vector2D>;
global using IFiniteElementFactory3D = BeamFEM.Models.FEM.Elements.IFiniteElementFactory<Telma.Vector3D>;

global using IBoundaryElement2D = BeamFEM.Models.FEM.Elements.IBoundaryElement<Telma.Vector2D, Telma.Vector1D>;
global using IBoundaryElement3D = BeamFEM.Models.FEM.Elements.IBoundaryElement<Telma.Vector3D, Telma.Vector2D>;
global using BoundaryElement2D =  BeamFEM.Models.FEM.Elements.FiniteElement<Telma.Vector2D, Telma.Vector1D>;
global using BoundaryElement3D =  BeamFEM.Models.FEM.Elements.FiniteElement<Telma.Vector3D, Telma.Vector2D>;

global using IMesh1D = BeamFEM.Models.FEM.Mesh.IMesh<Telma.Vector1D>;
global using IMesh2D = BeamFEM.Models.FEM.Mesh.IMesh<Telma.Vector2D>;
global using IMesh3D = BeamFEM.Models.FEM.Mesh.IMesh<Telma.Vector3D>;

global using Mesh2D = BeamFEM.Models.FEM.Mesh.Mesh<Telma.Vector2D, Telma.Vector1D>;
global using Mesh3D = BeamFEM.Models.FEM.Mesh.Mesh<Telma.Vector3D, Telma.Vector2D>;

global using IBoundaryElementFactory2D = BeamFEM.Models.FEM.Elements.IBoundaryElementFactory<Telma.Vector2D, Telma.Vector1D>;
global using IBoundaryElementFactory3D = BeamFEM.Models.FEM.Elements.IBoundaryElementFactory<Telma.Vector3D, Telma.Vector2D>;

global using IItegrator2D = BeamFEM.Models.FEM.Integrator.IIntegrator<Telma.Vector2D, Telma.Vector1D, BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops2X2>;
global using IItegrator3D = BeamFEM.Models.FEM.Integrator.IIntegrator<Telma.Vector3D, Telma.Vector2D, BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops3X3>;

global using NumericItegrator2D =  BeamFEM.Models.FEM.Integrator.NumericIntegrator<Telma.Vector2D, Telma.Vector1D, BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops2X2>;
global using NumbericItegrator3D = BeamFEM.Models.FEM.Integrator.NumericIntegrator<Telma.Vector3D, Telma.Vector2D, BeamFEM.Models.Core.CoordinateSystem.MatrixOperations.Ops3X3>;

global using IElementGeometry1D = BeamFEM.Models.FEM.Elements.IElementGeometry<Telma.Vector1D>;
global using IElementGeometry2D = BeamFEM.Models.FEM.Elements.IElementGeometry<Telma.Vector2D>;
global using IElementGeometry3D = BeamFEM.Models.FEM.Elements.IElementGeometry<Telma.Vector3D>;
global using VolumeElementGeometry1D = BeamFEM.Models.FEM.Elements.ElementGeometry<Telma.Vector1D, Telma.Vector1D>;
global using VolumeElementGeometry2D = BeamFEM.Models.FEM.Elements.ElementGeometry<Telma.Vector2D, Telma.Vector2D>;
global using VolumeElementGeometry3D = BeamFEM.Models.FEM.Elements.ElementGeometry<Telma.Vector3D, Telma.Vector3D>;

global using BoundaryElementGeometry2D = BeamFEM.Models.FEM.Elements.ElementGeometry<Telma.Vector2D, Telma.Vector1D>;
global using BoundaryElementGeometry3D = BeamFEM.Models.FEM.Elements.ElementGeometry<Telma.Vector3D, Telma.Vector2D>;