**MEMBRANE**

A solver and visualization tool for finding the equilibrium configurations of diffuse interface models on two-dimensional systems interacting with both intrinsic and extrinsic geometry.

The (local) minima are attained via a generalized non-local Allen-Cahn flow, i.e. by solving a reaction-diffusion equation with non-zero chemical potential, on a discrete mesh.

**Usage**: 

Run ```make``` to compile. 

```./membrane -h``` will list all options.

Run ```./msc``` for a sample typical execution.

The 2D mesh must be in the standard **GMSH** *.msh* format with a single physical surface. See the folder ```meshes/``` for a few examples.
