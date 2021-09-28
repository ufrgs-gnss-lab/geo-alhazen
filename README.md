Computation of specular reflections on a sphere: 
Assessment and validation based on the spherical horizon

V.H. Almeida Junior and F. Geremia-Nievinski

Functions to solve the reflection on a sphere using different algorithms.
Also, a function to solve reflection parameters on the spherical horizon. 

The main algorithm output parameters are:
- [matrix, with two columns] Reflection point position in two dimensions (x,y) (in meters)
- [vector] Grazing angle (degrees)
- [vector] Interferometric delay (in meters)
- [vector] Slant distance (in meters)
- [vector] Arc lenght (in meters)

Input parameters are:
- e [vector]: transmitting satellite elevation angle (in degrees)
- Ha [vector]: height of receiving antenna above reference surface (in meters)
- Ht [scalar, optional]: height of transmitting satellite above reference surface (in meters)
- R0 [scalar, optional]: rradius of the sphere
- algo [string, optional]: algorithm name
- traj [string, optional]: satellite trajectory

The output parameters size are dependent of number of points of the 
input parameters. On the spherical horizon, only one point is expected
as the number of point is elevation-angle dependent.
