**FEA Analysis: Two-Hole Plate**

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Calculate stresses in plate using MATLAB.

This was done in MATLAB using functions and classes from a [MATLAB package created by the ERSL](https://ersl.wisc.edu/research.html) at the University of Wisconsin-Madison.
While the results are not different from typical Finite Element Analysis methods, doing this in MATLAB allows for use of its powerful optimization tools.

The plate has two holes on the outside. One hole is where the plate is physically constrained and the othe is where the force is applied.
Three rectangular cutouts are found in the middle to reduce weight. Further optimization will be done to determine best sizing of these cutouts.

![Geometry](/MATLAB/FEA/PlatePieceFEAGeom.png)

A force of 



![Stress](/MATLAB/FEA/PlatePieceFEAStress.png)

The maority of the work in the code below is devoted to creating the nodes and edges that form the geometry. For more detail on the FEA
methods and routines, see the [ERSL](https://ersl.wisc.edu/research.html) web page.


```MATLAB
% Landon Droubay
% FEA: Two-Hole Plate

L = 6.5;        %length in meters
H = 2.5;        %height
R = 0.5;       % radius of outside holes
D = 2.0;       % distance from plate center to circles' center
W = 2;          % width of inner rectangles
h_r = 0.25;      % height of inner rectangles
d_r = 0.25;     % space between rectangles

PlatePiece.vertices = [0 H/2; 0 0; L 0; L H; 0 H; L/2-D-R H/2; L/2-D H/2+R; ...
    L/2-D+R H/2; L/2-W/2 H/2; L/2-W/2 H/2+h_r/2; L/2 H/2+h_r/2; L/2 H/2+h_r/2+d_r; ...
    L/2-W/2 H/2+h_r/2+d_r; L/2-W/2 H/2+3*h_r/2+d_r; L/2+W/2 H/2+3*h_r/2+d_r; ...
    L/2+W/2 H/2+h_r/2+d_r; L/2+W/2 H/2+h_r/2; L/2+W/2 H/2; L/2+D-R H/2; ...
    L/2+D H/2+R; L/2+D H/2-R; L/2+W/2 H/2-h_r/2; L/2 H/2-h_r/2; L/2 H/2-h_r/2-d_r; ...
    L/2+W/2 H/2-h_r/2-d_r; L/2+W/2 H/2-3*h_r/2-d_r; L/2-W/2 H/2-3*h_r/2-d_r; ...
    L/2-W/2 H/2-h_r/2-d_r; L/2-W/2 H/2-h_r/2; L/2-D H/2-R; L/2-D H/2; L/2+D H/2]';
PlatePiece.edges = [1 2; 2 3; 3 4; 4 5; 5 1; 1 6; 6 7; 7 8; 8 9; 9 10; ...
    10 11; 11 12; 12 13; 13 14; 14 15; 15 16; 16 12; 12 11; 11 17; 17 18; ...
    18 19; 19 20; 20 21; 21 19; 19 18; 18 22; 22 23; 23 24; 24 25; 25 26; ...
    26 27; 27 28; 28 24; 24 23; 23 29; 29 9; 9 8; 8 30; 30 6; 6 1]';
PlatePiece.arcs = [7 31; 8 31; 22 32; 23 32; 24 32; 38 31; 39 31]';
PlatePiece.virtual = [6; 9; 12; 18; 21; 25; 28; 34; 37; 40]';

fem = TriElasticity(PlatePiece,4000);
% fem.plotGeometry(); return;

fem = fem.setYoungsModulus(30e6);
fem = fem.setPoissonsRatio(0.28);
fem = fem.fixEdge([7 39]);
fem = fem.applyXForceOnEdge(23,100); 
fem = fem.assemble();
fem = fem.solve();
fem.plotStress();
fem.myMaxDelta;
stressed = fem.myMaxStress;
```
