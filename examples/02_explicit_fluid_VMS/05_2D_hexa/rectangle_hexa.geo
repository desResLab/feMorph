SetFactory("OpenCASCADE");

// Create a cylinder
length    = 30.0; // CM
height    = 4.0; // CM
thickness = 0.2; // CM

// Define 4 points
Point(1) = {0,0,0,0.0};
Point(2) = {length,0,0,0.0};
Point(3) = {length,height,0,0.0};
Point(4) = {0,height,0,0.0};

// Define the 4 circles
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(9) = {1,2,3,4};

// Create solid inlet surface
Plane Surface(1) = {9};

// Set mesh size
Transfinite Curve{1} = 100;
Transfinite Curve{2} = 20;
Transfinite Curve{3} = 100;
Transfinite Curve{4} = 20;
Transfinite Surface{1} = {1, 2, 3, 4};
Recombine Surface{1};

// Name all the surfaces and planes
Extrude {0,0,thickness} {Surface{1}; Layers{1};Recombine;}

// Name all the surfaces and planes
Physical Volume("fluid", 1) = {1};
Physical Surface("z_wall_1", 2) = {1};
Physical Surface("z_wall_2", 3) = {6};
Physical Surface("wall_top", 4) = {4};
Physical Surface("wall_bottom", 5) = {2};
Physical Surface("inlet", 6) = {5};
Physical Surface("outlet", 7) = {3};


