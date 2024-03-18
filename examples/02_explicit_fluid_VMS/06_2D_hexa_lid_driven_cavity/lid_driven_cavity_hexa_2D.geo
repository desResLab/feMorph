// Reynolds number 100: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.1 Pa·s
// Reynolds number 1,000: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.01 Pa·s
// Reynolds number 10,000: lid velocity U = 1 m·s-1, L = 1 m, ρ = 10 kg·m-3, µ = 0.001 Pa·s
// OR
// Reynolds number 100: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.1 Pa·s
// Reynolds number 1,000: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.01 Pa·s
// Reynolds number 10,000: lid velocity U = 10 m·s-1, L = 1 m, ρ = 1 kg·m-3, µ = 0.001 Pa·s

SetFactory("OpenCASCADE");

// Create a cylinder
length    = 1.0; // M
height    = 1.0; // M
thickness = 0.01; // M

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
Transfinite Curve{2} = 100;
Transfinite Curve{3} = 100;
Transfinite Curve{4} = 100;
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
Physical Surface("wall_left", 6) = {5};
Physical Surface("wall_right", 7) = {3};


