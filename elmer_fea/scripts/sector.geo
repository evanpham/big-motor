// Multi-pole segment of one phase: nsec stator C-cores (pole pitch apart), the ring coil cut to the
// segment, all rotor teeth that fall inside, and air. Torque is taken on the MIDDLE pole only;
// the neighbours supply the correct surroundings, and the cut faces are flux-parallel (n x A = 0).
// (A true periodic 1-pole sector was tried; Elmer's conforming periodic BC did not converge for edge elements.)
// The Onshape parts have the motor axis along Y; they are turned on import so the motor axis is Z.
// Frame: x = radial at the middle pole, y = tangential, z = axial.
// usage: gmsh -3 -setstring core <step> -setstring coil <step> -setstring bar <step>
//             -setnumber theta <deg> [-setnumber nsec 3 -setnumber h_gap ...] -o mesh.msh sector.geo
// Physical volumes : 1 = middle core, 2 = coil, 3 = middle rotor tooth, 4 = air, 5 = other cores, 6 = other rotor teeth
// Physical surface : 1 = outer
SetFactory("OpenCASCADE");
Geometry.OCCTargetUnit = "M";

If (!Exists(theta))  theta  = 0;      EndIf  // rotor angle [deg], 0 = aligned, + = toward +y
If (!Exists(pitch))  pitch  = 9;      EndIf  // pole pitch [deg] (40 poles)
If (!Exists(nsec))   nsec   = 3;      EndIf  // number of stator poles in the model (odd)
If (!Exists(h_gap))  h_gap  = 0.5e-3; EndIf
If (!Exists(h_iron)) h_iron = 1.2e-3; EndIf
If (!Exists(h_coil)) h_coil = 2.5e-3; EndIf
If (!Exists(h_air))  h_air  = 10e-3;  EndIf
If (!Exists(r_in))   r_in   = 40e-3;  EndIf  // air extent
If (!Exists(r_out))  r_out  = 160e-3; EndIf
If (!Exists(z_half)) z_half = 80e-3;  EndIf
half = nsec*pitch/2*Pi/180;
kmax = (nsec-1)/2;

stat() = ShapeFromFile(core);
cl()   = ShapeFromFile(coil);
b0()   = ShapeFromFile(bar);
// Onshape Y (axial) -> Z, Onshape Z (tangential) -> Y
Rotate {{1,0,0}, {0,0,0}, -Pi/2} { Volume{stat(), cl(), b0()}; }
vcore = Mass Volume{stat(0)}; vbar = Mass Volume{b0(0)};

// stator copies at k*pitch
cores() = stat();
For k In {1:kmax}
  cores() += Rotate {{0,0,1}, {0,0,0},  k*pitch*Pi/180} { Duplicata{ Volume{stat()}; } };
  cores() += Rotate {{0,0,1}, {0,0,0}, -k*pitch*Pi/180} { Duplicata{ Volume{stat()}; } };
EndFor
// trim any coil material that overlaps the core (CAD clash) and report it
vc0 = 0; For i In {0:#cl()-1} vc0 += Mass Volume{cl(i)}; EndFor
cl() = BooleanDifference{ Volume{cl()}; Delete; }{ Volume{cores()}; };
vc1 = 0; For i In {0:#cl()-1} vc1 += Mass Volume{cl(i)}; EndFor
If (vc0 - vc1 > 1e-12) Printf("WARNING: coil overlapped the cores by %.3g mm3 in total, trimmed", (vc0-vc1)*1e9); EndIf

// rotor teeth at theta + k*pitch (one extra each side; the cut keeps what is inside)
Rotate {{0,0,1}, {0,0,0}, theta*Pi/180} { Volume{b0()}; }
bars() = b0();
For k In {1:kmax+1}
  bars() += Rotate {{0,0,1}, {0,0,0},  k*pitch*Pi/180} { Duplicata{ Volume{b0()}; } };
  bars() += Rotate {{0,0,1}, {0,0,0}, -k*pitch*Pi/180} { Duplicata{ Volume{b0()}; } };
EndFor

// annular segment, |angle| <= half, extruded along z
pc = newp; Point(pc) = {0,0,-z_half};
p1 = newp; Point(p1) = {r_in*Cos(-half),  r_in*Sin(-half),  -z_half};
p2 = newp; Point(p2) = {r_out*Cos(-half), r_out*Sin(-half), -z_half};
p3 = newp; Point(p3) = {r_out*Cos(half),  r_out*Sin(half),  -z_half};
p4 = newp; Point(p4) = {r_in*Cos(half),   r_in*Sin(half),   -z_half};
l1 = newl; Line(l1) = {p1,p2};
c1 = newl; Circle(c1) = {p2,pc,p3};
l2 = newl; Line(l2) = {p3,p4};
c2 = newl; Circle(c2) = {p4,pc,p1};
cl0 = newcl; Curve Loop(cl0) = {l1,c1,l2,c2};
s0 = news; Plane Surface(s0) = {cl0};
ex[] = Extrude {0,0,2*z_half} { Surface{s0}; };
sec = ex[1];

cs() = BooleanIntersection{ Volume{cl()}; Delete; }{ Volume{sec}; };
bs() = BooleanIntersection{ Volume{bars()}; Delete; }{ Volume{sec}; };
// drop rotor-tooth slivers left by the sector cut (thin wedges at the segment edge break the solve).
// min_frag = fraction of a full tooth volume; 0 (default) keeps everything, as in all earlier runs.
If (!Exists(min_frag)) min_frag = 0; EndIf
keep() = {};
For i In {0:#bs()-1}
  m = Mass Volume{bs(i)};
  Printf("rotor piece %g: %.4g of a full tooth", i, m/vbar);
  If (m >= min_frag*vbar) keep() += bs(i); Else Recursive Delete { Volume{bs(i)}; } Printf("  dropped (below min_frag)"); EndIf
EndFor
bs() = keep();
all() = BooleanFragments{ Volume{sec}; Delete; }{ Volume{cores(), cs(), bs()}; Delete; };

// classify: air = largest; cores match the core volume; rotor = centroid beyond the gap; else coil.
// "middle" = the core at angle 0 and the rotor tooth at angle theta
vmax = 0;
For i In {0:#all()-1}
  m = Mass Volume{all(i)}; If (m > vmax) vmax = m; ia = all(i); EndIf
EndFor
c1v() = {}; c5v() = {}; pk() = {}; r3v() = {}; r6v() = {}; pa() = {ia};
For i In {0:#all()-1}
  If (all(i) != ia)
    m = Mass Volume{all(i)}; c() = CenterOfMass Volume{all(i)};
    ang = Atan2(c(1), c(0))*180/Pi; r = Sqrt(c(0)^2 + c(1)^2);
    If (Fabs(m - vcore) < 1e-3*vcore)
      If (Fabs(ang) < 1) c1v() += all(i); Else c5v() += all(i); EndIf
    ElseIf (r > 100e-3)
      If (Fabs(ang - theta) < 0.5 && Fabs(m - vbar) < 1e-3*vbar) r3v() += all(i); Else r6v() += all(i); EndIf
    Else
      pk() += all(i);
    EndIf
  EndIf
EndFor
Printf("volumes: middle core %g, other cores %g, coil %g, middle tooth %g, other teeth %g, air %g", #c1v(), #c5v(), #pk(), #r3v(), #r6v(), #pa());
If (#c1v() != 1 || #c5v() != nsec-1 || #pk() < 1 || #r3v() != 1 || #pa() != 1) Error("could not classify volumes"); EndIf

// coil cross-section area (for the current density): cut with the half-plane y = 0, x > 0
rq = news; Rectangle(rq) = {0.03, -0.1, 0, 0.15, 0.2};
Rotate {{1,0,0}, {0,0,0}, Pi/2} { Surface{rq}; }
sq() = BooleanIntersection{ Surface{rq}; Delete; }{ Volume{pk()}; };
coilA = 0; For i In {0:#sq()-1} coilA += Mass Surface{sq(i)}; EndFor
Recursive Delete { Surface{sq()}; }
Printf("coil_area=%.6e", coilA);

Physical Volume("core", 1)        = c1v();
Physical Volume("coil", 2)        = pk();
Physical Volume("rotor", 3)       = r3v();
Physical Volume("air", 4)         = pa();
Physical Volume("core_side", 5)   = c5v();
If (#r6v()) Physical Volume("rotor_side", 6) = r6v(); EndIf

// exterior faces = faces that bound exactly one volume
allf() = {};
For i In {0:#all()-1}
  f() = Abs(Boundary{ Volume{all(i)}; });
  allf() += f();
EndFor
ext() = {};
For i In {0:#allf()-1}
  n = 0;
  For j In {0:#allf()-1}
    If (allf(i) == allf(j)) n++; EndIf
  EndFor
  If (n == 1) ext() += allf(i); EndIf
EndFor
Physical Surface("outer", 1) = ext();
Printf("outer faces: %g", #ext());

// ---- mesh size ----
// gap: faces at the air-gap radius (stator pole faces r=99.75, rotor inner faces r=100.25)
iron() = {c1v(), c5v(), r3v(), r6v()};
gf() = {};
fi() = Abs(Boundary{ Volume{iron()}; });
For i In {0:#fi()-1}
  c() = CenterOfMass Surface{fi(i)}; r = Sqrt(c(0)^2 + c(1)^2);
  If (r > 99.0e-3 && r < 101.0e-3) gf() += fi(i); EndIf
EndFor
Field[1] = Distance;  Field[1].SurfacesList = {gf()};
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = h_gap; Field[2].SizeMax = h_air; Field[2].DistMin = 0.8e-3; Field[2].DistMax = 15e-3;
Field[3] = Constant;  Field[3].VolumesList = {iron()}; Field[3].VIn = h_iron; Field[3].VOut = h_air;
Field[4] = Constant;  Field[4].VolumesList = {pk()};   Field[4].VIn = h_coil; Field[4].VOut = h_air;
Field[5] = Distance;  Field[5].SurfacesList = {fi()};
Field[6] = Threshold; Field[6].InField = 5;
Field[6].SizeMin = h_iron; Field[6].SizeMax = h_air; Field[6].DistMin = 1e-3; Field[6].DistMax = 40e-3;
Field[7] = Min; Field[7].FieldsList = {2, 3, 4, 6};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm3D = 10;                // HXT, multithreaded
Mesh.Optimize = 1;
Mesh.MshFileVersion = 2.2;
