within Fahrzeuge;
model Fahrrad
  import MB = Modelica.Mechanics.MultiBody;
  PlanarMechanics.Parts.FixedTranslation chassisRear(r = {0, 0.5}) annotation (
    Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin={14,-60})));
  PlanarMechanics.Joints.Revolute revolute(useFlange = true,
    w(fixed=false, start=0),                                                              stateSelect = StateSelect.always,
    phi(start=0, fixed=false),
    z(start=0))                                                                                                                                            annotation (
    Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {14, 2})));
  inner PlanarMechanics.PlanarWorld planarWorld(defaultWidthFraction = 10, defaultZPosition = 0, constantGravity = {0, 0}) annotation (
    Placement(transformation(extent = {{-100, 80}, {-80, 100}})));
  PlanarMechanics.Parts.FixedTranslation chassisFront(r = {0, 0.5}) annotation (
    Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {14, -28})));
  PlanarMechanics.Parts.Body bodyCenter(I = 1.8, m = 1450, phi(fixed = false), w(fixed = false), v(fixed = true), r(fixed = true), enableGravity = false) annotation (
    Placement(transformation(extent = {{60, -54}, {80, -34}})));
  Modelica.Mechanics.Rotational.Sources.Position position(useSupport = false, exact = false, f_crit = 50) annotation (
    Placement(visible = true, transformation(extent = {{90, 0}, {70, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.Trapezoid trapezoid(amplitude = 0 * Modelica.Constants.pi / 4, rising = 0.1, width = 0.8, falling = 0.1, period = 1.8,
    offset=1*Modelica.Constants.pi/4)                                                                                                            annotation (
    Placement(transformation(extent={{42,60},{62,80}})));
  Rad rad2(radius = 0.3, r = {0, 5}, N = 0, vAdhesion_min = 0.05, vSlide_min = 0.15, sAdhesion = 0.04, sSlide = 0.12, mu_A = 0.8, mu_S = 0.4, B = 7, C = 1.6, D = 1, grenze = 1e-2, w_roll(start = 0)) annotation (
    Placement(transformation(extent = {{10, -70}, {-10, -90}})));
  //  h=0.4,
  // l_F=1.1,
  //  l_R=1.59,
  Fahrzeuge.Rad rad(radius = 0.3, r = {0, 5}, N = 0, vAdhesion_min = 0.05, vSlide_min = 0.15, sAdhesion = 0.04, sSlide = 0.12, mu_A = 0.8, mu_S = 0.4, B = 7, C = 1.6, D = 1, grenze = 1e-2, phi_roll(start = 0), w_roll(start = 1, fixed = true)) annotation (
    Placement(visible = true, transformation(extent={{10,60},{-10,80}},      rotation = 0)));
  //  h=0.4,
  //  l_F=1.1,
  //  l_R=1.59,
  Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 1.8) annotation (
    Placement(transformation(extent = {{-60, 60}, {-40, 80}})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia1(J = 1.8) annotation (
    Placement(transformation(extent = {{-50, -90}, {-30, -70}})));
  Fahrzeuge.Lastverteilung lastverteilung(g = 9.81, h = 0.4, l_F = 1.1, l_R = 1.59, m = bodyCenter.m) annotation (
    Placement(visible = true, transformation(extent={{-70,10},{-50,30}},      rotation = 0)));
  Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor annotation (
    Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 0, origin = {50, -20})));
  Modelica.Mechanics.Rotational.Sources.Torque torque
    annotation (Placement(transformation(extent={{-84,-88},{-64,-68}})));
  Modelica.Blocks.Sources.Constant const(k=100)
    annotation (Placement(transformation(extent={{-166,-74},{-146,-54}})));
equation
  connect(chassisFront.frame_b, revolute.frame_a) annotation (
    Line(points = {{14, -18}, {14, -8}}, color = {95, 95, 95}, thickness = 0.5));
  connect(chassisRear.frame_b, chassisFront.frame_a) annotation (
    Line(points={{14,-50},{14,-38}},      color = {95, 95, 95}, thickness = 0.5));
  connect(bodyCenter.frame_a, chassisRear.frame_b) annotation (
    Line(points={{60,-44},{38,-44},{38,-50},{14,-50}},
                                          color = {95, 95, 95}, thickness = 0.5));
  connect(position.flange, revolute.flange_a) annotation (
    Line(points = {{70, 10}, {30, 10}, {30, 2}, {24, 2}}));
  connect(rad2.frame_a, chassisRear.frame_a) annotation (
    Line(points={{4,-80},{14,-80},{14,-70}},        color = {95, 95, 95}, thickness = 0.5));
  connect(rad.frame_a, revolute.frame_b) annotation (
    Line(points={{4,70},{14,70},{14,12}},        color = {95, 95, 95}, thickness = 0.5));
  connect(inertia.flange_b, rad.flange_a) annotation (
    Line(points={{-40,70},{-10,70}}));
  connect(inertia1.flange_b, rad2.flange_a) annotation (
    Line(points = {{-30, -80}, {-10, -80}}, color = {0, 0, 0}));
  connect(lastverteilung.f_F, rad.dynamicLoad) annotation (
    Line(points={{-49,24},{0,24},{0,60}},      color = {0, 0, 127}));
  connect(lastverteilung.f_R, rad2.dynamicLoad) annotation (
    Line(points={{-49,17},{-2,17},{-2,-64},{0,-64},{0,-70}},
                                                  color = {0, 0, 127}));
  connect(angleSensor.flange, position.flange) annotation (
    Line(points = {{60, -20}, {60, 10}, {70, 10}}, color = {0, 0, 0}));
  connect(angleSensor.phi, lastverteilung.delta) annotation (
    Line(points={{39,-20},{34,-20},{34,16},{0,16},{0,14},{-40,14},{-40,2},{-76,
          2},{-76,20},{-70,20}},                                 color = {0, 0, 127}));
  connect(position.phi_ref, trapezoid.y) annotation (Line(points={{92,10},{98,
          10},{98,70},{63,70}}, color={0,0,127}));
  connect(inertia1.flange_a, torque.flange)
    annotation (Line(points={{-50,-80},{-50,-78},{-64,-78}}, color={0,0,0}));
  connect(lastverteilung.f_lo, rad.f_long_out) annotation (Line(points={{-56,30},
          {-56,56},{-64,56},{-64,88},{0,88},{0,81}},
                                   color={0,0,127}));
  connect(lastverteilung.f_la, rad.f_lat_out) annotation (Line(points={{-64,30},
          {-64,38},{-28,38},{-28,84},{-5,84},{-5,78}}, color={0,0,127}));
  connect(lastverteilung.f_R_long, rad2.f_long_out) annotation (Line(points={{-65.4,
          12},{-64,12},{-64,-66},{-20,-66},{-20,-96},{0,-96},{0,-91}},
        color={0,0,127}));
  connect(const.y, torque.tau) annotation (Line(points={{-145,-64},{-94,-64},{
          -94,-78},{-86,-78}}, color={0,0,127}));
  connect(lastverteilung.M, const.y) annotation (Line(points={{-56,9.2},{-56,0},
          {-120,0},{-120,-64},{-145,-64}},    color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio = false)),
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    uses(Modelica(version = "3.2.3")));
end Fahrrad;
