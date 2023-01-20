within Fahrzeuge;
model Lastverteilung
  import      Modelica.Units.SI;

  parameter SI.Mass m;
  parameter SI.Acceleration g = 9.81;

  parameter SI.Length h;
  parameter SI.Length l_R;
  parameter SI.Length l_F;

  Modelica.Blocks.Interfaces.RealInput delta annotation (
    Placement(visible = true, transformation(origin={-80,-52},   extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin={-100,0},     extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput f_lo annotation (Placement(
      visible=true,
      transformation(
        origin={-80,80},
        extent={{-20,-20},{20,20}},
        rotation=0),
      iconTransformation(
        origin={40,100},
        extent={{-20,-20},{20,20}},
        rotation=270)));
  Modelica.Blocks.Interfaces.RealInput f_la annotation (Placement(
      visible=true,
      transformation(
        origin={-80,50},
        extent={{-20,-20},{20,20}},
        rotation=0),
      iconTransformation(
        origin={-40,100},
        extent={{-20,-20},{20,20}},
        rotation=-90)));
  Modelica.Blocks.Interfaces.RealInput f_R_long annotation (Placement(
      visible=true,
      transformation(
        origin={-80,-86},
        extent={{-20,-20},{20,20}},
        rotation=0),
      iconTransformation(
        origin={-54,-80},
        extent={{-20,-20},{20,20}},
        rotation=90)));


  Modelica.Blocks.Interfaces.RealOutput f_F
    annotation (Placement(transformation(extent={{100,30},{120,50}}),
        iconTransformation(extent={{100,30},{120,50}})));

  Modelica.Blocks.Interfaces.RealOutput f_R
    annotation (Placement(transformation(extent={{100,-40},{120,-20}}),
        iconTransformation(extent={{100,-40},{120,-20}})));


  Modelica.Blocks.Interfaces.RealInput M annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,-100}), iconTransformation(extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={40,-108})));
equation

  f_R =( h*(f_lo * cos(delta)- f_la * sin(delta)+ f_R_long)+f_F*l_F)/l_R;

  M = f_F+f_R-m*g;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")));
end Lastverteilung;
