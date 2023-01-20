within ;
package Fahrzeuge "Modell"
  extends Modelica.Icons.Package;
  import Modelica.Units.SI;
  import MB = Modelica.Mechanics.MultiBody;



  annotation (
    version="0.0.3",
    versionDate="2021-11-30",
    dateModified = "2021-11-30 12:00:00Z",
    uses(                                  Modelica(version="4.0.0"),
      PlanarMechanics(version="1.5.1")),
  conversion(noneFromVersion="0.0.1", noneFromVersion="0.0.2"));
end Fahrzeuge;
