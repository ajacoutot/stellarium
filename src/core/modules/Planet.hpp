/*
 * Stellarium
 * Copyright (C) 2002 Fabien Chereau
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */

#ifndef _PLANET_HPP_
#define _PLANET_HPP_

#include "StelObject.hpp"
#include "StelProjector.hpp"
#include "VecMath.hpp"
#include "StelFader.hpp"
#include "StelTextureTypes.hpp"
#include "StelProjectorType.hpp"

#include <QString>

// The callback type for the external position computation function
// The last variable is the userData pointer, which is NULL for Planets, but used in derived classes. E.g. points to the CometOrbit for Comets.
typedef void (*posFuncType)(double, double*, void*);

// GZ2016: new axis functions for external computation of axis orientations for selected objects.
// The last variable is a pointer to the Planet object.
typedef void (*axisFuncType)(double, double*, void*);


typedef void (OsculatingFunctType)(double jde0,double jde,double xyz[3]);

// epoch J2000: 12 UT on 1 Jan 2000
#define J2000 2451545.0
#define ORBIT_SEGMENTS 360

class StelFont;
class StelPainter;
class StelTranslator;
class QOpenGLShaderProgram;

// Class used to store rotational elements, i.e. axis orientation for the planetary body.
// Data are read from ssystem.ini in SolarSystem.cpp::loadPlanets().
//
// GZ2016-01: This seems to be an outdated model.
// TODO: For historical reasons, some of the values are given in hours, others in earth days, making life harder. Most moons have coupled rotation, they should be auto-converted from
//       Finding a consistent scheme would be helpful. Also, IAU has planet axes in ICRF, while Stellarium has them coded in VSOP87.
//       Except if rot_pole... values are given, then they are ICRF and transformed on the fly to VSOP87, stored in here.
//
// IAU standards like the Report of the IAU Working Group on Cartographic Coordinates and Rotational elements 2009 (DOI:10.1007/s10569-010-9320-4)
// has axes given w.r.t. J2000 ICRF.
// Before 0.15, if the planet elements in ssystem.ini had elements "rot_pole_ra" and "rot_pole_de" given, the poles were transformed to
// ecliptically-based directions (VSOP87) obliquity/ascendingNode. But these were only ra/de_J2000 poles. Some axes have precession, which need to be modelled/updated.
// The new way (0.15): We still use the previous elements obliquity and ascendingNode in Planet::computeTransMatrix(.,.)
// and Planet::getSiderealTime(), but can update them periodically from:
// ra=ra0+T*ra1
// de=de0+T*de1         ( --> obliquity, ascendingNode)
// rot_rotation_offset [degrees]     =W0
// rot_periode  [hours] =  computed from rot_pole_W1[deg/day] if that exists.  360 [deg] / _rot_ [hours] --> 360 * _rot_ / 24 [deg/hours]
// In addition, the objects with more complicated element behaviour must be updated with special functions in those two functions.
// New keys in ssystem.ini, their storage in RotationElements, and their equivalents in the IAU report:
// rot_pole_ra  [degrees]     re.ra0      constant term for alpha_0
// rot_pole_de  [degrees]     re.ra1      constant term for delta_0
// rot_pole_ra1 [degrees/cy]  re.de0      T factor for alpha_0
// rot_pole_de1 [degrees/cy]  re.de1      T factor for  delta_0
// rot_pole_W0  [degrees]     re.offset   constant term fo W. Will be stored in re.offset, and is signal for using the IAU model (causes setting re.useICRF=true)
// rot_pole_W1  [degrees/day] re.period   d factor for W      Will be converted to rot_periode=360*24/W1.

class RotationElements
{
public:
	RotationElements(void) : period(1.), offset(0.), epoch(J2000), obliquity(0.), ascendingNode(0.), //precessionRate(0.),
		siderealPeriod(0.),
		useICRF(false), ra0(0.), ra1(0.), de0(0.), de1(0.), W0(0.), W1(0.) {}
	float period;          // (sidereal) rotation period [earth days]   CURRENTLY NOT:  If useICRF, this is from the time term of W.
	float offset;          // rotation at epoch  [degrees]              CURRENTLY NOT:  If useICRF, this is the constant term of W
	double epoch;          // JDE (JD TT) of epoch for these elements
	float obliquity;       // tilt of rotation axis w.r.t. ecliptic [radians]
	float ascendingNode;   // long. of ascending node of equator on the ecliptic [radians]
	// Field rot_precession_rate in ssystem.ini is no longer used. We still keep Earth's value as it is evaluated in older versions (until 0.13.*).
//	float precessionRate;  // rate of precession of rotation axis in [rads/JulianCentury(36525d)] [ NO LONGER USED WITH 0.14 (was used for Earth only, and that was too simple.) ]
	double siderealPeriod; // sidereal period (Planet year in earth days) [earth days]
	// GZ for 0.16: I propose changes here: The 6 new entries after the switch are enough for many objects. Else, design special_functions.
	bool useICRF;          // Use values w.r.t. ICRF (should ultimately be true for all objects!) This can be set when rot_pole_w0 is given. Updating the axis is required if ra1<>0
	double ra0;            // [rad] RA_0 right ascension of north pole. ssystem.ini: rot_pole_ra    /180*M_PI
	double ra1;            // [rad/century] rate of change in axis ra   ssystem.ini: rot_pole_ra1   /180*M_PI
	double de0;            // [rad] DE_0 declination of north pole      ssystem.ini: rot_pole_de    /180*M_PI
	double de1;            // [rad/century] rate of change in axis de   ssystem.ini: rot_pole_de1   /180*M_PI
	// These values are only in the modern algorithms. invalid if W0=0.
	double W0;             // [deg] mean longitude of prime meridian along equator measured from intersection with ICRS plane at epoch.
	double W1;             // [deg/d] mean longitude motion. W=W0+d*W1.
};

// Class to manage rings for planets like saturn
class Ring
{
public:
	Ring(float radiusMin, float radiusMax,const QString &texname);
	double getSize(void) const {return radiusMax;}
	const float radiusMin;
	const float radiusMax;
	StelTextureSP tex;
};


class Planet : public StelObject
{
public:
	friend class SolarSystem;

	Q_ENUMS(PlanetType)
	Q_ENUMS(PlanetOrbitColorStyle)
	Q_ENUMS(ApparentMagnitudeAlgorithm)
	//! numeric typecodes for the type descriptions in ssystem.ini
	// GZ: Until 0.13 QStrings were used for types.
	// GZ: Enums are slightly faster than string comparisons in time-critical comparisons.
	// GZ: If other types are introduced, add here and the string in init().
	enum PlanetType
	{
		isStar,			// ssystem.ini: type="star"
		isPlanet,		// ssystem.ini: type="planet"
		isMoon,			// ssystem.ini: type="moon"
		isAsteroid,		// ssystem.ini: type="asteroid"
		isPlutino,		// ssystem.ini: type="plutino"
		isComet,		// ssystem.ini: type="comet"
		isDwarfPlanet,		// ssystem.ini: type="dwarf planet"
		isCubewano,		// ssystem.ini: type="cubewano"
		isSDO,			// ssystem.ini: type="sdo"
		isOCO,			// ssystem.ini: type="oco"
		isUNDEFINED		// ssystem.ini: type=<anything else>
	};

	enum PlanetOrbitColorStyle
	{
		ocsOneColor,		// One color for all orbits
		ocsGroups,		// Separate colors for each group of Solar system bodies
		ocsMajorPlanets		// Separate colors for each of major planets of Solar system
	};

	enum ApparentMagnitudeAlgorithm
	{
		Mueller_1893,	// G. Mueller, based on visual observations 1877-91. [Expl.Suppl.1961]
		Astr_Alm_1984,	// Astronomical Almanac 1984 and later. These give V (instrumental) magnitudes (allegedly from D.L. Harris, but this is wrong!)
		Expl_Sup_1992,	// Algorithm provided by Pere Planesas (Observatorio Astronomico Nacional) (Was called "Planesas")
		Expl_Sup_2013,	// Explanatory Supplement to the Astronomical Almanac, 3rd edition 2013
//		Planesas,		// Algorithm provided by Pere Planesas (Observatorio Astronomico Nacional)
//		Mueller,		// G. Mueller, based on visual observations 1877-91. [Expl.Suppl.1961]
//		Harris,			// Astronomical Almanac 1984 and later. These give V (instrumental) magnitudes (D.L. Harris)
		UndefinedAlgorithm,
		Generic		// Visual magnitude based on phase angle and albedo. The formula source for this is totally unknown!
	};

	Planet(const QString& englishName,
	       int flagLighting,
	       double radius,
	       double oblateness,
	       Vec3f halocolor,
	       float albedo,
	       const QString& texMapName,
	       const QString& normalMapName,
	       posFuncType _coordFunc,
	       void* userDataPtr,
	       OsculatingFunctType *osculatingFunc,
	       bool closeOrbit,
	       bool hidden,
	       bool hasAtmosphere,
	       bool hasHalo,
	       const QString &pTypeStr);

	virtual ~Planet();

	//! Initializes static vars. Must be called before creating first planet.
	// Currently ensured by SolarSystem::init()
	static void init();

	///////////////////////////////////////////////////////////////////////////
	// Methods inherited from StelObject
	//! Get a string with data about the Planet.
	//! Planets support the following InfoStringGroup flags:
	//! - Name
	//! - Magnitude
	//! - RaDec
	//! - AltAzi
	//! - Distance
	//! - Size
	//! - PlainText
	//! - Extra: Heliocentric Ecliptical Coordinates & Observer-planetocentric Ecliptical Coordinates, Phase, illumination, phase angle & elongation from the Sun
	//! @param core the StelCore object
	//! @param flags a set of InfoStringGroup items to include in the return value.
	//! @return a QString containing an HMTL encoded description of the Planet.
	virtual QString getInfoString(const StelCore *core, const InfoStringGroup& flags) const;
	virtual double getCloseViewFov(const StelCore* core) const;
	virtual double getSatellitesFov(const StelCore* core) const;
	virtual double getParentSatellitesFov(const StelCore* core) const;
	virtual float getVMagnitude(const StelCore* core) const;
	virtual float getSelectPriority(const StelCore* core) const;
	virtual Vec3f getInfoColor(void) const;
	virtual QString getType(void) const {return "Planet";}
	virtual Vec3d getJ2000EquatorialPos(const StelCore *core) const;
	virtual QString getEnglishName(void) const;
	virtual QString getNameI18n(void) const;
	//! Get angular semidiameter, degrees. If planet display is artificially enlarged (e.g. Moon upscale), value will also be increased.
	virtual double getAngularSize(const StelCore* core) const;
	virtual bool hasAtmosphere(void) {return atmosphere;}
	virtual bool hasHalo(void) {return halo;}

	///////////////////////////////////////////////////////////////////////////
	// Methods of SolarSystem object
	//! Translate planet name using the passed translator
	virtual void translateName(const StelTranslator &trans);

	// Draw the Planet
	// GZ Made that virtual to allow comets having their own draw().
	virtual void draw(StelCore* core, float maxMagLabels, const QFont& planetNameFont);

	///////////////////////////////////////////////////////////////////////////
	// Methods specific to Planet
	//! Get the radius of the planet in AU.
	//! @return the radius of the planet in astronomical units.
	double getRadius(void) const {return radius;}
	//! Get the value (1-f) for oblateness f.
	double getOneMinusOblateness(void) const {return oneMinusOblateness;}
	//! Get duration of sidereal day (earth days, may come from rot_periode or orbit_period (for moons) from ssystem.ini)
	double getSiderealDay(void) const { if (re.W1) return 360.0/re.W1; else return re.period;} // I assume the more modern values are better.
	//! Get duration of sidereal year
	// must be virtual for Comets.
	virtual double getSiderealPeriod(void) const { return re.siderealPeriod; }
	//! Get duration of mean solar day, in earth days.
	double getMeanSolarDay(void) const;

	const QString& getTextMapName() const {return texMapName;}	
	const QString getPlanetTypeString() const {return pTypeMap.value(pType);}
	PlanetType getPlanetType() const {return pType;}

	void setNativeName(QString planet) { nativeName = planet; }

	//! Return the absolute magnitude (read from file ssystem.ini)
	float getAbsoluteMagnitude() const {return absoluteMagnitude;}
	//! Return the mean opposition magnitude, defined as V(1,0)+5log10(a(a-1))
	//! A return value of 100 signals invalid result.
	float getMeanOppositionMagnitude() const;
	ApparentMagnitudeAlgorithm getApparentMagnitudeAlgorithm() const { return vMagAlgorithm; }
	const QString getApparentMagnitudeAlgorithmString() const { return vMagAlgorithmMap.value(vMagAlgorithm); }
	void setApparentMagnitudeAlgorithm(QString algorithm);

	//! Compute the z rotation to use from equatorial to geographic coordinates. For general applicability we need both time flavours:
	//! @param JD is JD(UT) for Earth
	//! @param JDE is used for other locations
	double getSiderealTime(double JD, double JDE) const;
	Mat4d getRotEquatorialToVsop87(void) const;
	void setRotEquatorialToVsop87(const Mat4d &m);

	const RotationElements &getRotationElements(void) const {return re;}
	// Set the -o-r-b-i-t-a-l- ROTATIONAL elements
	// _period: duration of sidereal rotation [Julian days]
	// _offset: [angle at _epoch. ]
	// _epoch: [JDE]
	// _obliquity [rad]
	// _ascendingNode of equator on ecliptic[rad]
	// ra_pole=_ra0 + T*_ra1. ra_pole and de_pole must be computed more than for initialisation for J2000
	// de_pole=_de0 + T*_de1. ra and de values to be stored in [rad]
	// _w0, _w1 to be given in degrees!
	// _precessionRate [rad/JulCt] (was only given for earth, and is no longer used!)
	// _siderealPeriod [earth days] orbital duration. THIS DOES NOT BELONG HERE!
	void setRotationElements(const float _period, const float _offset, const double _epoch,
				 const float _obliquity, const float _ascendingNode,
				 const double _ra0,
				 const double _ra1,
				 const double _de0,
				 const double _de1,
				 const double _w0,
				 const double _w1,
				 //float _precessionRate,
				 const double _siderealPeriod);
	double getRotAscendingnode(void) const {return re.ascendingNode;}
	// return angle between axis and normal of ecliptic plane (or, for a moon, equatorial/reference plane defined by parent).
	// For Earth, this is the angle between axis and normal to current ecliptic of date, i.e. the ecliptic obliquity of date JDE.
	// TODO: decide if this is always angle between axis and J2000 ecliptic, or should be axis//current ecliptic!
	double getRotObliquity(double JDE) const;


	// Compute the position in the parent Planet coordinate system
	void computePositionWithoutOrbits(const double dateJDE);
	void computePosition(const double dateJDE);

	// Compute the transformation matrix from the local Planet coordinate to the parent Planet coordinate.
	// This requires both flavours of JD in cases involving Earth.
	void computeTransMatrix(double JD, double JDE);

	//! Get the phase angle (radians) for an observer at pos obsPos in heliocentric coordinates (in AU)
	double getPhaseAngle(const Vec3d& obsPos) const;
	//! Get the elongation angle (radians) for an observer at pos obsPos in heliocentric coordinates (in AU)
	double getElongation(const Vec3d& obsPos) const;
	//! Get the angular radius (degrees) of the planet spheroid (i.e. without the rings)
	double getSpheroidAngularSize(const StelCore* core) const;
	//! Get the planet phase (illuminated fraction of the planet disk, 0..1) for an observer at pos obsPos in heliocentric coordinates (in AU)
	float getPhase(const Vec3d& obsPos) const;

	//! Get the Planet position in the parent Planet ecliptic coordinate in AU
	Vec3d getEclipticPos() const;

	// Return the heliocentric ecliptical position
	Vec3d getHeliocentricEclipticPos() const {return getHeliocentricPos(eclipticPos);}

	//! Return the heliocentric transformation for local (parentocentric) coordinate
	//! @arg p planetocentric rectangular ecliptical coordinate (J2000)
	//! @return heliocentric rectangular ecliptical coordinates (J2000)
	Vec3d getHeliocentricPos(Vec3d p) const;
	//! Propagate the heliocentric coordinates to parentocentric coordinates
	//! @arg pos heliocentric rectangular ecliptical coordinate (J2000)
	void setHeliocentricEclipticPos(const Vec3d &pos);

	//! Compute and return the distance to the given position in heliocentric ecliptical (J2000) coordinate (in AU)
	double computeDistance(const Vec3d& obsHelioPos);
	//! Return the last computed distance to the given position in heliocentric ecliptical (J2000) coordinate (in AU)
	double getDistance(void) const {return distance;}

	void setRings(Ring* r) {rings = r;}

	void setSphereScale(float s) {sphereScale = s;}
	float getSphereScale(void) const {return sphereScale;}

	const QSharedPointer<Planet> getParent(void) const {return parent;}

	static void setLabelColor(const Vec3f& lc) {labelColor = lc;}
	static const Vec3f& getLabelColor(void) {return labelColor;}

	// update displayed elements. @param deltaTime: ms (since last call)
	virtual void update(int deltaTime);

	void setFlagHints(bool b){hintFader = b;}
	bool getFlagHints(void) const {return hintFader;}

	void setFlagLabels(bool b){flagLabels = b;}
	bool getFlagLabels(void) const {return flagLabels;}

	bool flagNativeName;
	void setFlagNativeName(bool b) { flagNativeName = b; }
	bool getFlagNativeName(void) { return flagNativeName; }

	bool flagTranslatedName;
	void setFlagTranslatedName(bool b) { flagTranslatedName = b; }
	bool getFlagTranslatedName(void) { return flagTranslatedName; }

	///////////////////////////////////////////////////////////////////////////
	// DEPRECATED
	///// Orbit related code
	// Should move to an OrbitPath class which works on a SolarSystemObject, not a Planet
	void setFlagOrbits(bool b){orbitFader = b;}
	bool getFlagOrbits(void) const {return orbitFader;}
	LinearFader orbitFader;
	// draw orbital path of Planet
	void drawOrbit(const StelCore*);
	Vec3d orbit[ORBIT_SEGMENTS+1];   // store heliocentric coordinates for drawing the orbit
	Vec3d orbitP[ORBIT_SEGMENTS+1];  // store local coordinate for orbit
	double lastOrbitJDE;
	double deltaJDE;                 // time difference between positional updates.
	double deltaOrbitJDE;
	bool orbitCached;                // whether orbit calculations are cached for drawing orbit yet
	bool closeOrbit;                 // whether to connect the beginning of the orbit line to
					 // the end: good for elliptical orbits, bad for parabolic
					 // and hyperbolic orbits

	static Vec3f orbitColor;
	static void setOrbitColor(const Vec3f& oc) {orbitColor = oc;}
	static const Vec3f& getOrbitColor() {return orbitColor;}

	static Vec3f orbitMajorPlanetsColor;
	static void setMajorPlanetOrbitColor(const Vec3f& oc) { orbitMajorPlanetsColor = oc;}
	static const Vec3f& getMajorPlanetOrbitColor() {return orbitMajorPlanetsColor;}

	static Vec3f orbitMoonsColor;
	static void setMoonOrbitColor(const Vec3f& oc) { orbitMoonsColor = oc;}
	static const Vec3f& getMoonOrbitColor() {return orbitMoonsColor;}

	static Vec3f orbitMinorPlanetsColor;
	static void setMinorPlanetOrbitColor(const Vec3f& oc) { orbitMinorPlanetsColor = oc;}
	static const Vec3f& getMinorPlanetOrbitColor() {return orbitMinorPlanetsColor;}

	static Vec3f orbitDwarfPlanetsColor;
	static void setDwarfPlanetOrbitColor(const Vec3f& oc) { orbitDwarfPlanetsColor = oc;}
	static const Vec3f& getDwarfPlanetOrbitColor() {return orbitDwarfPlanetsColor;}

	static Vec3f orbitCubewanosColor;
	static void setCubewanoOrbitColor(const Vec3f& oc) { orbitCubewanosColor = oc;}
	static const Vec3f& getCubewanoOrbitColor() {return orbitCubewanosColor;}

	static Vec3f orbitPlutinosColor;
	static void setPlutinoOrbitColor(const Vec3f& oc) { orbitPlutinosColor = oc;}
	static const Vec3f& getPlutinoOrbitColor() {return orbitPlutinosColor;}

	static Vec3f orbitScatteredDiscObjectsColor;
	static void setScatteredDiscObjectOrbitColor(const Vec3f& oc) { orbitScatteredDiscObjectsColor = oc;}
	static const Vec3f& getScatteredDiscObjectOrbitColor() {return orbitScatteredDiscObjectsColor;}

	static Vec3f orbitOortCloudObjectsColor;
	static void setOortCloudObjectOrbitColor(const Vec3f& oc) { orbitOortCloudObjectsColor = oc;}
	static const Vec3f& getOortCloudObjectOrbitColor() {return orbitOortCloudObjectsColor;}

	static Vec3f orbitCometsColor;
	static void setCometOrbitColor(const Vec3f& oc) { orbitCometsColor = oc;}
	static const Vec3f& getCometOrbitColor() {return orbitCometsColor;}

	static Vec3f orbitMercuryColor;
	static void setMercuryOrbitColor(const Vec3f& oc) { orbitMercuryColor = oc;}
	static const Vec3f& getMercuryOrbitColor() {return orbitMercuryColor;}

	static Vec3f orbitVenusColor;
	static void setVenusOrbitColor(const Vec3f& oc) { orbitVenusColor = oc;}
	static const Vec3f& getVenusOrbitColor() {return orbitVenusColor;}

	static Vec3f orbitEarthColor;
	static void setEarthOrbitColor(const Vec3f& oc) { orbitEarthColor = oc;}
	static const Vec3f& getEarthOrbitColor() {return orbitEarthColor;}

	static Vec3f orbitMarsColor;
	static void setMarsOrbitColor(const Vec3f& oc) { orbitMarsColor = oc;}
	static const Vec3f& getMarsOrbitColor() {return orbitMarsColor;}

	static Vec3f orbitJupiterColor;
	static void setJupiterOrbitColor(const Vec3f& oc) { orbitJupiterColor = oc;}
	static const Vec3f& getJupiterOrbitColor() {return orbitJupiterColor;}

	static Vec3f orbitSaturnColor;
	static void setSaturnOrbitColor(const Vec3f& oc) { orbitSaturnColor = oc;}
	static const Vec3f& getSaturnOrbitColor() {return orbitSaturnColor;}

	static Vec3f orbitUranusColor;
	static void setUranusOrbitColor(const Vec3f& oc) { orbitUranusColor = oc;}
	static const Vec3f& getUranusOrbitColor() {return orbitUranusColor;}

	static Vec3f orbitNeptuneColor;
	static void setNeptuneOrbitColor(const Vec3f& oc) { orbitNeptuneColor = oc;}
	static const Vec3f& getNeptuneOrbitColor() {return orbitNeptuneColor;}

	static bool permanentDrawingOrbits;
	static PlanetOrbitColorStyle orbitColorStyle;

	//! Return the list of planets which project some shadow on this planet
	QVector<const Planet*> getCandidatesForShadow() const;
	
protected:
	static StelTextureSP texEarthShadow;     // for lunar eclipses

	void computeModelMatrix(Mat4d &result) const;

	Vec3f getCurrentOrbitColor();
	
	// Return the information string "ready to print" :)
	QString getSkyLabel(const StelCore* core) const;

	// Draw the 3d model. Call the proper functions if there are rings etc..
	void draw3dModel(StelCore* core, StelProjector::ModelViewTranformP transfo, float screenSz, bool drawOnlyRing=false);

	// Draw the 3D sphere
	void drawSphere(StelPainter* painter, float screenSz, bool drawOnlyRing=false);

	// Draw the circle and name of the Planet
	void drawHints(const StelCore* core, const QFont& planetNameFont);

	QString englishName;             // english planet name
	QString nameI18;                 // International translated name
	QString nativeName;              // Can be used in a skyculture
	QString texMapName;              // Texture file path
	QString normalMapName;              // Texture file path
	int flagLighting;                // Set whether light computation has to be proceed
	RotationElements re;             // Rotation param
	double radius;                   // Planet equatorial radius in AU
	double oneMinusOblateness;       // (polar radius)/(equatorial radius)
	Vec3d eclipticPos;               // Position in AU in the rectangular ecliptic coordinate system, (GZ2016: presumably equinox J2000)
					 // centered on the parent Planet
	Vec3d screenPos;                 // Used to store temporarily the 2D position on screen
	Vec3d previousScreenPos;         // The position of this planet in the previous frame.
	Vec3f haloColor;                 // exclusively used for drawing the planet halo

	float albedo;                    // Planet albedo. Used for magnitude computation (but formula dubious!)
	float absoluteMagnitude;         // since 2017: V(1,0) from Explanatory Supplement or WGCCRE2009 paper for the planets, H in the H,G magnitude system for Minor planets, H10 for comets.
					 // This is the apparent visual magnitude when 1AU from sun and observer, with zero phase angle.
	Mat4d rotLocalToParent;          // GZ2015: was undocumented.
					 // Apparently this is the axis orientation with respect to the parent body. For planets, this is axis orientation w.r.t. VSOP87A/J2000 ecliptical system.
	float axisRotation;              // Rotation angle of the Planet on its axis, degrees.
					 // For Earth, this should be Greenwich Mean Sidereal Time GMST.
					 // For V0.15+, and for planets computed after the IAU2009 paper this is angle W (rotDeg),
					 // i.e. angle between ascending node of body equator w.r.t. ICRF equator and its prime meridian.
	StelTextureSP texMap;            // Planet map texture
	StelTextureSP normalMap;         // Planet normal map texture

	Ring* rings;                     // Planet rings
	double distance;                 // Temporary variable used to store the distance to a given point
					 // it is used for sorting while drawing
	float sphereScale;               // Artificial scaling for better viewing
	double lastJDE;                  // caches JDE of last positional computation
	// The callback for the calculation of the equatorial rect heliocentric position at time JDE.
	posFuncType coordFunc;
	void* userDataPtr;               // this is always set to an Orbit object usable for positional computations. For the major planets, it is NULL.

	OsculatingFunctType *const osculatingFunc;
	QSharedPointer<Planet> parent;           // Planet parent i.e. sun for earth
	QList<QSharedPointer<Planet> > satellites;      // satellites of the Planet
	LinearFader hintFader;
	LinearFader labelsFader;         // Store the current state of the label for this planet
	bool flagLabels;                 // Define whether labels should be displayed
	bool hidden;                     // useful for fake planets used as observation positions - not drawn or labeled
	bool atmosphere;                 // Does the planet have an atmosphere?
	bool halo;                       // Does the planet have a halo?	
	PlanetType pType;                // Type of body

	ApparentMagnitudeAlgorithm vMagAlgorithm;

	static Vec3f labelColor;
	static StelTextureSP hintCircleTex;	
	static QMap<PlanetType, QString> pTypeMap; // Maps fast type to english name.
	static QMap<ApparentMagnitudeAlgorithm, QString> vMagAlgorithmMap;

	static bool flagCustomGrsSettings;	// Is enabled usage of custom settings for calculation of position of Great Red Spot?
	static double customGrsJD;		// Initial JD for calculation of position of Great Red Spot
	static int customGrsLongitude;		// Longitude of Great Red Spot (System II, degrees)
	static double customGrsDrift;		// Annual drift of Great Red Spot position (degrees)
	
	// 0.16: Axes of planets and moons require terms depending on T=(jde-J2000)/36525, described in Explanatory Supplement 2013, Tables 10.1 and 10.10-14.
	// Others require frequent updates, depending on jde-J2000. (Moon etc.)
	// These should be updated as frequently as needed, optimally with the planet. Light time correction should be applied when needed.
	// beat place to call update is the SolarSystem::computePlanets()
	struct PlanetCorrections {
		double JDE_E; // keep record of when these values are valid: Earth
		double JDE_J; // keep record of when these values are valid: Jupiter
		double JDE_S; // keep record of when these values are valid: Saturn
		double JDE_U; // keep record of when these values are valid: Uranus
		double JDE_N; // keep record of when these values are valid: Neptune

		double E1; // Earth corrections. These are from WGCCRE2009.
		double E2;
		double E3;
		double E4;
		double E5;
		double E6;
		double E7;
		double E8;
		double E9;
		double E10;
		double E11;
		double E12;
		double E13;
		double Ja1; // Jupiter axis terms, Table 10.1
		double Ja2;
		double Ja3;
		double Ja4;
		double Ja5;
		double Na; // Neptune axix term
		double J1; // corrective terms for Jupiter' moons, Table 10.10
		double J2;
		double J3;
		double J4;
		double J5;
		double J6;
		double J7;
		double J8;
		double S1; // corrective terms for Saturn's moons, Table 10.12
		double S2;
		double S3;
		double S4;
		double S5;
		double S6;
		//double U1; // corrective terms for Uranus's moons, Table 10.14.
		//double U2;
		//double U3;
		//double U4;
		//double U5;
		//double U6;
		//double U7;
		//double U8;
		//double U9;
		//double U10;
		double U11;
		double U12;
		double U13;
		double U14;
		double U15;
		double U16;
		double N1; // corrective terms for Neptune's moons, Table 10.15 (N=Na!)
		double N2;
		double N3;
		double N4;
		double N5;
		double N6;
		double N7;
	};
	static PlanetCorrections planetCorrections;
	// Update the planet corrections. planet=3(Moon), 5(Jupiter), 6(Saturn), 7(Uranus), 8(Neptune).
	// The values are immediately converted to radians.
	static void updatePlanetCorrections(const double JDE, const int planet);

	// Shader-related variables
	struct PlanetShaderVars {
		int projectionMatrix;
		int texCoord;
		int unprojectedVertex;
		int vertex;
		int texture;
		int lightDirection;
		int eyeDirection;
		int diffuseLight;
		int ambientLight;
		int shadowCount;
		int shadowData;
		int sunInfo;
		int skyBrightness;
		
		void initLocations(QOpenGLShaderProgram*);
	};
	static PlanetShaderVars planetShaderVars;
	static QOpenGLShaderProgram* planetShaderProgram;

	// Shader-related variables
	struct RingPlanetShaderVars : public PlanetShaderVars {
		// Rings-specific variables
		int isRing;
		int ring;
		int outerRadius;
		int innerRadius;
		int ringS;
	};
	static RingPlanetShaderVars ringPlanetShaderVars;
	static QOpenGLShaderProgram* ringPlanetShaderProgram;
	
	struct MoonShaderVars : public PlanetShaderVars {
		// Moon-specific variables
		int earthShadow;
		int normalMap;
	};
	static MoonShaderVars moonShaderVars;
	static QOpenGLShaderProgram* moonShaderProgram;
	
	static void initShader();
	static void deinitShader();
};

#endif // _PLANET_HPP_

