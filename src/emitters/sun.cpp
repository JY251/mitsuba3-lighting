/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// #include <mitsuba/render/scene.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/transform.h> // ScalarTransform4f

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/qmc.h>
// #include "sunsky/sunmodel.h"

NAMESPACE_BEGIN(mitsuba)

/* Apparent radius of the sun as seen from the earth (in degrees).
   This is an approximation--the actual value is somewhere between
   0.526 and 0.545 depending on the time of year */
#define SUN_APP_RADIUS 0.5358

// Revise the following as we only focussing on RGB color space
// #if SPECTRUM_SAMPLES == 3
# define SUN_PIXELFORMAT Bitmap::PixelFormat::RGB // in include/mitsuba/core/bitmap.h, we have RGB instead of ERGB; We have no ESpectrum
// #else
// # define SUN_PIXELFORMAT Bitmap::ESpectrum
// #endif

/*!\plugin{sun}{Sun emitter}
 * \icon{emitter_sun}
 * \order{7}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of aerosol present in the atmosphere.
 *         Valid range: 2-10. \default{3, corresponding to a clear sky in a temperate climate}
 *     }
 *     \parameter{year, month, day}{\Integer}{Denote the date of the
 *      observation \default{2010, 07, 10}}
 *     \parameter{hour,minute,\showbreak second}{\Float}{Local time
 *       at the location of the observer in 24-hour format\default{15, 00, 00,
 *       i.e. 3PM}}
 *     \parameter{latitude, longitude, timezone}{\Float}{
 *       These three parameters specify the oberver's latitude and longitude
 *       in degrees, and the local timezone offset in hours, which are required
 *       to compute the sun's position. \default{35.6894, 139.6917, 9 --- Tokyo, Japan}
 *     }
 *     \parameter{sunDirection}{\Vector}{Allows to manually
 *       override the sun direction in world space. When this value
 *       is provided, parameters pertaining to the computation
 *       of the sun direction (\code{year, hour, latitude,} etc.
 *       are unnecessary. \default{none}
 *     }
 *     \parameter{resolution}{\Integer}{Specifies the horizontal resolution of the precomputed
 *         image that is used to represent the sun environment map \default{512, i.e. 512$\times$256}}
 *     \parameter{scale}{\Float}{
 *         This parameter can be used to scale the amount of illumination
 *         emitted by the sun emitter. \default{1}
 *     }
 *     \parameter{sunRadiusScale}{\Float}{
 *         Scale factor to adjust the radius of the sun, while preserving its power.
 *         Set to \code{0} to turn it into a directional light source.
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 * This plugin implements the physically-based sun model proposed by
 * Preetham et al. \cite{Preetham1999Practical}. Using the provided position
 * and time information (see \pluginref{sky} for details), it can determine the
 * position of the sun as seen from the position of the observer.
 * The radiance arriving at the earth surface is then found based on the spectral
 * emission profile of the sun and the extinction cross-section of the
 * atmosphere (which depends on the \code{turbidity} and the zenith angle of the sun).
 *
 * Like the \code{blackbody} emission profile (Page~\pageref{sec:blackbody}),
 * the sun model introduces physical units into the rendering process.
 * The radiance values computed by this plugin have units of power ($W$) per
 * unit area ($m^{-2}$) per steradian ($sr^{-1}$) per unit wavelength ($nm^{-1}$).
 * If these units are inconsistent with your scene description, you may use the
 * optional \texttt{scale} parameter to adjust them.
 *
 * This plugin supplies proper spectral power distributions when Mitsuba is
 * compiled in spectral rendering mode. Otherwise, they are simply projected onto
 * a linear RGB color space.
 *
 * \remarks{
 *   \item The sun is an intense light source that subtends a tiny solid angle.
 *   This can be a problem for certain rendering techniques (e.g. path
 *   tracing), which produce high variance output (i.e. noise in renderings)
 *   when the scene also contains specular or glossy or materials.
 * }
 */

#define EARTH_MEAN_RADIUS 6371.01   // In km
#define ASTRONOMICAL_UNIT 149597890 // In km


// Newly added as they are needed for sunEmitter
template <typename Float>
using Point2f = Point<Float, 2>; // I am not sure where it is defined, so I use this...
//// Convert radians to degrees
template <typename Float>
Float radToDeg(Float value) { return value * (static_cast<Float>(180.0f) / static_cast<Float>(M_PI)); }

template <typename Float>
/// Convert degrees to radians
Float degToRad(Float value) { return value * (static_cast<Float>(M_PI) / static_cast<Float>(180.0f)); }

template <typename Float>
struct SphericalCoordinates {
    Float elevation;
    Float azimuth;

    inline SphericalCoordinates() { }

    inline SphericalCoordinates(Float elevation, Float azimuth)
        : elevation(elevation), azimuth(azimuth) { }

		// As stream is not required in mitsuba3, the following constructor is not required
    // inline SphericalCoordinates(Stream *stream) {
    //     elevation = stream->readFloat();
    //     azimuth = stream->readFloat();
    // }

    // void serialize(Stream *stream) const {
    //     stream->writeFloat(elevation);
    //     stream->writeFloat(azimuth);
    // }

    // Not to add line such as "Float radToDeg(Float value);" before this line will cause an error.
    // Whether you add <Float> or not, it works. 
    // radToDeg<Float> does not work. But make this as just "radToDeg" works.
    std::string toString() const {
        std::ostringstream oss;
        oss << "SphericalCoordinates[elevation = " << radToDeg<Float>(elevation)
            << ", azimuth = " << radToDeg<Float>(azimuth) << "]";
        return oss.str();
    }
};

template <typename Float>
optix::Vector3f toSphere(const SphericalCoordinates<Float> coords) {
    Float sinTheta, cosTheta, sinPhi, cosPhi;

    dr::sincos(coords.elevation, &sinTheta, &cosTheta);
    dr::sincos(coords.azimuth, &sinPhi, &cosPhi);

    return optix::Vector3f(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta);
}

template <typename Float>
struct DateTimeRecord {
    int year;
    int month;
    int day;
    Float hour;
    Float minute;
    Float second;

    std::string toString() const {
        std::ostringstream oss;
        oss << "DateTimeRecord[year = " << year
            << ", month= " << month
            << ", day = " << day
            << ", hour = " << hour
            << ", minute = " << minute
            << ", second = " << second << "]";
        return oss.str();
    }
};

template <typename Float>
struct LocationRecord {
    Float longitude;
    Float latitude;
    Float timezone;

    std::string toString() const {
        std::ostringstream oss;
        oss << "LocationRecord[latitude = " << latitude
            << ", longitude = " << longitude
            << ", timezone = " << timezone << "]";
        return oss.str();
    }
};

template <typename Float>
SphericalCoordinates<Float> computeSunCoordinates(const DateTimeRecord<Float> &dateTime, const LocationRecord<Float> &location) {
    // Main variables
    double elapsedJulianDays, decHours;
    double eclipticLongitude, eclipticObliquity;
    double rightAscension, declination;
    double elevation, azimuth;

    // Auxiliary variables
    double dY;
    double dX;

    /* Calculate difference in days between the current Julian Day
       and JD 2451545.0, which is noon 1 January 2000 Universal Time */
    {
        // Calculate time of the day in UT decimal hours
        decHours = dateTime.hour - location.timezone +
            (dateTime.minute + dateTime.second / 60.0 ) / 60.0;

        // Calculate current Julian Day
        int liAux1 = (dateTime.month-14) / 12;
        int liAux2 = (1461*(dateTime.year + 4800 + liAux1)) / 4
            + (367 * (dateTime.month - 2 - 12 * liAux1)) / 12
            - (3 * ((dateTime.year + 4900 + liAux1) / 100)) / 4
            + dateTime.day - 32075;
        double dJulianDate = (double) liAux2 - 0.5 + decHours / 24.0;

        // Calculate difference between current Julian Day and JD 2451545.0
        elapsedJulianDays = dJulianDate - 2451545.0;
    }

    /* Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
       ecliptic in radians but without limiting the angle to be less than 2*Pi
       (i.e., the result may be greater than 2*Pi) */
    {
        double omega = 2.1429 - 0.0010394594 * elapsedJulianDays;
        double meanLongitude = 4.8950630 + 0.017202791698 * elapsedJulianDays; // Radians
        double anomaly = 6.2400600 + 0.0172019699 * elapsedJulianDays;

        eclipticLongitude = meanLongitude + 0.03341607 * std::sin(anomaly)
            + 0.00034894 * std::sin(2*anomaly) - 0.0001134
            - 0.0000203 * std::sin(omega);

        eclipticObliquity = 0.4090928 - 6.2140e-9 * elapsedJulianDays
            + 0.0000396 * std::cos(omega);
    }

    /* Calculate celestial coordinates ( right ascension and declination ) in radians
       but without limiting the angle to be less than 2*Pi (i.e., the result may be
       greater than 2*Pi) */
    {
        double sinEclipticLongitude = std::sin(eclipticLongitude);
        dY = std::cos(eclipticObliquity) * sinEclipticLongitude;
        dX = std::cos(eclipticLongitude);
        rightAscension = std::atan2(dY, dX);
        if (rightAscension < 0.0)
            rightAscension += 2*M_PI;
        declination = std::asin(std::sin(eclipticObliquity) * sinEclipticLongitude);
    }

    // Calculate local coordinates (azimuth and zenith angle) in degrees
    {
        double greenwichMeanSiderealTime = 6.6974243242
            + 0.0657098283 * elapsedJulianDays + decHours;

        double localMeanSiderealTime = degToRad((Float) ((greenwichMeanSiderealTime * 15
            + location.longitude)));

        double latitudeInRadians = degToRad(location.latitude);
        double cosLatitude = std::cos(latitudeInRadians);
        double sinLatitude = std::sin(latitudeInRadians);

        double hourAngle = localMeanSiderealTime - rightAscension;
        double cosHourAngle = std::cos(hourAngle);

        elevation = std::acos(cosLatitude * cosHourAngle
            * std::cos(declination) + std::sin(declination) * sinLatitude);

        dY = -std::sin(hourAngle);
        dX = std::tan(declination) * cosLatitude - sinLatitude * cosHourAngle;

        azimuth = std::atan2(dY, dX);
        if (azimuth < 0.0)
            azimuth += 2*M_PI;

        // Parallax Correction
        elevation += (EARTH_MEAN_RADIUS / ASTRONOMICAL_UNIT) * std::sin(elevation);
    }

    return SphericalCoordinates((Float) elevation, (Float) azimuth);
}

/// Van der Corput radical inverse in base 2 with single precision
inline float radicalInverse2Single(uint32_t n, uint32_t scramble = 0U) {
    /* Efficiently reverse the bits in 'n' using binary operations */
#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))) || defined(__clang__)
    n = __builtin_bswap32(n);
#else
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
#endif
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);

    // Account for the available precision and scramble
    n = (n >> (32 - 24)) ^ (scramble & ~-(1 << 24));

    return (float) n / (float) (1U << 24);
}

/// Van der Corput radical inverse in base 2 with double precision
inline double radicalInverse2Double(uint64_t n, uint64_t scramble = 0ULL) {
    /* Efficiently reverse the bits in 'n' using binary operations */
#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))) || defined(__clang__)
    n = __builtin_bswap64(n);
#else
    n = (n << 32) | (n >> 32);
    n = ((n & 0x0000ffff0000ffffULL) << 16) | ((n & 0xffff0000ffff0000ULL) >> 16);
    n = ((n & 0x00ff00ff00ff00ffULL) << 8)  | ((n & 0xff00ff00ff00ff00ULL) >> 8);
#endif
    n = ((n & 0x0f0f0f0f0f0f0f0fULL) << 4)  | ((n & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
    n = ((n & 0x3333333333333333ULL) << 2)  | ((n & 0xccccccccccccccccULL) >> 2);
    n = ((n & 0x5555555555555555ULL) << 1)  | ((n & 0xaaaaaaaaaaaaaaaaULL) >> 1);

    // Account for the available precision and scramble
    n = (n >> (64 - 53)) ^ (scramble & ~-(1LL << 53));

    return (double) n / (double) (1ULL << 53);
}

// sample02 and sample02Double are included in qmc.h in mitsuba1, but not in mitsuba3. Therefore, I include them here.
/// Generate an element from a (0, 2) sequence (without scrambling)

template <typename Float>
inline Point2f<Float> sample02(size_t n) {

    #if defined(SINGLE_PRECISION)
        return Point2f<Float>(
            radicalInverse2Single((uint32_t) n),
            sobol_2((uint32_t) n)
        );
    #else
        return Point2f<Float>(
            radicalInverse2Double((uint64_t) n),
            sobol_2((uint64_t) n)
        );
    #endif
}


template <typename Float>
SphericalCoordinates<Float> computeSunCoordinates(const Properties &props) {
    /* configure position of sun */
    if (props.has_property("sunDirection")) {
        if (props.has_property("latitude") || props.has_property("longitude")
            || props.has_property("timezone") || props.has_property("day")
            || props.has_property("time"))
            Throw("Both the 'sunDirection' parameter and time/location "
                    "information were provided -- only one of them can be specified at a time!");

        return computeSunCoordinates(
            props.get<optix::Vector3f>("sunDirection"),
            (ScalarTransform4f) props.get<ScalarTransform4f>("to_world", ScalarTransform4f()))
    } else {
        LocationRecord<Float> location;
        DateTimeRecord<Float> dateTime;

        location.latitude  = props.get<Float>("latitude", 35.6894f);
        location.longitude = props.get<Float>("longitude", 139.6917f);
        location.timezone  = props.get<Float>("timezone", 9);
        dateTime.year      = props.get<int>("year", 2010);
        dateTime.day       = props.get<int>("day", 10);
        dateTime.month     = props.get<int>("month", 7);
        dateTime.hour      = props.get<Float>("hour", 15.0f);
        dateTime.minute    = props.get<Float>("minute", 0.0f);
        dateTime.second    = props.get<Float>("second", 0.0f);

        SphericalCoordinates coords = computeSunCoordinates(dateTime, location);

        // SLog(EDebug, "Computed sun position for %s and %s: %s",
        //     location.toString().c_str(), dateTime.toString().c_str(),
        //     coords.toString().c_str());

        return coords;
    }
}

/* The following is from the implementation of "A Practical Analytic Model for
   Daylight" by A.J. Preetham, Peter Shirley, and Brian Smits */

/* All data lifted from MI. Units are either [] or cm^-1. refer when in doubt MI */

// k_o Spectrum table from pg 127, MI.
template <typename Float> 
Float k_oWavelengths[64] = {
    300, 305, 310, 315, 320, 325, 330, 335, 340, 345,
    350, 355, 445, 450, 455, 460, 465, 470, 475, 480,
    485, 490, 495, 500, 505, 510, 515, 520, 525, 530,
    535, 540, 545, 550, 555, 560, 565, 570, 575, 580,
    585, 590, 595, 600, 605, 610, 620, 630, 640, 650,
    660, 670, 680, 690, 700, 710, 720, 730, 740, 750,
    760, 770, 780, 790
};

template <typename Float>
Float k_oAmplitudes[65] = {
    10.0, 4.8, 2.7, 1.35, .8, .380, .160, .075, .04, .019, .007,
    .0, .003, .003, .004, .006, .008, .009, .012, .014, .017,
    .021, .025, .03, .035, .04, .045, .048, .057, .063, .07,
    .075, .08, .085, .095, .103, .110, .12, .122, .12, .118,
    .115, .12, .125, .130, .12, .105, .09, .079, .067, .057,
    .048, .036, .028, .023, .018, .014, .011, .010, .009,
    .007, .004, .0, .0
};

// k_g Spectrum table from pg 130, MI.
template <typename Float>
Float k_gWavelengths[4] = {
    759, 760, 770, 771
};

template <typename Float>
Float k_gAmplitudes[4] = {
    0, 3.0, 0.210, 0
};

// k_wa Spectrum table from pg 130, MI.
template <typename Float>
Float k_waWavelengths[13] = {
    689, 690, 700, 710, 720,
    730, 740, 750, 760, 770,
    780, 790, 800
};

template <typename Float>
Float k_waAmplitudes[13] = {
    0, 0.160e-1, 0.240e-1, 0.125e-1,
    0.100e+1, 0.870, 0.610e-1, 0.100e-2,
    0.100e-4, 0.100e-4, 0.600e-3,
    0.175e-1, 0.360e-1
};

/* Wavelengths corresponding to the table below */
template <typename Float>
Float solWavelengths[38] = {
    380, 390, 400, 410, 420, 430, 440, 450,
    460, 470, 480, 490, 500, 510, 520, 530,
    540, 550, 560, 570, 580, 590, 600, 610,
    620, 630, 640, 650, 660, 670, 680, 690,
    700, 710, 720, 730, 740, 750
};

/* Solar amplitude in watts / (m^2 * nm * sr) */
template <typename Float>
Float solAmplitudes[38] = {
    16559.0, 16233.7, 21127.5, 25888.2, 25829.1,
    24232.3, 26760.5, 29658.3, 30545.4, 30057.5,
    30663.7, 28830.4, 28712.1, 27825.0, 27100.6,
    27233.6, 26361.3, 25503.8, 25060.2, 25311.6,
    25355.9, 25134.2, 24631.5, 24173.2, 23685.3,
    23212.1, 22827.7, 22339.8, 21970.2, 21526.7,
    21097.9, 20728.3, 20240.4, 19870.8, 19427.2,
    19072.4, 18628.9, 18259.2
};

template <typename Float, typename Spectrum>
Spectrum computeSunRadiance(Float theta, Float turbidity) {
    Spectrum k_oCurve(k_oWavelengths<Float>, k_oAmplitudes<Float>, 64);
    Spectrum k_gCurve(k_gWavelengths<Float>, k_gAmplitudes<Float>, 4);
    Spectrum k_waCurve(k_waWavelengths<Float>, k_waAmplitudes<Float>, 13);
    Spectrum solCurve(solWavelengths<Float>, solAmplitudes<Float>, 38);
    Float data[91], wavelengths[91];  // (800 - 350) / 5  + 1

    Float beta = 0.04608365822050f * turbidity - 0.04586025928522f;

    // Relative Optical Mass
    Float m = 1.0f / (std::cos(theta) + 0.15f *
        std::pow(93.885f - theta/M_PI*180.0f, (Float) -1.253f));

    Float lambda;
    int i = 0;
    for (i = 0, lambda = 350; i < 91; i++, lambda += 5) {
        // Rayleigh Scattering
        // Results agree with the graph (pg 115, MI) */
        Float tauR = dr::exp(-m * 0.008735f * std::pow(lambda/1000.0f, (Float) -4.08));

        // Aerosol (water + dust) attenuation
        // beta - amount of aerosols present
        // alpha - ratio of small to large particle sizes. (0:4,usually 1.3)
        // Results agree with the graph (pg 121, MI)
        const Float alpha = 1.3f;
        Float tauA = dr::exp(-m * beta * std::pow(lambda/1000.0f, -alpha));  // lambda should be in um

        // Attenuation due to ozone absorption
        // lOzone - amount of ozone in cm(NTP)
        // Results agree with the graph (pg 128, MI)
        const Float lOzone = .35f;
        Float tauO = dr::exp(-m * k_oCurve.eval(lambda) * lOzone);

        // Attenuation due to mixed gases absorption
        // Results agree with the graph (pg 131, MI)
        Float tauG = dr::exp(-1.41f * k_gCurve.eval(lambda) * m / std::pow(1 + 118.93f
            * k_gCurve.eval(lambda) * m, (Float) 0.45f));

        // Attenuation due to water vapor absorbtion
        // w - precipitable water vapor in centimeters (standard = 2)
        // Results agree with the graph (pg 132, MI)
        const Float w = 2.0;
        Float tauWA = dr::exp(-0.2385f * k_waCurve.eval(lambda) * w * m /
                std::pow(1 + 20.07f * k_waCurve.eval(lambda) * w * m, (Float) 0.45f));

        data[i] = solCurve.eval(lambda) * tauR * tauA * tauO * tauG * tauWA;
        wavelengths[i] = lambda;
    }

    Spectrum interpolated(wavelengths, data, 91);
    Spectrum discretized;
    discretized.fromContinuousSpectrum(interpolated);
    discretized.clampNegative();

    return discretized;
}


template <typename Float, typename Spectrum>
class SunEmitter final : public Emitter<Float, Spectrum> {
public:
    // Without this, error on line 116: 'Base' was not declared in this scope
    // I have recognized that "m_to_world" is used in mitsuba3 instead of "m_worldtransform" in mitsuba
    // I have also recognized that "sample_wavelengths" is used in mitsuba3 instead of "m_samplingWeight" in mitsuba
    MI_IMPORT_BASE(Emitter, m_to_world, sample_wavelengths) 
    MI_IMPORT_TYPES(Scene, Texture)

    SunEmitter(const Properties &props)
            : Base(props) {
        m_scale = props.get<Float>("scale", 1.0f);
        m_resolution = props.get<int>("resolution", 512);
        m_sun = computeSunCoordinates<Float>(props);
        m_sunRadiusScale = props.get<Float>("sunRadiusScale", 1.0f);
        m_turbidity = props.get<Float>("turbidity", 3.0f);
        m_stretch = props.get<Float>("stretch", 1.0f);
    }

    // SunEmitter(Stream *stream, InstanceManager *manager)
    //         : Emitter(stream, manager) {
    //     m_scale = stream->readFloat();
    //     m_sunRadiusScale = stream->readFloat();
    //     m_turbidity = stream->readFloat();
    //     m_resolution = stream->readInt();
    //     m_sun = SphericalCoordinates(stream);
    //     configure();
    // }

    // `serialize` is not required in mitsuba3
    // Stream and InstanceManager (defined in `include/mitsuba/core/serialization.h` in mitsuba1) is also related to serialization, so they are not required in mitsuba3
    // void serialize(Stream *stream, InstanceManager *manager) const {
    //     Emitter::serialize(stream, manager);
    //     stream->writeFloat(m_scale);
    //     stream->writeFloat(m_sunRadiusScale);
    //     stream->writeFloat(m_turbidity);
    //     stream->writeInt(m_resolution);
    //     m_sun.serialize(stream);
    // }

    void configure() {
        SphericalCoordinates sun(m_sun);
        sun.elevation *= m_stretch;
        m_sunDir = toSphere(sun);

        /* Solid angle covered by the sun */
        m_theta = degToRad(SUN_APP_RADIUS * 0.5f);
        m_solidAngle = 2 * M_PI * (1 - std::cos(m_theta));
        m_radiance = computeSunRadiance(m_sun.elevation, m_turbidity) * m_scale;
    }

    bool isCompound() const {
        return true;
    }

    ref<Emitter<Float, Spectrum>> getElement(size_t i) {
        if (i != 0)
            return NULL;

        if (m_sunRadiusScale == 0) {
            Properties props("directional");
            const Transform3f &trafo = m_to_world.value(); // NOTE: we have transform3f and transform4f, so I am not sure which one to use
            props.set_array3f("direction", -trafo(m_sunDir)); // setVector is not available in mitsuba3; As trafo seems to be 3D, set_array3f is used
            props.set_float("samplingWeight", sample_wavelengths);

            props.set_array3f("irradiance", m_radiance * m_solidAngle); // as we only use RGB as the spectrum, setSpectrum is defined as set_array3f in mitsuba3

            // Emitter<Float, Spectrum> *emitter = static_cast<Emitter *>(
            //     PluginManager::instance()->create_object<Emitter>(props));

            ref<Emitter<Float, Spectrum>> emitter = PluginManager::instance()->create_object<Emitter<Float, Spectrum>>(props); // ref to line 116 of ior.h

            emitter->configure();
            return emitter;
        }

        /* Rasterizing the sphere to an environment map and checking the
           individual pixels for coverage (which is what Mitsuba 0.3.0 did)
           was slow and not very effective; for instance the power varied
           dramatically with resolution changes. Since the sphere generally
           just covers a few pixels, the code below rasterizes it much more
           efficiently by generating a few thousand QMC samples.

           Step 1: compute a *very* rough estimate of how many
           pixel in the output environment map will be covered
           by the sun */
        size_t pixelCount = m_resolution*m_resolution/2;
        Float cosTheta = std::cos(m_theta * m_sunRadiusScale);

        /* Ratio of the sphere that is covered by the sun */
        Float coveredPortion = 0.5f * (1 - cosTheta);

        /* Approx. number of samples that need to be generated,
           be very conservative */
        size_t nSamples = (size_t) std::max((Float) 100,
            (pixelCount * coveredPortion * 1000));

        ref<Bitmap> bitmap = new Bitmap(SUN_PIXELFORMAT, m_component_format,
            Vector2i(m_resolution, m_resolution/2));
        bitmap->clear();
        Frame3f frame(m_sunDir);

        Point2f factor(bitmap->width() / (2*M_PI), // getWidth/Height is mitsuba1 is replaced by width/height in mitsuba3 (in include/mitsuba/core/bitmap.h)
            bitmap->height() / M_PI);

        Spectrum *target = (Spectrum *) bitmap->data(); // data() is used in mitsuba3 instead of getFloatData() in mitsuba1
        Spectrum value =
            m_radiance * (2 * M_PI * (1-std::cos(m_theta))) *
            static_cast<Float>(bitmap->width() * bitmap->height())
            / (2 * M_PI * M_PI * nSamples);

        for (size_t i=0; i<nSamples; ++i) {
            Vector dir = frame.toWorld(
                warp::square_to_uniform_cone(cosTheta, sample02(i)));

            Float sinTheta = drjit::safe_sqrt(1-dir.y*dir.y);
            SphericalCoordinates sphCoords = fromSphere(dir);

            Point2i pos(
                std::min(std::max(0, (int) (sphCoords.azimuth * factor.x)), bitmap->width()-1),
                std::min(std::max(0, (int) (sphCoords.elevation * factor.y)), bitmap->height()-1));

            target[pos.x + pos.y * bitmap->width()] += value / std::max((Float) 1e-3f, sinTheta);
        }

        /* Instantiate a nested envmap plugin */
        Properties props("envmap");
        Properties::Data bitmapData;
        bitmapData.ptr = (uint8_t *) bitmap.get();
        bitmapData.size = sizeof(Bitmap);
        props.setData("bitmap", bitmapData);
        props.setAnimatedTransform("toWorld", m_to_world);
        props.set_float("samplingWeight", sample_wavelengths);
        ref<Emitter<Float, Spectrum>> emitter = PluginManager::instance()->create_object<Emitter<Float, Spectrum>>(props);
        emitter->configure();
        return emitter;
    }

    // The following functions are not necesasry so they are commented out
    // AABB getAABB() const {
    //     NotImplementedError("getAABB");
    // }

    // std::string toString() const {
    //     std::ostringstream oss;
    //     oss << "SunEmitter[" << std::endl
    //         << "  sunDir = " << m_sunDir.toString() << "," << std::endl
    //         << "  sunRadiusScale = " << m_sunRadiusScale << "," << std::endl
    //         << "  turbidity = " << m_turbidity << "," << std::endl
    //         << "  scale = " << m_scale << std::endl
    //         << "]";
    //     return oss.str();
    // }

    MI_DECLARE_CLASS()
protected:
    /// Environment map resolution
    int m_resolution;
    /// Constant scale factor applied to the model
    Float m_scale;
    /// Scale factor that can be applied to the sun radius
    Float m_sunRadiusScale;
    /// Angle cutoff for the sun disk (w/o scaling)
    Float m_theta;
    /// Solid angle covered by the sun (w/o scaling)
    Float m_solidAngle;
    /// Position of the sun in spherical coordinates
    SphericalCoordinates<Float> m_sun;
    /// Direction of the sun (untransformed)
    optix::Vector3f m_sunDir;
    /// Turbidity of the atmosphere
    Float m_turbidity;
    /// Radiance arriving from the sun disk
    Spectrum m_radiance;
    /// Stretch factor to extend to the bottom hemisphere
    Float m_stretch;
};

MI_IMPLEMENT_CLASS_VARIANT(SunEmitter, Emitter)
MI_EXPORT_PLUGIN(SunEmitter, "Sun emitter");
NAMESPACE_END(mitsuba)

