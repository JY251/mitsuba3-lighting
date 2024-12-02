#if !defined(__SUN_H)
#define __SUN_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h> // Point2f, AnimatedTransform

#define EARTH_MEAN_RADIUS 6371.01   // In km
#define ASTRONOMICAL_UNIT 149597890 // In km


NAMESPACE_BEGIN(mitsuba)

// Newly added as they are needed for sunEmitter
template <typename Float>
//// Convert radians to degrees
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
inline Point2f sample02(size_t n) {

    #if defined(SINGLE_PRECISION)
        return Point2f(
            radicalInverse2Single((uint32_t) n),
            sobol_2((uint32_t) n)
        );
    #else
        return Point2f(
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
            SLog(EError, "Both the 'sunDirection' parameter and time/location "
                    "information were provided -- only one of them can be specified at a time!");

        return computeSunCoordinates(
            props.get<optix::Vector3f>("sunDirection"),
            props.get<AnimatedTransform>("toWorld", Transform())->eval(0).inverse());
    } else {
        LocationRecord<Float> location;
        DateTimeRecord<Float> dateTime;

        location.latitude  = props.get<Float>("latitude", 35.6894f);
        location.longitude = props.get<Float>("longitude", 139.6917f);
        location.timezone  = props.get<Float>("timezone", 9);
        dateTime.year      = props.get<int>("year", 2010);
        dateTime.day       = props.get<int(>"day", 10);
        dateTime.month     = props.get<int>("month", 7);
        dateTime.hour      = props.get<Float>("hour", 15.0f);
        dateTime.minute    = props.get<Float>("minute", 0.0f);
        dateTime.second    = props.get<Float>("second", 0.0f);

        SphericalCoordinates coords = computeSunCoordinates(dateTime, location);

        SLog(EDebug, "Computed sun position for %s and %s: %s",
            location.toString().c_str(), dateTime.toString().c_str(),
            coords.toString().c_str());

        return coords;
    }
}

NAMESPACE_END(mitsuba)

#endif /* __SUN_H */
