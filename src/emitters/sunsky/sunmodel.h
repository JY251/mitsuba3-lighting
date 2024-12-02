#if !defined(__SUN_H)
#define __SUN_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/util.h> // For radToDeg

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


NAMESPACE_END(mitsuba)

#endif /* __SUN_H */
