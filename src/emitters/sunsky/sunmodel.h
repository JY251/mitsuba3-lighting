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



NAMESPACE_END(mitsuba)

#endif /* __SUN_H */
