#include <mitsuba/core/fwd.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/sensor.h>

NAMESPACE_BEGIN(mitsuba)


// MI_VARIANT class: docs/src/developer_guide/variants_cpp.rst
MI_VARIANT class MyIrradianceMeter final: public Sensor<Float, Spectrum> {
	public:
		MI_IMPORT_BASE(Sensor, m_film, m_shape, sample_wavelengths, m_to_world)
		MI_IMPORT_TYPES(Shape)

		MyIrradianceMeter(const Properties &props) : Base(props) {
        // if (props.has_property("to_world"))
        //     Throw("Found a 'to_world' transformation -- this is not allowed. "
        //           "The irradiance meter inherits this transformation from its parent "
        //           "shape.");

        // if (m_film->rfilter()->radius() > .5f + math::RayEpsilon<Float>)
        //     Log(Warn, "This sensor should only be used with a reconstruction filter"
        //        "of radius 0.5 or lower (e.g. default 'box' filter)");
            // m_center = props.get<Vector3f>("center");
            // m_width = props.get<Float32>("width"); // width in x-direction
            // m_height = props.get<Float32>("height"); // height in y-direction
            // m_normal = props.get<Vector3f>("normal");
		}
    
		// sample_ray_differential
    // NOTE: not to be confused with std::pair<RayDifferential3f, Spectrum>
    std::pair<Ray3f, Spectrum>
    sample_ray(Float time, Float wavelength_sample,
                            const Point2f & sample2,
                            const Point2f & sample3,
                            Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Sample spatial component
        PositionSample3f ps = dr::zeros<PositionSample3f>(); // not dr::zeros is also OK (as currently we do not need differentials)

        // m_center = Vector3f(0, 0, 0); // As `m_center` is protected, we can't access it directly
        // ps.p = m_center + Vector3f(1 - 2*sample2.x(), 1 - 2*sample2.y(), 0);
        ps.p = Vector3f(1 - 2*sample2.x(), 1 - 2*sample2.y(), 0);
        ps.p = -1*ps.p; // As the sensor is facing downwards, we need to flip the direction of the sensor (in this way, we obtain the result with no right-left flip and no upside-down flip)
        ps.n = Vector3f(0, 0, 1);

        // value(): GPU type data; scalar(): CPU type data
        ps.p = m_to_world.value().transform_affine(ps.p);
        ps.n = m_to_world.value().transform_affine(ps.n);
        
        // 2. Sample directional component
        // Vector3f local = warp::square_to_cosine_hemisphere(sample3) * M_PI; # Multiply by M_PI here does not make sense
        Vector3f local = warp::square_to_cosine_hemisphere(sample3);

        // 3. Sample spectrum
        auto [wavelengths, wav_weight] =
            sample_wavelengths(dr::zeros<SurfaceInteraction3f>(),
                               wavelength_sample,
                               active);

        Vector3f d = Frame3f(ps.n).to_world(local);
        Point3f o = ps.p + d * math::RayEpsilon<Float>; // RayEpsilon is done to avoid zero-division error by adding a small value to the origin

        Ray3f ray;
        ray.time = time;
        ray.wavelengths = wavelengths;
        ray.o = o;
        ray.d = d;

        return { ray, wav_weight*M_PI };
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        return { m_shape->sample_direction(it, sample, active), dr::Pi<ScalarFloat> };
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        return m_shape->pdf_direction(it, ds, active);
    }

    Spectrum eval(const SurfaceInteraction3f &/*si*/, Mask /*active*/) const override {
        return dr::Pi<ScalarFloat> / m_shape->surface_area();
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "IrradianceMeter[" << std::endl << "  surface_area = ";

        if (m_shape)
            oss << m_shape->surface_area();
        else
            oss << " <no shape attached!>";
        oss << "," << std::endl;

        oss << "  film = " << indent(m_film) << "," << std::endl << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

protected:
    // Vector2f m_height;
    // Vector2f m_width;

    // In 3D space, a rectangle is determined if we know 3 points
    Vector3f m_center; // center of the sensor
    Float32 m_width; // width of the sensor
    Float32 m_height; // height of the sensor
    Vector3f m_normal; // normal of the sensor
};

MI_IMPLEMENT_CLASS_VARIANT(MyIrradianceMeter, Sensor)
MI_EXPORT_PLUGIN(MyIrradianceMeter, "MyIrradianceMeter");
NAMESPACE_END(mitsuba)
