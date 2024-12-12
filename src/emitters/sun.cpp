#include <mitsuba/mitsuba.h> // std::string
#include <math.h> // M_PI
#include <mitsuba/render/fwd.h> // MI_EXPORT, MI_IMPORT
// #include <mitsuba/core/platform.h> // MI_EXPORT, MI_IMPORT
// #include <string.h> // This lines does not appear in mitsuba3 except for ext files 

NAMESPACE_BEGIN(mitsuba)

/* Apparent radius of the sun as seen from the earth (in degrees).
   This is an approximation--the actual value is somewhere between
   0.526 and 0.545 depending on the time of year */
#define SUN_APP_RADIUS 0.5358

#if SPECTRUM_SAMPLES == 3
# define SUN_PIXELFORMAT Bitmap::ERGB
#else
# define SUN_PIXELFORMAT Bitmap::ESpectrum
#endif

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

/////////////////////////////////////// Header ///////////////////////////////////////
// MI_IMPORT_BASE()
// MI_IMPORT_TYPES(Frame3f, Spectrum)

template <typename Float>
struct SphericalCoordinates {
    Float elevation;
    Float azimuth;

    inline SphericalCoordinates() { }

		// Q. Not sure about the following constructor
    inline SphericalCoordinates(Float elevation, Float azimuth)
        : elevation(elevation), azimuth(azimuth) { }

		// As stream is not used in mitsuba3, this function is not required
    // inline SphericalCoordinates(Stream *stream) {
    //     elevation = stream->readFloat();
    //     azimuth = stream->readFloat();
    // }

    // void serialize(Stream *stream) const {
    //     stream->writeFloat(elevation);
    //     stream->writeFloat(azimuth);
    // }

    // Q. This does not work well and this is not essential, so I commented it out
    // std::string toString() const {
    //     std::ostringstream oss;
    //     oss << "SphericalCoordinates[elevation = " << radToDeg(elevation)
    //         << ", azimuth = " << radToDeg(azimuth) << "]";
    //     return oss.str();
    // }
};

/**
 * \brief Abstract continous spectral power distribution data type,
 * which supports evaluation at arbitrary wavelengths.
 *
 * Here, the term 'continous' doesn't necessarily mean that the
 * underlying spectrum is continous, but rather emphasizes the fact
 * that it is a function over the reals (as opposed to the discrete
 * spectrum, which only stores samples for a discrete set of wavelengths).
 *
//  * \ingroup libpython
//  * \ingroup libcore
 */

// I did not add MI_EXPORT/MI_IMPORT as this is used only inside the file
// This is a base class; the functions are "virtual" (they are "empty") and will be overridden in the derived class
template <typename Float>
class ContinuousSpectrum {
public:
    /**
     * Evaluate the value of the spectral power distribution
     * at the given wavelength.
     *
     * \param lambda  A wavelength in nanometers
     */
    virtual Float eval(Float lambda) const = 0;

    /**
     * \brief Integrate the spectral power distribution
     * over a given interval and return the average value
     *
     * Unless overridden in a subclass, the integration is done
     * using adaptive Gauss-Lobatto quadrature.
     *
     * \param lambdaMin
     *     The lower interval bound in nanometers
     *
     * \param lambdaMax
     *     The upper interval bound in nanometers
     *
     * \remark If \c lambdaMin >= \c lambdaMax, the
     *     implementation will return zero.
     */
    virtual Float average(Float lambdaMin, Float lambdaMax) const;

    /// \brief Return a string representation
    virtual std::string toString() const = 0;

    /// Virtual destructor
    virtual ~ContinuousSpectrum() { }
};

/**
 * \brief This spectral power distribution is defined as the
 * product of two other continuous spectra.
 */
template <typename Float>
class ProductSpectrum : public ContinuousSpectrum<Float> {
public:
    /** \brief Return the value of the spectral power distribution
     * at the given wavelength.
     */
    ProductSpectrum(const ContinuousSpectrum<Float> &s1,
        const ContinuousSpectrum<Float> &s2) : m_spec1(s1),
        m_spec2(s2) { }

    /** \brief Return the value of the spectral power distribution
     * at the given wavelength.
     */
    virtual Float eval(Float lambda) const;

    /// Virtual destructor
    virtual ~ProductSpectrum() { }

    /// Return a string representation
    std::string toString() const;
private:
    const ContinuousSpectrum<Float> &m_spec1;
    const ContinuousSpectrum<Float> &m_spec2;
};

/**
 * \brief Linearly interpolated spectral power distribution
 *
 * This class implements a linearly interpolated spectral
 * power distribution that is defined over a discrete set of
 * measurements at different wavelengths. Outside of the
 * specified range, the spectrum is assumed to be zero. Hence,
 * at least two entries are required to produce a nonzero
 * spectrum.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
template <typename Float>
class InterpolatedSpectrum : public ContinuousSpectrum<Float> {
public:
    /**
     * \brief Create a new interpolated spectrum with space
     * for the specified number of samples
     */
    InterpolatedSpectrum(size_t size = 0);

    /**
     * \brief Create a interpolated spectrum instance from
     * a float array
     */
    InterpolatedSpectrum(const Float *wavelengths,
        const Float *values, size_t nEntries);

    /**
     * \brief Read an interpolated spectrum from a simple
     * ASCII format.
     *
     * Each line of the file should contain an entry of the form
     * \verbatim
     * <wavelength in nm> <value>
     * \endverbatim
     * Comments preceded by '#' are also valid.
     */
    InterpolatedSpectrum(const fs::path &path);

    /**
     * \brief Append an entry to the spectral power distribution.
     *
     * Entries must be added in order of increasing wavelength
     */
    void append(Float lambda, Float value);

    /**
     * \brief This function adds a zero entry before and after
     * the stored wavelength range.
     *
     * This is useful when handling datasets that don't fall
     * off to zero at the ends. The spacing of the added entries
     * is determined by computing the average spacing of the
     * existing samples.
     */
    void zeroExtend();

    /// Clear all stored entries
    void clear();

    /**
     * \brief Return the value of the spectral power distribution
     * at the given wavelength.
     */
    Float eval(Float lambda) const;

    /**
     * \brief Integrate the spectral power distribution
     * over a given interval and return the average value
     *
     * This method overrides the implementation in
     * \ref ContinousSpectrum, since the integral can be
     * analytically computed for linearly interpolated spectra.
     *
     * \param lambdaMin
     *     The lower interval bound in nanometers
     *
     * \param lambdaMax
     *     The upper interval bound in nanometers
     *
     * \remark If \c lambdaMin >= \c lambdaMax, the
     *     implementation will return zero.
     */
    Float average(Float lambdaMin, Float lambdaMax) const;

    /// \brief Return a string representation
    std::string toString() const;

		// Q. I don't really understand the meaning of the following function
    /// Virtual destructor
    virtual ~InterpolatedSpectrum() { }
protected:
    std::vector<Float> m_wavelengths, m_values;
};

/////////////////////////////////////// Value ///////////////////////////////////////
static const int CIE_samples = 471;
template <typename Float>
const Float CIE_wavelengths[CIE_samples] = {
    360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373,
    374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
    388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401,
    402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415,
    416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429,
    430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443,
    444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457,
    458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471,
    472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485,
    486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499,
    500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513,
    514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527,
    528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541,
    542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555,
    556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569,
    570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583,
    584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597,
    598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611,
    612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625,
    626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639,
    640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653,
    654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667,
    668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681,
    682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695,
    696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709,
    710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723,
    724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737,
    738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751,
    752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765,
    766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779,
    780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793,
    794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807,
    808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821,
    822, 823, 824, 825, 826, 827, 828, 829, 830 };
template <typename Float>
const Float CIE_X_entries[CIE_samples];
template <typename Float>
const Float CIE_Y_entries[CIE_samples];
template <typename Float>
const Float CIE_Z_entries[CIE_samples];
/// @{ \name Interpolated versions of the XYZ color matching functions
template <typename Float>
static InterpolatedSpectrum CIE_X_interp(CIE_wavelengths<Float>, CIE_X_entries<Float>, CIE_samples);
template <typename Float>
static InterpolatedSpectrum CIE_Y_interp(CIE_wavelengths<Float>, CIE_Y_entries<Float>, CIE_samples);
template <typename Float>
static InterpolatedSpectrum CIE_Z_interp(CIE_wavelengths<Float>, CIE_Z_entries<Float>, CIE_samples);
/// @}


/////////////////////////////////////// Supplymentary Function ///////////////////////////////////////
template <typename Spectrum>
constexpr auto extract_color_dimension() {
    return std::tuple_size<typename std::decay<Spectrum>::type>::value;
}

template <typename T>
struct ExtractFirstTemplateArg;

template <template <typename, size_t> class Color, typename T, size_t N>
struct ExtractFirstTemplateArg<Color<T, N>> {
    using type = T; // First template argument
};

template <typename Float, typename Spectrum>
void fromContinuousSpectrum(Spectrum &spectrum, const ContinuousSpectrum<Float> &smooth) {
#if SPECTRUM_SAMPLES == 3
    /* Convolve with the XYZ matching functions and convert to RGB */
    Float start = CIE_wavelengths[0], end = CIE_wavelengths[CIE_samples-1];

    Float X = ProductSpectrum(smooth, CIE_X_interp).average(start, end);
    Float Y = ProductSpectrum(smooth, CIE_Y_interp).average(start, end);
    Float Z = ProductSpectrum(smooth, CIE_Z_interp).average(start, end);
    Float normalization = 1.0f / CIE_Y_interp.average(start, end);

    X *= normalization; Y *= normalization; Z *= normalization;

    fromXYZ(X, Y, Z);
#else
    /* Spectral rendering mode -- average over each bin */

    constexpr auto SPECTRUM_SAMPLES = extract_color_dimension<Spectrum>();
    using SpectrumScalarType = typename ExtractFirstTemplateArg<Spectrum>::type;

    SpectrumScalarType s[SPECTRUM_SAMPLES]; // line 603 of spectrum.h (We assume Spectrum = Color<Float, *>. If it is not Float, we need to extract the first argument of Spectrum = Color(*, *) and use it as type )
    // smooth.average is stored in s; s stores "smooth" i.e., ContinuousSpectrum; This is something like "Color<Float, 3>"
    // Type is something like "Float"

    for (int i=0; i<SPECTRUM_SAMPLES; i++)
        // s[i] is used to store the spectral values corresponding to different wavelength intervals for spectral rendering
        // m_wavelength belong to "ContinuousSpectrum" (protected)
        s[i] = smooth.average(spectrum.m_wavelengths[i], spectrum.m_wavelengths[i+1]);
#endif
}


template <typename Float>
Float radToDeg(Float value) { return value * (180.0f / M_PI); }

template <typename Float>
Float degToRad(Float value) { return value * (M_PI / 180.0f); }
		



/////////////////////////////////////// Main ///////////////////////////////////////
template <typename Float, typename Spectrum>
class SunEmitter final : public Emitter<Float, Spectrum> {
public:
		// MI_IMPORT_BASE();
		// MI_IMPORT_TYPES();
        MI_IMPORT_CORE_TYPES() // Frame3f


		SunEmitter(const Properties &props): Base(props) {
				m_scale = props.get<Float>("scale", 1.0f);
				m_resolution = props.get<int>("resolution", 512);
				m_sun = computeSunCoordinates<Float>(props);
				m_sunRadiusScale = props.get<Float>("sunRadiusScale", 1.0f);
				m_turbidity = props.get<Float>("turbidity", 3.0f);
				m_stretch = props.get<Float>("stretch", 1.0f);

				m_sun = SphericalCoordinates<Float>(props);
				
				// The followings are done in "configure" function in mituba1

				SphericalCoordinates<Float> sun(m_sun);
				sun.elevation *= m_stretch;
				m_sunDir = toSphere<Float>(sun);

        /* Solid angle covered by the sun */
        m_theta = degToRad<Float>(SUN_APP_RADIUS * 0.5f);
        m_solidAngle = 2 * M_PI * (1 - dr::cos(m_theta));
        m_radiance = computeSunRadiance(m_sun.elevation, m_turbidity) * m_scale;
		}


        Frame3f toSphere(const SphericalCoordinates<Float> coords) {
                // converts spherical coordinates to cartesian coordinates
            Float sinTheta, cosTheta, sinPhi, cosPhi;

            dr::sincos(coords.elevation, &sinTheta, &cosTheta);
            dr::sincos(coords.azimuth, &sinPhi, &cosPhi);

            return Frame3f(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta);
        }

        Spectrum computeSunRadiance(Float theta, Float turbidity) {
            InterpolatedSpectrum k_oCurve(k_oWavelengths, k_oAmplitudes, 64);
            InterpolatedSpectrum k_gCurve(k_gWavelengths, k_gAmplitudes, 4);
            InterpolatedSpectrum k_waCurve(k_waWavelengths, k_waAmplitudes, 13);
            InterpolatedSpectrum solCurve(solWavelengths, solAmplitudes, 38);
                Float data[91], wavelengths[91];  // (800 - 350) / 5  + 1
            Float beta = 0.04608365822050f * turbidity - 0.04586025928522f;

            // Relative Optical Mass
            Float m = 1.0f / (std::cos(theta) + 0.15f *
                dr::pow(93.885f - theta/M_PI*180.0f, (Float) -1.253f));

            Float lambda;
            int i = 0;
            for (i = 0, lambda = 350; i < 91; i++, lambda += 5) {
                // Rayleigh Scattering
                // Results agree with the graph (pg 115, MI) */
                Float tauR = dr::exp(-m * 0.008735f * dr::pow(lambda/1000.0f, (Float) -4.08));

                // Aerosol (water + dust) attenuation
                // beta - amount of aerosols present
                // alpha - ratio of small to large particle sizes. (0:4,usually 1.3)
                // Results agree with the graph (pg 121, MI)
                const Float alpha = 1.3f;
                Float tauA = dr::exp(-m * beta * dr::pow(lambda/1000.0f, -alpha));  // lambda should be in um

                // Attenuation due to ozone absorption
                // lOzone - amount of ozone in cm(NTP)
                // Results agree with the graph (pg 128, MI)
                const Float lOzone = .35f;
                Float tauO = dr::exp(-m * k_oCurve.eval(lambda) * lOzone);

                // Attenuation due to mixed gases absorption
                // Results agree with the graph (pg 131, MI)
                Float tauG = dr::exp(-1.41f * k_gCurve.eval(lambda) * m / dr::pow(1 + 118.93f
                    * k_gCurve.eval(lambda) * m, (Float) 0.45f));

                // Attenuation due to water vapor absorbtion
                // w - precipitable water vapor in centimeters (standard = 2)
                // Results agree with the graph (pg 132, MI)
                const Float w = 2.0;
                Float tauWA = dr::exp(-0.2385f * k_waCurve.eval(lambda) * w * m /
                        dr::pow(1 + 20.07f * k_waCurve.eval(lambda) * w * m, (Float) 0.45f));

                data[i] = solCurve.eval(lambda) * tauR * tauA * tauO * tauG * tauWA;
                wavelengths[i] = lambda;
            }

            InterpolatedSpectrum interpolated(wavelengths, data, 91);
            Spectrum discretized;
            discretized.fromContinuousSpectrum(interpolated);
            discretized.clampNegative();

            return discretized;
        }

};
