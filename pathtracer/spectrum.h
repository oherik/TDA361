#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include <functional>
#include <float.h>
#include <cmath>
#include <algorithm>    // std::sort
#include <glm/glm.hpp> //vec3

#ifdef PBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // PBRT_FLOAT_AS_DOUBLE

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif

#define EPSILON 0.0001f

class SampledSpectrum;
class RGBSpectrum;
template <int nSamples> class CoefficientSpectrum;
//typedef RGBSpectrum Spectrum;
typedef SampledSpectrum Spectrum;
enum class SpectrumType { Reflectance, Illuminant };

extern bool SpectrumSamplesSorted(const float *lambda, const float *vals, int n);
extern void SortSpectrumSamples(float *lambda, float *vals, int n);
extern float InterpolateSpectrumSamples(const float *lambda, const float *vals, int n, float l);
extern float AverageSpectrumSamples(const float *lambda, const float *vals, int n, float lambdaStart, float lambdaEnd);

// General functions
static  Float Lerp(Float t, Float v1, Float v2) {
	return (1 - t) * v1 + t * v2;
}

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60; //Should be enough for the visible spectrum mvh Erik

//XYZ stuff
static const int nCIESamples = 471;
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
extern const Float CIE_lambda[nCIESamples];
static const Float CIE_Y_integral = 106.856895; // integral over Y(lambda) * d(lambda)

//Smith's method for converting RGB to SPD
static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];

// D65 spectral power, sunlight-ish
extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];


///////////////////////////////////////////////////////////////////////////////
//	Conversion stuff. Uses wavelengths defined by some kind of HD-TV
///////////////////////////////////////////////////////////////////////////////
inline void XYZToRGB(const Float xyz[3], Float rgb[3]) {
	rgb[0] = 3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
	rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
	rgb[2] = 0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}

///////////////////////////////////////////////////////////////////////////////
//	The inverse of the above
///////////////////////////////////////////////////////////////////////////////
inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
	xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
	xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
	xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}


template <int nSpectrumSamples> class CoefficientSpectrum {  //Based on pbrt
public:
	///////////////////////////////////////////////
	//	Initialize a spectrum with constant value v over all samples
	///////////////////////////////////////////////
	CoefficientSpectrum(Float v = 0.f) {
		for (int i = 0; i < nSpectrumSamples; ++i) {
			c[i] = v;
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Add two spectral distributions
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
		////assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] += s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		//assert(!ret.HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] += s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Subtranct one spectral dsitribution from another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator-=(const CoefficientSpectrum &s2) {
		//assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] -= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		assert(!ret.HasNaNs());
		assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] -= s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Unary negation of a distribution
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum operator-() const {
		CoefficientSpectrum ret = *this;
		//assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = -ret.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Multiply one spectral dsitribution with another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator*=(const CoefficientSpectrum &s2) {
		//assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] *= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator*(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		//assert(!ret.HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] *= s2.c[i];
		return ret;
	}


	///////////////////////////////////////////////////////////////////////////////
	//	Multiply one spectral dsitribution with a float
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator*=(Float f) {
		//assert(!HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] *= f;
		return *this;
	}
	CoefficientSpectrum operator*(Float f) const {
		CoefficientSpectrum ret = *this;
		//assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] *= f;
		return ret;
	}
	friend inline CoefficientSpectrum operator*(Float f, const CoefficientSpectrum &s1) {
		return s1 * f;
	}


	///////////////////////////////////////////////////////////////////////////////
	//	Divide one spectral distribution with another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator/=(const CoefficientSpectrum &s2) {
		//assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] /= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		//assert(!ret.HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] /= s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Equality test
	///////////////////////////////////////////////////////////////////////////////
	bool operator== (const CoefficientSpectrum &s2) const {
		//assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != s2.c[i])
				return false;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Inequality test
	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const CoefficientSpectrum &s2) const {
		//assert(!HasNaNs());
		//assert(!s2.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i){
			if (c[i] == s2.c[i]){
				return false;
            }
        }
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	See if it's completely black
	///////////////////////////////////////////////////////////////////////////////
	bool IsBlack() const {
		//assert(!HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i) {
			if (c[i] > EPSILON) return false;
		}
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Get a mean diff between two spectrums
	///////////////////////////////////////////////////////////////////////////////
	Float diff(const CoefficientSpectrum &s) const {
		Float sum = 0;
		for (int i = 0; i < nSpectralSamples; ++i) {
			sum += c[i]-s.c[i];
		}
		return sum/ nSpectralSamples;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Square root of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
		CoefficientSpectrum ret;
		//////assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::sqrt(s.c[i]);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Power of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s, float exp) {
		CoefficientSpectrum ret;
		//assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::pow(s.c[i], exp);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	The base-e exponential of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
		CoefficientSpectrum ret;
		//assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::exp(s.c[i]);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Clamping a spectrum to defined range
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum Clamp(Float low = 0, Float high = INFINITY) const {
		CoefficientSpectrum ret;
		//assert(!ret.HasNaNs());
		for (int i = 0; i < nSpectrumSamples; ++i) {
			float before = c[i];
			ret.c[i] = ::Clamp(c[i], low, high);
		}

		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Check if a spectrum is erronous, eg it has a zero-division in it somewhere
	///////////////////////////////////////////////////////////////////////////////
	bool HasNaNs() const {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (std::isnan(c[i])) return true;
		return false;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Access sample value
	///////////////////////////////////////////////////////////////////////////////
	Float &operator[](int i) {
		return c[i];
	}


	///////////////////////////////////////////////////////////////////////////////
	//	Constant sample size
	///////////////////////////////////////////////////////////////////////////////
	static const int nSamples = nSpectrumSamples;


	//TODO make more general

protected:
	///////////////////////////////////////////////////////////////////////////////
	//	The spectrum
	///////////////////////////////////////////////////////////////////////////////
	Float c[nSpectrumSamples];


};



///////////////////////////////////////////////////////////////////////////////
//	An SPD with uniformly spaced samples in the visible spectrum
///////////////////////////////////////////////////////////////////////////////


class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
public:
	SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) { }
	SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v) : CoefficientSpectrum<nSpectralSamples>(v) { }
	///////////////////////////////////////////////////////////////////////////////
	//	Takes arrays of SPD sample values v at given wavelengths lambda and uses them to define a piecewise linear function to represent the SPD
	///////////////////////////////////////////////////////////////////////////////
	static SampledSpectrum FromSampled(const Float *lambda, const Float *v, int n) {
		// Sort samples if unordered
		if (!SpectrumSamplesSorted(lambda, v, n)) {
			std::vector<Float> slambda(&lambda[0], &lambda[n]);
			std::vector<Float> sv(&v[0], &v[n]);
			SortSpectrumSamples(&slambda[0], &sv[0], n);
			return FromSampled(&slambda[0], &sv[0], n);
		}

		SampledSpectrum r;
		for (int i = 0; i < nSpectralSamples; ++i) {
			Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
				sampledLambdaStart, sampledLambdaEnd);
			r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
		}
		return r;
	}


	

	///////////////////////////////////////////////////////////////////////////////
	//	Compute using Reimann sum
	///////////////////////////////////////////////////////////////////////////////
	void ToXYZ(Float xyz[3]) const {
		xyz[0] = xyz[1] = xyz[2] = 0.f;
		for (int i = 0; i < nSpectralSamples; ++i) {
			xyz[0] += X.c[i] * c[i];
			xyz[1] += Y.c[i] * c[i];
			xyz[2] += Z.c[i] * c[i];
		}
		Float scale = Float(sampledLambdaEnd - sampledLambdaStart) / Float(CIE_Y_integral * nSpectralSamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	y is sort of the luminance, good to have
	///////////////////////////////////////////////////////////////////////////////
	Float y() const {
		Float yy = 0.f;
		for (int i = 0; i < nSpectralSamples; ++i) {
			yy += Y.c[i] * c[i];
		}	
		return yy * Float(sampledLambdaEnd - sampledLambdaStart) / Float(CIE_Y_integral * nSpectralSamples);
	}
	
	///////////////////////////////////////////////////////////////////////////////
	//	Convert the samples to rgb values, and stores them in a provided Float[3]
	///////////////////////////////////////////////////////////////////////////////
	void ToRGB(Float rgb[3]) const {
		//assert(!HasNaNs());
		Float xyz[3];
		ToXYZ(xyz);
		XYZToRGB(xyz, rgb);
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Convert the samples to a vec3 containing RGB values
	///////////////////////////////////////////////////////////////////////////////
	glm::vec3 ToRGB() {
		Float rgb[3];
		ToRGB(rgb);
		return glm::vec3(rgb[0], rgb[1], rgb[2]);
	}
	
	RGBSpectrum ToRGBSpectrum() const;
	static SampledSpectrum FromRGB(const float rgb[3], SpectrumType type);
	static SampledSpectrum FromRGB(glm::vec3, SpectrumType type);

	static SampledSpectrum FromXYZ(const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
		//assert(!HasNaNs());
		Float rgb[3];
		XYZToRGB(xyz, rgb);
		return FromRGB(rgb, type);
	}

	SampledSpectrum(const RGBSpectrum &r, SpectrumType type = SpectrumType::Reflectance);

	static void Init() {
		for (int i = 0; i < nSpectralSamples; ++i) {
			Float wl0 = Lerp(Float(i) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0, wl1);
			Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0, wl1);
			Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0, wl1);
		}
		//	Compute RGB to spectrum functions for SampledSpectrum, code taken from PBRT
		for (int i = 0; i < nSpectralSamples; ++i) {
			Float wl0 = Lerp(Float(i) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			rgbRefl2SpectWhite.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectCyan.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,	nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectYellow.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectGreen.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen, nRGB2SpectSamples, wl0, wl1);
			rgbRefl2SpectBlue.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue, nRGB2SpectSamples, wl0, wl1);
			
			rgbIllum2SpectWhite.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectCyan.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectMagenta.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectYellow.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectRed.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectGreen.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen, nRGB2SpectSamples, wl0, wl1);
			rgbIllum2SpectBlue.c[i] = AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue, nRGB2SpectSamples, wl0, wl1);
		}
	}
	

private:
	static SampledSpectrum X, Y, Z;
	static float yint;
	static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
	static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
	static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
	static SampledSpectrum rgbRefl2SpectBlue;
	static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
	static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
	static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
	static SampledSpectrum rgbIllum2SpectBlue;
};

///////////////////////////////////////////////////////////////////////////////
//	RGB spectrum
///////////////////////////////////////////////////////////////////////////////
class RGBSpectrum : public CoefficientSpectrum<3> {
	public:
		RGBSpectrum(Float v = 0.f) : CoefficientSpectrum<3>(v) { }
		RGBSpectrum(const CoefficientSpectrum<3> &v) : CoefficientSpectrum<3>(v) { }
	
		void ToRGB(Float *rgb) const {
			rgb[0] = c[0];
			rgb[1] = c[1];
			rgb[2] = c[2];
		}
		const RGBSpectrum &ToRGBSpectrum() const {
			return *this;
		}

		void ToXYZ(Float xyz[3]) const {
			RGBToXYZ(c, xyz);
		}

		static RGBSpectrum FromXYZ(const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
			Float rgb[3];
			RGBSpectrum rgbSpec;
			XYZToRGB(xyz, rgbSpec.c);
			return rgbSpec;
		}

		Float y() const {
			Float xyz[3];
			RGBToXYZ(c, xyz);
			return xyz[1];
		}

		static RGBSpectrum FromRGB(const float rgb[3], SpectrumType type = SpectrumType::Reflectance) {
			RGBSpectrum s;
			s.c[0] = rgb[0];
			s.c[1] = rgb[1];
			s.c[2] = rgb[2];
			return s;
		}

		static RGBSpectrum FromSampled(const Float *lambda, const Float *v,	int n) {
			if (!SpectrumSamplesSorted(lambda, v, n)) {
				std::vector<Float> slambda(&lambda[0], &lambda[n]);
				std::vector<Float> sv(&v[0], &v[n]);
				SortSpectrumSamples(&slambda[0], &sv[0], n);
				return FromSampled(&slambda[0], &sv[0], n);
			}
			Float xyz[3] = { 0, 0, 0 };
			for (int i = 0; i < nCIESamples; ++i) {
				Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
				xyz[0] += val * CIE_X[i];
				xyz[1] += val * CIE_Y[i];
				xyz[2] += val * CIE_Z[i];
			}
			Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /	Float(CIE_Y_integral * nCIESamples);
			xyz[0] *= scale;
			xyz[1] *= scale;
			xyz[2] *= scale;
			return FromXYZ(xyz);
		}

		
};



///////////////////////////////////////////////////////////////////////////////
//	Linear interpolation between two spectrums
///////////////////////////////////////////////////////////////////////////////
inline Spectrum Lerp(Float t, const Spectrum &s1, const Spectrum &s2) {
	return s1 * (1-t) + s2 * t;
}

///////////////////////////////////////////////////////////////////////////////
//	Some finding stuff
///////////////////////////////////////////////////////////////////////////////
template <typename Predicate> int FindInterval(int size,
	const Predicate &pred) {
	int first = 0, len = size;
	while (len > 0) {
		int half = len >> 1, middle = first + half;
		if (pred(middle)) {
			first = middle + 1;
			len -= half + 1;
		}
		else
			len = half;
	}
	return Clamp(first - 1, 0, size - 2);
}

