#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include <functional>
#include <algorithm>    // std::sort

#ifdef PBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // PBRT_FLOAT_AS_DOUBLE


class SampledSpectrum;
template <int nSamples> class CoefficientSpectrum;
typedef RGBSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;

static  Float Lerp(Float t, Float v1, Float v2) {
	return (1 - t) * v1 + t * v2;
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
static const Float CIE_Y_integral = 106.856895;


Float AverageSpectrumSamples(const Float *lambda, const Float *vals,
	int n, Float lambdaStart, Float lambdaEnd) {
	if (lambdaEnd <= lambda[0]) return vals[0];
	if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
	if (n == 1) return vals[0];
	Float sum = 0;
	if (lambdaStart < lambda[0])
		sum += vals[0] * (lambda[0] - lambdaStart);
	if (lambdaEnd > lambda[n - 1])
		sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);
	int i = 0;
	while (lambdaStart > lambda[i + 1]) ++i;
	auto interp = [lambda, vals](Float w, int i) {
		return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]),
			vals[i], vals[i + 1]);
	};
	for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
		Float segLambdaStart = std::max(lambdaStart, lambda[i]);
		Float segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
		sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
			(segLambdaEnd - segLambdaStart);
	}
	return sum / (lambdaEnd - lambdaStart);
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
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] += s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] += s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Subtranct one spectral dsitribution from another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator-=(const CoefficientSpectrum &s2) {
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] -= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] -= s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Unary negation of a distribution
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum operator-() const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = -ret.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Multiply one spectral dsitribution with another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator*=(const CoefficientSpectrum &s2) {
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] *= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator*(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] *= s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Divide one spectral distribution with another
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum &operator/=(const CoefficientSpectrum &s2) {
		for (int i = 0; i < nSpectrumSamples; ++i)
			c[i] /= s2.c[i];
		return *this;
	}
	CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
		CoefficientSpectrum ret = *this;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] /= s2.c[i];
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Equality test
	///////////////////////////////////////////////////////////////////////////////
	bool operator== (const CoefficientSpectrum &s2) const {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != s2.c[i])
				return false;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Inequality test
	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const CoefficientSpectrum &s2) const {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] == s2.c[i];)
				return false;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	See if it's completely black
	///////////////////////////////////////////////////////////////////////////////
	bool IsBlack() const {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != 0.) return false;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Square root of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::sqrt(s.c[i]);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Power of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s, float exp) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::pow(s.c[i], exp);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	The base-e exponential of a spectrum
	///////////////////////////////////////////////////////////////////////////////
	friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::exp(s.c[i]);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Clamping a spectrum to defined range
	///////////////////////////////////////////////////////////////////////////////
	CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = ::Clamp(c[i], low, high);
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
	//	Check if the samples are sorted by wavelength
	///////////////////////////////////////////////////////////////////////////////
	static bool SpectrumSamplesSorted(const Float *lambda, const Float *v, int n) {
		std::vector<Float> slambda(&lambda[0], &lambda[n]);
		for (int i = 1; i < nSpectralSamples; ++i) {
			if (slambda[i] < slambda[i - 1]) return false;
		}
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Sort samples by wavelength
	///////////////////////////////////////////////////////////////////////////////
	static int SortSpectrumSamples(const Float *lambda, const Float *v, int n) {
		std::vector<Float> slambda(&lambda[0], &lambda[n]);
		std::vector<Float> sv(&v[0], &v[n]);

		sortVecPair(slambda, sv, less_than());
	}


	struct less_than {
		inline bool operator() (const std::pair<Float, Float>& a, const std::pair<Float, Float>& b)
		{
			return (a.first < b.first);
		}
	};

	template <typename T, typename R, typename Compare>
	static int sortVecPair(std::vector<T>& vecA, std::vector<R>& vecB, Compare cmp) {

		std::vector<pair<T, R>> vecC;
		vecC.reserve(vecA.size());
		for (int i = 0; i<vecA.size(); i++)
		{
			vecC.push_back(std::make_pair(vecA[i], vecB[i]));
		}

		std::sort(vecC.begin(), vecC.end(), cmp);

		vecA.clear();
		vecB.clear();
		vecA.reserve(vecC.size());
		vecB.reserve(vecC.size());
		for (int i = 0; i<vecC.size(); i++)
		{
			vecA.push_back(vecC[i].first);
			vecB.push_back(vecC[i].second);
		}
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
		Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
			Float(CIE_Y_integral * nSpectralSamples);
		xyz[0] *= scale;
		xyz[1] *= scale;
		xyz[2] *= scale;
	}

	static void Init() {
		for (int i = 0; i < nSpectralSamples; ++i) {
			Float wl0 = Lerp(Float(i) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
			X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0, wl1);
			Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0, wl1);
			Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0, wl1);
		}
		//Compute XYZ matching functions for SampledSpectrum 324
		//	Compute RGB to spectrum functions for SampledSpectrum
	}

private:

	static SampledSpectrum X, Y, Z;
};




///////////////////////////////////////////////////////////////////////////////
//	Linear interpolation between two spectrums
///////////////////////////////////////////////////////////////////////////////
inline Spectrum Lerp(Float t, const Spectrum &s1, const Spectrum &s2) {
	return s1 * (1 - t) + s2 * t;
}
