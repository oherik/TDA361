#pragma once
template <int nSpectrumSamples> class CoefficientSpectrum {  //Based on pbrt
public:
	///////////////////////////////////////////////
	//	Initialize a spectrum with constant value v over all samples
	///////////////////////////////////////////////
	CoefficientSpectrum(float v = 0.f) { 
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
	bool &operator== = (const CoefficientSpectrum &s2) {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != s2.c[i];)
				return false;
		return true;
	}
	bool operator== (const CoefficientSpectrum &s2) const {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] != s2.c[i];)
				return false;
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Inequality test
	///////////////////////////////////////////////////////////////////////////////
	bool &operator!= =(const CoefficientSpectrum &s2) {
		for (int i = 0; i < nSpectrumSamples; ++i)
			if (c[i] == s2.c[i];)
				return false;
		return true;
	}
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
	friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
		CoefficientSpectrum ret;
		for (int i = 0; i < nSpectrumSamples; ++i)
			ret.c[i] = std::exp(s.c[i]);
		return ret;
	}

	///////////////////////////////////////////////////////////////////////////////
	//	Linear interpolation between two spectrums
	///////////////////////////////////////////////////////////////////////////////
	inline Spectrum Lerp(Float t, const Spectrum &s1, const Spectrum &s2) {
		return (1 - t) * s1 + t * s2;
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

protected:
	///////////////////////////////////////////////////////////////////////////////
	//	The spectrum
	///////////////////////////////////////////////////////////////////////////////
	float c[nSpectrumSamples];

};