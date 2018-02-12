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
	//	The spectrum
	///////////////////////////////////////////////////////////////////////////////
	Float c[nSpectrumSamples];


//protected:
//	CoefficientSpectrum Protected Data 316
};