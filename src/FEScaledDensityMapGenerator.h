#pragma once
#include "FECore/FEDataGenerator.h"
#include "FECore/FEModelParam.h"

class FEScaledDensityMapGenerator : public FEDataGenerator
{
public:
	FEScaledDensityMapGenerator(FEModel* fem);
	~FEScaledDensityMapGenerator();

	bool Init() override;

	// generate the data array for the given element set
	bool Generate(FEDomainMap& data) override;

private:
	FEParamDouble	m_E;	    //!< Young's modulus
	FEParamDouble	m_v;	    //!< Poisson's ratio
	FEParamDouble	m_density0;	//!< Unscaled density
    double		    m_dt;	    //!< Target time step

	DECLARE_FECORE_CLASS();
};
