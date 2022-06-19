#pragma once
#include "FEBioMech/FEElasticSolidDomain.h"


class FEBIOMECH_API FEExplicitSolidDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FEExplicitSolidDomain(FEModel* pfem);
	
	double StrainEnergy( ){ return m_se; }

public:  // overrides from ElasticSolidDomain

	// update the element stress
	virtual void UpdateElementStress(int iel, const FETimeInfo& tp);

public:  // overrides from MeshPartition

    virtual void Update(const FETimeInfo& tp);

private:

	// energies
	double m_se;

	// bulk viscosity damping parameters
	double m_dummy;  // to avoid memory corruption TODO figure out what happens
	double m_bv_c1;  // somehow this get value 0 when #pragma omp parallel for shared(NE, berr) in FEElasticSolidDomain::Update wtf
	double m_bv_c2;

	DECLARE_FECORE_CLASS();

};
