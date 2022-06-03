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

protected:

	// energies
	double m_se;

	// bulk viscosity damping parameters
	double m_bv_c1;
	double m_bv_c2;

};
