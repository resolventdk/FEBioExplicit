#pragma once
#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Class for storing material data relevant for explicit dynamics 
class FEExplicitMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEExplicitMaterialPoint(FEMaterialPoint* pt);

	//! copy
	FEMaterialPoint* Copy();

	//! initialization
	void Init(bool bflag);

	//! Serialization
	void Serialize(DumpStream& dmp);

public:

	double m_h;  // characteristic lenght of element

};

//-----------------------------------------------------------------------------
class FEExplicitMaterial
{
public:
	virtual FEElasticMaterial* GetElasticMaterial() = 0;
};

//-----------------------------------------------------------------------------
//! This material adds a bulk viscosity pressure to stress tensor 
class FEExplicitElastic : public FEElasticMaterial, public FEExplicitMaterial
{
public:
	//! constructor
	FEExplicitElastic(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;

	//! return the elastic material
	FEElasticMaterial* GetElasticMaterial() override { return m_mat; }

public:
	//! Cauchy stress 
	mat3ds Stress(FEMaterialPoint& mp) override;

	tens4ds Tangent(FEMaterialPoint& mp) override { return m_mat->Tangent(mp); }

	double WaveSpeed(FEMaterialPoint& mp);

protected:

private:
	FEElasticMaterial*	m_mat;	// elastic base material
	FEParamDouble       m_bv_c1;  // bulk viscosity damping linear parameter
	FEParamDouble       m_bv_c2;  // bulk viscosisty damping quadratic parameter

	DECLARE_FECORE_CLASS();
};
