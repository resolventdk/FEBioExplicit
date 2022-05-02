#include "FEExplicitElastic.h"
#include "FEBioMech/FEIsotropicElastic.h"

#include "FECore/log.h"

//! constructor
FEExplicitMaterialPoint::FEExplicitMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) 
{ 

}

//! copy
FEMaterialPoint* FEExplicitMaterialPoint::Copy()
{
	return new FEExplicitMaterialPoint((m_pNext?m_pNext->Copy():0));
}

//! initialization
void FEExplicitMaterialPoint::Init(bool bflag)
{
	m_h = 0.0; // characteristic element lenght
}

//! Serialization
void FEExplicitMaterialPoint::Serialize(DumpStream& dmp)
{
	FEMaterialPoint::Serialize(dmp);
	dmp;// & m_bv_c1 & m_bv_c2;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEExplicitElastic, FEElasticMaterial)
	ADD_PROPERTY(m_mat, "elastic");
	ADD_PARAMETER(m_bv_c1, "bulk_viscosity_c1");
	ADD_PARAMETER(m_bv_c2, "bulk_viscosity_c2");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEExplicitElastic::FEExplicitElastic(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_mat = nullptr;
	m_bv_c1 = 0.06; // default value of linear coefficient
	m_bv_c2 = 1.44; // default value of quadratic coefficient
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEExplicitElastic::CreateMaterialPointData()
{ 
		FEMaterialPoint* pm = m_mat->CreateMaterialPointData();
		return new FEExplicitMaterialPoint(pm);
}

//-----------------------------------------------------------------------------
// Find maximum wavespeed
double FEExplicitElastic::WaveSpeed(FEMaterialPoint& mp)
{

	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());

	double dens = Density(mp);

	// evaluate the base material stress derivative wrt strain
	mat3ds E = (ep.RightCauchyGreen() - mat3dd(1.0))/2.0;
	tens4dmm dsde = m_mat->MaterialTangent(ep, E);  // what Strain to give as input
	double dsde_xx = dsde(0,0,0,0); // x-dir
	double dsde_yy = dsde(1,1,1,1); // y-dir
	double dsde_zz = dsde(2,2,2,2); // z-dir
	
	// get maximum in main directions
	// TODO use largest eigenvalue
	if ((dsde_xx >= dsde_yy) && (dsde_xx >= dsde_zz)) {
		return sqrt(dsde_xx / dens);
	} else if ( dsde_yy >= dsde_zz) {
		return sqrt(dsde_yy / dens);
	} else {
		return sqrt(dsde_zz / dens);
	}

}

//-----------------------------------------------------------------------------
// Take elastic stress and add bulk viscosity
// https://abaqus-docs.mit.edu/2017/English/SIMACAEGSARefMap/simagsa-c-ovwbulkvisc.htm
mat3ds FEExplicitElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEExplicitMaterialPoint& exp = *(mp.ExtractData<FEExplicitMaterialPoint>());

	// add bulk viscosity pressure
	double dens = Density(mp);
	double volumetricStrainRate = ep.RateOfDeformation().tr(); // tr(dEdt), E = 0.5*(L+L^T), L = duidxj
	double h = exp.m_h;
	double waveSpd = exp.m_c;
	double bv_pres = m_bv_c1(mp)*dens*h*waveSpd*volumetricStrainRate;  // linear
	if (volumetricStrainRate < 0.0) {
		bv_pres -= m_bv_c2(mp)*dens*h*h*volumetricStrainRate*volumetricStrainRate;  // quadratic
	}

	// return the base material stress + bulk viscosity pressure
	return m_mat->Stress(mp) + mat3dd(bv_pres);
}