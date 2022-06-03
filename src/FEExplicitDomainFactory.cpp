#include "FEExplicitDomainFactory.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FEElasticShellDomain.h"
#include "FEBioMech/FEElasticTrussDomain.h"
#include "FEBioMech/FERigidSolidDomain.h"
#include "FEBioMech/FERigidShellDomain.h"
#include "FEBioMech/FERemodelingElasticDomain.h"
#include "FEBioMech/FEUDGHexDomain.h"
#include "FEBioMech/FEUT4Domain.h"
#include "FEBioMech/FE3FieldElasticSolidDomain.h"
#include "FEBioMech/FEDiscreteElasticDomain.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FEBioMech/FEDiscreteElementMaterial.h"
#include "FEBioMech/FESRIElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FEExplicitDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	const char* sztype = 0;
	FE_Element_Class eclass = spec.eclass;
	FE_Element_Shape eshape = spec.eshape;
	FE_Element_Type etype = spec.etype;
	if (dynamic_cast<FERigidMaterial*>(pmat))
	{
		// rigid elements
		if      (eclass == FE_ELEM_SOLID) sztype = "rigid-solid";
		else if (eclass == FE_ELEM_SHELL) 
		{
			if (spec.m_shell_formulation == OLD_SHELL) sztype = "rigid-shell-old";
			else sztype = "rigid-shell";
		}
		else return 0;
	}
	else if (dynamic_cast<FEElasticMaterial*>(pmat))
	{
		// flexible elements
		sztype = "explicit-solid";
	} else {  // all other configurations are currently not supported!
		return 0;
	}

	if (sztype)
	{
		FEDomain* pd = fecore_new<FEDomain>(sztype, pfem);
		if (pd) pd->SetMaterial(pmat);
		return pd;
	}
	else return 0;
}
