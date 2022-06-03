#include "FEExplicitSolidDomain.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEIsotropicElastic.h"
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioMech/FEBioMech.h"

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEExplicitSolidDomain::FEExplicitSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) //, FESolidDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dofSU(pfem), m_dofV(pfem), m_dofSV(pfem), m_dofSA(pfem), m_dof(pfem)
{

	// bulk viscosity damping coefficients
	m_bv_c1 = 0.06;
	m_bv_c2 = 1.44;

	// zero energy sums
	m_se = 0.0;  

}

void FEExplicitSolidDomain::Update(const FETimeInfo& tp){

	// zero energy sums before summation
	m_se = 0.0;  

	// update ancestor (will update element stresses)
	FEElasticSolidDomain::Update(tp);

}


//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FEExplicitSolidDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double dt = tp.timeIncrement;
    
	FEIsotropicElastic* pmi = dynamic_cast<FEIsotropicElastic*>(m_pMat);
	if (!pmi){  // the actual elastic material is isotropic elastic
		feLogError("Only isotropic elastic material supported!");
	}

	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d r[NELN], v[NELN], a[NELN];
	GetCurrentNodalCoordinates(el, r);

	// update dynamic quantities
	if (m_update_dynamic)
	{
		for (int j = 0; j<neln; ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			v[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
			a[j] = node.m_at;
		}
	}

	// get characteristic element lenght in reference frame
	double h0 = pow(Volume(el), 0.333);

	// get gauss weights
	double* gw = el.GaussWeights();

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// material point coordinates
		pt.m_rt = el.Evaluate(r, n);

		// get the deformation gradient and determinant at intermediate time
        double Jt;
        mat3d Ft, Fp;
        Jt = defgrad(el, Ft, n);
        defgradp(el, Fp, n);

		pt.m_F = Ft;
		pt.m_J = Jt;


        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (dt > 1e-9) ? (Ft - Fp)*Fi / dt : mat3d(0.0); // avoid divison by zero at beginning of analysis
		if (m_update_dynamic)
		{
			pt.m_v = el.Evaluate(v, n);
			pt.m_a = el.Evaluate(a, n);
		}

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, tp);

		// get jacobian determinant in reference frame
		double J0 = detJ0(el, n);								
		
		// sum strain energy
		m_se += pmi->StrainEnergyDensity(mp)*J0*gw[n];		

		// get material parameters
		double mE = pmi->m_E(mp);
		double mv = pmi->m_v(mp);
		double dens = pmi->Density(mp);

		// lame parameters
		double lam = (mv*mE/((1+mv)*(1-2*mv)));
		double mu  = (0.5*mE/(1+mv));						
		
		// wavespeed
		double waveSpd = sqrt((lam+2.0*mu)/dens);
        
		// calculate the stress at this material point
//		pt.m_s = m_pMat->Stress(mp);
		pt.m_s = m_pMat->SolidStress(mp);
        
		// add bulk viscosity pressure
		double volumetricStrainRate = pt.RateOfDeformation().tr(); // tr(dEdt), E = 0.5*(L+L^T), L = duidxj
		double bv_pres = m_bv_c1*dens*h0*waveSpd*volumetricStrainRate;  // linear
		if (volumetricStrainRate < 0.0) {
			bv_pres -= m_bv_c2*dens*h0*h0*volumetricStrainRate*volumetricStrainRate;  // quadratic
		}
		pt.m_s += mat3dd(bv_pres);

    }
}