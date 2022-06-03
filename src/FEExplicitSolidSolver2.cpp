#include "FEExplicitSolidSolver2.h"
#include "FEExplicitSolidDomain.h"
#include "FEExplicitUtils.h"

#include "FEBioMech/FEBioMech.h"
#include "FEBioMech/FEBodyForce.h"
#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEElasticShellDomain.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FEIsotropicElastic.h"
#include "FEBioMech/FEMechModel.h"
#include "FEBioMech/FEResidualVector.h"
#include "FEBioMech/FEResidualVector.h"
#include "FEBioMech/FERigidBody.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/RigidBC.h"

#include "FECore/FEAnalysis.h"
#include "FECore/FELinearConstraintManager.h"
#include "FECore/FEMesh.h"
#include "FECore/FEModel.h"
#include "FECore/FENodalLoad.h"
#include "FECore/FESurfaceElementShape.h"
#include "FECore/FESurfaceLoad.h"
#include "FECore/log.h"

#include <math.h>       /* pow */

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEExplicitSolidSolver2, FESolver)
	ADD_PARAMETER(m_mass_lumping, "mass_lumping");
	ADD_PARAMETER(m_dyn_damping, "dyn_damping");
	ADD_PARAMETER(m_print_stride, "print_stride");
	ADD_PARAMETER(m_adjust_step_size_safety, FE_RANGE_RIGHT_OPEN(0.0, 1.0), "adjust_step_size_safety");
	ADD_PARAMETER(m_adjust_step_size_stride, "adjust_step_size_stride");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
//! ensure stable time step
bool time_step_limiter_cb(FEModel* pfem, unsigned int nwen, void* pd)
{
	// TODO: fix dt if no solids
	// https://abaqus-docs.mit.edu/2017/English/SIMACAEANLRefMap/simaanl-c-expdynamic.htm#simaanl-c-expdynamic-stability

	FEExplicitSolidSolver2* psolver = (FEExplicitSolidSolver2*) pd;

	// get current analysis step
	FEAnalysis* pstep = pfem->GetCurrentStep();

	// get number of ntimesteps completed
	int ntimesteps = pfem->GetCurrentStep()->m_ntimesteps;

	// get time step size safety and latest value of dtcrit
	double safety = psolver->m_adjust_step_size_safety;
	double dtcrit = psolver->m_dtcrit;

	// avoid printing to screen all the time, here beacuse UPDATE_TIME_CB is called before status print
	if (ntimesteps % psolver->m_print_stride == 0){
		pfem->UnBlockLog();
	} else {
		pfem->BlockLog();  // TODO: still print errors
	}

	// check if we want to re-estimate dtcrit
	if (ntimesteps % psolver->m_adjust_step_size_stride != 0){
		// enforce time step size unless next point is a must point and dt < dtcrit*safety
		if ((pstep->m_dt < safety*dtcrit) && (pstep->m_timeController) && (pstep->m_timeController->m_nmust >= 0)) {
			return true;
		} else {	
			pstep->m_dt = safety*dtcrit;	
			feLogEx(pfem, "Setting time step size: %7.2e\n", safety*dtcrit);
		}
		return true;  // return without re-estimating dtcrit
	}

	// estimate global dtcrit
	FEMesh& mesh = pfem->GetMesh();
	double globalMinDT = 1e99;
	int globalCritElId = 0;

	// print header
	feLogEx(pfem, "\n");
	feLogEx(pfem, "Stable time step size estimation\n");
	feLogEx(pfem, "=================================================\n");
	feLogEx(pfem, "Domain              Wavespeed Stepsize CritElemId\n");

	// loop all domains
	bool any_evaluated_domains = false;
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{		

		double domMinWaveSpd = 1e99;
		double domMinDT = 1e99;	
		int domCritElId = 0;

		// check whether it is a solid domain
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
		if (pbd)  // it is an elastic solid domain
		{
				FEIsotropicElastic* pmi = dynamic_cast<FEIsotropicElastic*>(pbd->GetMaterial());
				if (pmi){  // the actual elastic material is isotropic elastic	

					any_evaluated_domains = true;

					// loop over all the elements	
					for (int iel = 0; iel < pbd->Elements(); ++iel)
					{
						FESolidElement& el = pbd->Element(iel);

						// get minimum edge lenght (squared)    
						double h = get_elem_hmin_alt(mesh, el);

						// loop Gauss points				
						int nint = el.GaussPoints();
						double* gw = el.GaussWeights();
						for (int j = 0; j < nint; ++j)
						{

							// get material parameters
							FEMaterialPoint& mp = *el.GetMaterialPoint(j);
							double mE = pmi->m_E(mp);
							double mv = pmi->m_v(mp);
							double dens = pmi->Density(mp);

							// lame parameters
							double lam = (mv*mE/((1+mv)*(1-2*mv)));
							double mu  = (0.5*mE/(1+mv));

							// wavespeed
							double waveSpd = sqrt((lam+2*mu) / dens);
							if (waveSpd < domMinWaveSpd) domMinWaveSpd = waveSpd;

							// time step
							double dt = 2.0 * h / waveSpd; // 2 -> scheme factor
							if (dt < domMinDT) {
								domMinDT = dt;
								domCritElId = el.GetID();
							}

						} // for gauss point						
					} // for element
					if (domMinDT < globalMinDT) {
						globalMinDT = domMinDT;
						globalCritElId = domCritElId;
					}
					// print domain stats
					feLogEx(pfem, "%-19s %-7.2e  %-7.2e %d\n", 
							mesh.Domain(nd).GetName().c_str(), 
							domMinWaveSpd, domMinDT, domCritElId);
				} else {  // if actual elastic material is isotropic elastic
					feLogErrorEx(pfem, "Wavespeed calculation of other materials than 'isotropic elastic' not unspoorted!");
			} // if explicit elastic material
		} // if elastic domain
	} // for domain
	feLogEx(pfem, "\n");

	// do not adjust if no domains were evaluated
	if (!any_evaluated_domains) {
		feLogEx(pfem, "No supported domains were evaluated, not adjusting dt...\n");
		return true;
	}

	// print global stats
	feLogEx(pfem, "Global stable time step size: %7.2e\n", globalMinDT);
	feLogEx(pfem, "Global critical element id: %d\n", globalCritElId);

	// update critical time step size
	dtcrit = globalMinDT;
	psolver->m_dtcrit = dtcrit;			

	// enforce time step size unless next point is a must point and dt < dtcrit*safety
	if ((pstep->m_dt < safety*dtcrit) && (pstep->m_timeController) && (pstep->m_timeController->m_nmust >= 0)) {
		return true;
	} else {
		pstep->m_dt = safety*dtcrit;
		feLogEx(pfem, "Setting time step size: %7.2e\n", safety*dtcrit);
	}

	return true;
}



//-----------------------------------------------------------------------------
FEExplicitSolidSolver2::FEExplicitSolidSolver2(FEModel* pfem) : 
	FESolver(pfem), m_rigidSolver(pfem),
	m_dofU(pfem), m_dofV(pfem), m_dofSQ(pfem), m_dofRQ(pfem),
	m_dofSU(pfem), m_dofSV(pfem), m_dofSA(pfem)
{
	m_dyn_damping = 0.99;
	m_niter = 0;
	m_nreq = 0;

	m_mass_lumping = HRZ_LUMPING;

	m_adjust_step_size_safety = 0.8;  // default safety factor for limiting time step size
	m_adjust_step_size_stride = 1;  // by default we adjust time step every step
	m_print_stride = 1;  // by default we print to screen every step

	// Allocate degrees of freedom
	DOFS& dofs = pfem->GetDOFS();
	int varD = dofs.AddVariable("displacement", VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable("shell rotation", VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
	int varSU = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT), VAR_VEC3);
	dofs.SetDOFName(varSU, 0, "sx");
	dofs.SetDOFName(varSU, 1, "sy");
	dofs.SetDOFName(varSU, 2, "sz");
	int varSV = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY), VAR_VEC3);
	dofs.SetDOFName(varSV, 0, "svx");
	dofs.SetDOFName(varSV, 1, "svy");
	dofs.SetDOFName(varSV, 2, "svz");
	int varSA = dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION), VAR_VEC3);
	dofs.SetDOFName(varSA, 0, "sax");
	dofs.SetDOFName(varSA, 1, "say");
	dofs.SetDOFName(varSA, 2, "saz");

	// get the DOF indices
	m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCTIY));
	m_dofSQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION));
	m_dofRQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
	m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	m_dofSV.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY));
	m_dofSA.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION));

	// add time step controller
	pfem->AddCallback(time_step_limiter_cb, CB_UPDATE_TIME, this);

}

//-----------------------------------------------------------------------------
void FEExplicitSolidSolver2::Clean()
{
	// if we were blocking logging, unblock	
	GetFEModel()->UnBlockLog();
}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver2::CalculateMassMatrix()
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();

	vector<double> dummy(m_Mi);
	FEGlobalVector Mi(fem, m_Mi, dummy);
	matrix me;
	vector <int> lm;
	vector <double> el_lumped_mass;


	// loop over all domains
	if (m_mass_lumping == NO_MASS_LUMPING)
	{
		// use consistent mass matrix.
		// TODO: implement this
		assert(false);
		return false;
	}
	else if (m_mass_lumping == ROW_SUM_LUMPING)
	{

		feLog("\n");
		feLog("Domain inertia\n");
		feLog("===========================\n");
		feLog("Domain              Mass   \n");
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			// check whether it is a solid domain
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
			if (pbd && !dynamic_cast<FERigidMaterial*>(pbd->GetMaterial()))  // it is an elastic solid domain			
			{
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());

				// nodal coordinates
				vec3d r0[FEElement::MAX_NODES];

				// loop over all the elements
				double Md = 0.0;  // domain total mass
				for (int iel = 0; iel < pbd->Elements(); ++iel)
				{
					FESolidElement& el = pbd->Element(iel);
					pbd->UnpackLM(el, lm);

					int nint = el.GaussPoints();
					int neln = el.Nodes();

					me.resize(neln, neln);
					me.zero();

					// create the element mass matrix
					for (int n = 0; n < nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];

						double* H = el.H(n);
						for (int i = 0; i < neln; ++i)
							for (int j = 0; j < neln; ++j)
							{
								double kab = H[i] * H[j] * detJ0*d;
								me[i][j] += kab;
							}
					}

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						for (int j = 0; j < neln; ++j)
						{
							double kab = me[i][j];
							el_lumped_mass[3 * i    ] += kab;
							el_lumped_mass[3 * i + 1] += kab;
							el_lumped_mass[3 * i + 2] += kab;
						}
						Md += el_lumped_mass[3 * i];
					}

					// assemble element matrix into inv_mass vector 
					Mi.Assemble(el.m_node, lm, el_lumped_mass);

					// hjs: account for mass of rigid nodes
					if (fem.RigidBodies() > 0)
					{
						for (int i = 0; i < neln; ++i)
						{
							FENode& node = mesh.Node(el.m_node[i]);
							if (node.m_rid >= 0)
							{
									FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
																								
									// add mass
									double node_lumped_mass = el_lumped_mass[3 * i];
									RB.m_mass += node_lumped_mass;

									// add moment of inertia (symmetric tensor)
									vec3d dr = node.m_r0 - RB.m_r0;
									RB.m_moi.xx() += node_lumped_mass * (dr.y*dr.y + dr.z*dr.z);
									RB.m_moi.yy() += node_lumped_mass * (dr.x*dr.x + dr.z*dr.z);
									RB.m_moi.zz() += node_lumped_mass * (dr.x*dr.x + dr.y*dr.y);								
									RB.m_moi.xy() -= node_lumped_mass * dr.x*dr.y;
									RB.m_moi.yz() -= node_lumped_mass * dr.y*dr.z;
									RB.m_moi.xz() -= node_lumped_mass * dr.x*dr.z;
							}
						}
					}

				} // loop over elements
			
				// log domain mass
				feLog("%-19s %-7.2e\n", mesh.Domain(nd).GetName().c_str(), Md);
			
			}
			else if (dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd)))
			{
				FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(psd->GetMaterial());
				
				if (!dynamic_cast<FERigidMaterial*>(psd->GetMaterial()))  // if non-rigid material
				{

				// loop over all the elements
				for (int iel = 0; iel < psd->Elements(); ++iel)
				{
					FEShellElement& el = psd->Element(iel);
					psd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					int ndof = 6 * el.Nodes();
					ke.resize(ndof, ndof);
					ke.zero();

					// calculate inertial stiffness
					psd->ElementMassMatrix(el, ke, 1.0);

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(ndof, 0.0);
					for (int i = 0; i < ndof; ++i)
					{
						for (int j = 0; j < ndof; ++j)
						{
							double kab = ke[i][j];
							el_lumped_mass[i] += kab;
						}
					}

					// assemble element matrix into inv_mass vector 
					Mi.Assemble(el.m_node, lm, el_lumped_mass);

					// hjs: account for mass of rigid nodes
					if (fem.RigidBodies() > 0)
					{
						for (int i = 0; i < neln; ++i)
						{
							FENode& node = mesh.Node(el.m_node[i]);
							if (node.m_rid >= 0)
							{
									FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
																								
									// add mass
									double node_lumped_mass = el_lumped_mass[i];
									RB.m_mass += node_lumped_mass;

									// add moment of inertia (symmetric tensor)
									vec3d dr = node.m_r0 - RB.m_r0;
									RB.m_moi.xx() += node_lumped_mass * (dr.y*dr.y + dr.z*dr.z);
									RB.m_moi.yy() += node_lumped_mass * (dr.x*dr.x + dr.z*dr.z);
									RB.m_moi.zz() += node_lumped_mass * (dr.x*dr.x + dr.y*dr.y);								
									RB.m_moi.xy() -= node_lumped_mass * dr.x*dr.y;
									RB.m_moi.yz() -= node_lumped_mass * dr.y*dr.z;
									RB.m_moi.xz() -= node_lumped_mass * dr.x*dr.z;
							}
						}
					}

				}
				}  // if non-rigid
			}
			else
			{
//				return false;
			}
		}
	}
	else if (m_mass_lumping == HRZ_LUMPING)
	{
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			// check whether it is a solid domain
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
			if (pbd && !dynamic_cast<FERigidMaterial*>(pbd->GetMaterial()))  // it is an elastic solid domain
			{
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());

				// loop over all the elements
				double Md = 0.0;  // domain total mass
				for (int iel = 0; iel < pbd->Elements(); ++iel)
				{
					FESolidElement& el = pbd->Element(iel);
					pbd->UnpackLM(el, lm);

					int nint = el.GaussPoints();
					int neln = el.Nodes();

					me.resize(neln, neln);
					me.zero();

					// calculate the element mass matrix (and element mass).
					double Me = 0.0;
					double* w = el.GaussWeights();
					for (int n = 0; n < nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = pbd->detJ0(el, n)*el.GaussWeights()[n];
						Me += d * detJ0 * w[n];

						double* H = el.H(n);
						for (int i = 0; i < neln; ++i)
							for (int j = 0; j < neln; ++j)
							{
								double kab = H[i] * H[j] * detJ0*d;
								me[i][j] += kab;
							}
					}
					Md += Me;

					// calculate sum of diagonals
					double S = 0.0;
					for (int i = 0; i < neln; ++i) S += me[i][i];

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(3 * neln, 0.0);
					for (int i = 0; i < neln; ++i)
					{
						double mab = me[i][i] * Me / S;
						el_lumped_mass[3 * i    ] = mab;
						el_lumped_mass[3 * i + 1] = mab;
						el_lumped_mass[3 * i + 2] = mab;
					}

					// assemble element matrix into inv_mass vector 
					Mi.Assemble(el.m_node, lm, el_lumped_mass);

					// hjs: account for mass of rigid nodes
					if (fem.RigidBodies() > 0)
					{
						for (int i = 0; i < neln; ++i)
						{
							FENode& node = mesh.Node(el.m_node[i]);
							if (node.m_rid >= 0)
							{
									FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
																								
									// add mass
									double node_lumped_mass = el_lumped_mass[3 * i];
									RB.m_mass += node_lumped_mass;

									// add moment of inertia (symmetric tensor)
									vec3d dr = node.m_r0 - RB.m_r0;
									RB.m_moi.xx() += node_lumped_mass * (dr.y*dr.y + dr.z*dr.z);
									RB.m_moi.yy() += node_lumped_mass * (dr.x*dr.x + dr.z*dr.z);
									RB.m_moi.zz() += node_lumped_mass * (dr.x*dr.x + dr.y*dr.y);								
									RB.m_moi.xy() -= node_lumped_mass * dr.x*dr.y;
									RB.m_moi.yz() -= node_lumped_mass * dr.y*dr.z;
									RB.m_moi.xz() -= node_lumped_mass * dr.x*dr.z;
							}
						}
					}

				} // loop over elements

				// log domain mass
				feLog("%-19s %-7.2e\n", mesh.Domain(nd).GetName().c_str(), Md);

			}
			else if(dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd)))
			{
				FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(psd->GetMaterial());
				if (!dynamic_cast<FERigidMaterial*>(psd->GetMaterial()))  // if non-rigid material
				{

				// loop over all the elements
				double Md = 0.0;  // domain total mass
				for (int iel = 0; iel < psd->Elements(); ++iel)
				{
					FEShellElement& el = psd->Element(iel);
					psd->UnpackLM(el, lm);

					// create the element's stiffness matrix
					FEElementMatrix ke(el);
					int neln = el.Nodes();
					int ndof = 6 * el.Nodes();
					ke.resize(ndof, ndof);
					ke.zero();

					// calculate inertial stiffness
					psd->ElementMassMatrix(el, ke, 1.0);

					// calculate the element mass
					double Me = 0.0;
					double* w = el.GaussWeights();
					for (int n = 0; n < el.GaussPoints(); ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						double d = pme->Density(mp);
						double detJ0 = psd->detJ0(el, n) * el.GaussWeights()[n];
						Me += d * detJ0 * w[n];
					}
					Md += Me;

					// calculate sum of diagonals
					double S = 0.0;
					for (int i = 0; i < ndof; ++i) S += ke[i][i] / 3.0;

					// reduce to a lumped mass vector and add up the total
					el_lumped_mass.assign(ndof, 0.0);
					for (int i = 0; i < ndof; ++i)
					{
						double mab = ke[i][i] * Me / S;
						el_lumped_mass[i] = mab;
					}

					// assemble element matrix into inv_mass vector 
					Mi.Assemble(el.m_node, lm, el_lumped_mass);

					// hjs: account for mass of rigid nodes
					if (fem.RigidBodies() > 0)
					{
						for (int i = 0; i < neln; ++i)
						{
							FENode& node = mesh.Node(el.m_node[i]);
							if (node.m_rid >= 0)
							{
									FERigidBody& RB = *fem.GetRigidBody(node.m_rid);
																								
									// add mass
									double node_lumped_mass = el_lumped_mass[i];
									RB.m_mass += node_lumped_mass;

									// add moment of inertia (symmetric tensor)
									vec3d dr = node.m_r0 - RB.m_r0;
									RB.m_moi.xx() += node_lumped_mass * (dr.y*dr.y + dr.z*dr.z);
									RB.m_moi.yy() += node_lumped_mass * (dr.x*dr.x + dr.z*dr.z);
									RB.m_moi.zz() += node_lumped_mass * (dr.x*dr.x + dr.y*dr.y);								
									RB.m_moi.xy() -= node_lumped_mass * dr.x*dr.y;
									RB.m_moi.yz() -= node_lumped_mass * dr.y*dr.z;
									RB.m_moi.xz() -= node_lumped_mass * dr.x*dr.z;
							}
						}
					}

				}

				// log domain mass
				feLog("%-19s %-7.2e\n", mesh.Domain(nd).GetName().c_str(), Md);
				
				}  // if non-rigid material
			}
			else
			{
				// TODO: we can only do solid domains right now.
				return false;
			}
		}
	}
	else
	{
		assert(false);
		return false;
	}

	// we need the inverse of the lumped masses later
	// Also, make sure the lumped masses are positive.
	for (int i = 0; i < m_Mi.size(); ++i)
	{
//		if (m_Mi[i] <= 0.0) return false;
		if (m_Mi[i] > 1.0e-12) m_Mi[i] = 1.0 / m_Mi[i];
	}

	// also log rigid body inertia now that we have added that of rigified nodes
	const int NRB = fem.RigidBodies();
	if (NRB > 0){
		feLog("\n");
		feLog("Rigid body inertia\n");
		feLog("================================================================================\n");
		feLog("Domain              Mass     MoI_xx   MoI_yy   MoI_zz   MoI_xy   MoI_yz   MoI_xz\n");
		for (int i=0; i<NRB; ++i)
		{
			FERigidBody& RB = *fem.GetRigidBody(i);
			feLog("%-19s %-7.2e %-7.2e %-7.2e %-7.2e %-7.2e %-7.2e %-7.2e\n", 
			RB.GetName().c_str(), RB.m_mass, 
			RB.m_moi.xx(), RB.m_moi.yy(), RB.m_moi.zz(),
			RB.m_moi.xy(), RB.m_moi.yz(), RB.m_moi.xz());
		}
	}
	feLog("\n");

	return true;
}

//-----------------------------------------------------------------------------
//! Determine the number of linear equations and assign equation numbers
//!	This function initializes the equation system.
//! It is assumed that all free dofs up until now have been given an ID >= 0
//! and the fixed or rigid dofs an ID < 0.
//! After this operation the nodal ID array will contain the equation
//! number assigned to the corresponding degree of freedom. To distinguish
//! between free or unconstrained dofs and constrained ones the following rules
//! apply to the ID array:
//!
//!           /
//!          |  >=  0 --> dof j of node i is a free dof
//! ID[i][j] <  == -1 --> dof j of node i is a fixed (no equation assigned too) (hjs: no reactions stored)
//!          |  <  -1 --> dof j of node i is constrained and has equation nr = -ID[i][j]-2
//!           \
//!

bool FEExplicitSolidSolver2::InitEquations()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FESolver::InitEquations() == false) return false;

	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq); // hjs: avoid dependency of rigidSolver (only used here)
	if (neq == -1) return false; 
	else m_neq = neq;

	feLog("number of regular equations %d, number of rigid equations %d\n", m_nreq, m_neq-m_nreq);

	return true;

}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver2::Init()
{
	if (FESolver::Init() == false) return false;

	// get nr of equations including rigids
	int neq = m_neq;

	// allocate vectors
	m_Fn.assign(neq, 0);
	m_Fr.assign(neq, 0);
	m_ui.assign(neq, 0);
	m_Ut.assign(neq, 0);
	m_Mi.assign(neq, 0.0); // hjs: use nreq

	GetFEModel()->Update();

	// we need to fill the total displacement vector m_Ut
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofSQ[0]);
	gather(m_Ut, mesh, m_dofSQ[1]);
	gather(m_Ut, mesh, m_dofSQ[2]);
	gather(m_Ut, mesh, m_dofSU[0]);
	gather(m_Ut, mesh, m_dofSU[1]);
	gather(m_Ut, mesh, m_dofSU[2]);

	// calculate the inverse mass vector for the explicit analysis
	if (CalculateMassMatrix() == false)
	{
		feLogError("Failed building mass matrix.");
		return false;
	}

	// Calculate initial residual to be used on the first time step
	if (Residual(m_R0) == false) return false;

	// calculate the initial acceleration
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) node.m_at.x = m_R0[n] * m_Mi[n];
		if ((n = node.m_ID[m_dofU[1]]) >= 0) node.m_at.y = m_R0[n] * m_Mi[n];
		if ((n = node.m_ID[m_dofU[2]]) >= 0) node.m_at.z = m_R0[n] * m_Mi[n];

		if ((n = node.m_ID[m_dofSU[0]]) >= 0) node.set(m_dofSA[0], m_R0[n] * m_Mi[n]);
		if ((n = node.m_ID[m_dofSU[1]]) >= 0) node.set(m_dofSA[1], m_R0[n] * m_Mi[n]);
		if ((n = node.m_ID[m_dofSU[2]]) >= 0) node.set(m_dofSA[2], m_R0[n] * m_Mi[n]);
	}

	// set the dynamic update flag only if we are running a dynamic analysis
	bool b = (fem.GetCurrentStep()->m_nanalysis == FE_DYNAMIC ? true : false);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEElasticSolidDomain* d = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
		FEElasticShellDomain* s = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(i));
		if (d) d->SetDynamicUpdateFlag(b);
		if (s) s->SetDynamicUpdateFlag(b);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEExplicitSolidSolver2::Update(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// update kinematics
	UpdateKinematics(ui);

	// update element stresses
	fem.Update();

}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEExplicitSolidSolver2::UpdateKinematics(vector<double>& ui)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	UpdateRigidBodies(ui);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofU[0]);
	scatter(U, mesh, m_dofU[1]);
	scatter(U, mesh, m_dofU[2]);
	// rotational dofs
	scatter(U, mesh, m_dofSQ[0]);
	scatter(U, mesh, m_dofSQ[1]);
	scatter(U, mesh, m_dofSQ[2]);
	// shell displacement
	scatter(U, mesh, m_dofSU[0]);
	scatter(U, mesh, m_dofSU[1]);
	scatter(U, mesh, m_dofSU[2]);

	// make sure the prescribed displacements are fullfilled
	int ndis = fem.BoundaryConditions();
	for (int i=0; i<ndis; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
	}

	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints() > 0)
	{
		LCM.Update();
	}

	// Update the spatial nodal positions
	// Don't update rigid nodes since they are already updated
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid == -1)
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
	}
}

//-----------------------------------------------------------------------------
//! Updates the rigid body data
void FEExplicitSolidSolver2::UpdateRigidBodies(vector<double>& ui)
{
	// get the number of rigid bodies
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	const int NRB = fem.RigidBodies();

	// first calculate the rigid body displacement increments
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		int *lm = RB.m_LM;
		double* du = RB.m_du;

		if (RB.m_prb == 0)
		{
			for (int j=0; j<6; ++j)
			{
				du[j] = (lm[j] >=0 ? ui[lm[j]] : 0);
			}
		}
	}

	// for prescribed displacements, the displacement increments are evaluated differently
	// TODO: Is this really necessary? Why can't the ui vector contain the correct values?
	const int NRD = fem.RigidPrescribedBCs();
	for (int i=0; i<NRD; ++i)
	{
		FERigidBodyDisplacement& dc = *fem.GetRigidPrescribedBC(i);
		if (dc.IsActive())
		{
			FERigidBody& RB = *fem.GetRigidBody(dc.GetID());
			int I = dc.GetBC();
			RB.m_du[I] = dc.Value() - RB.m_Up[I];
		}
	}

	// update the rigid bodies
	for (int i=0; i<NRB; ++i)
	{
		// get the rigid body
		FERigidBody& RB = *fem.GetRigidBody(i);
		double* du = RB.m_du;

		// This is the "old" update algorithm which has some issues. It does not produce the correct
		// rigid body orientation when the rotational degrees of freedom are prescribed.
		RB.m_rt.x = RB.m_rp.x + du[0];
		RB.m_rt.y = RB.m_rp.y + du[1];
		RB.m_rt.z = RB.m_rp.z + du[2];

		vec3d r = vec3d(du[3], du[4], du[5]);
		double w = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
		quatd dq = quatd(w, r);

		quatd Q = dq*RB.m_qp;
		Q.MakeUnit();
		RB.SetRotation(Q);

		if (RB.m_prb) du = RB.m_dul;
		RB.m_Ut[0] = RB.m_Up[0] + du[0];
		RB.m_Ut[1] = RB.m_Up[1] + du[1];
		RB.m_Ut[2] = RB.m_Up[2] + du[2];
		RB.m_Ut[3] = RB.m_Up[3] + du[3];
		RB.m_Ut[4] = RB.m_Up[4] + du[4];
		RB.m_Ut[5] = RB.m_Up[5] + du[5];
	}

	// we need to update the position of rigid nodes
	fem.UpdateRigidMesh();

	// Since the rigid nodes are repositioned we need to update the displacement DOFS
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		if (node.m_rid >= 0)
		{
			vec3d ut = node.m_rt - node.m_r0;
			node.set_vec3d(m_dofU[0], m_dofU[1], m_dofU[2], ut);
		}
	}
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEExplicitSolidSolver2::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);
	ar & m_nrhs & m_niter & m_nref & m_ntotref & m_naug & m_neq & m_nreq;
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the DoSolve routine 
//!  and deals with exceptions that require the immediate termination of
//!	 the solution eg negative Jacobians.
bool FEExplicitSolidSolver2::SolveStep()
{
	bool bret;

	try
	{
		// let's try to solve the step
		bret = DoSolve();
	}
	catch (NegativeJacobian e)
	{
		// A negative jacobian was detected
		feLogError("Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
		return false;
	}
	catch (MaxStiffnessReformations) // shouldn't happen for an explicit analysis!
	{
		// max nr of reformations is reached
		feLogError("Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		feLogWarning("User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		feLogWarning("User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize) // shouldn't happen for an explicit analysis!
	{
		// a zero line step size was detected
		feLogError("Zero line step size.");
		return false;
	}
	catch (EnergyDiverging) // shouldn't happen for an explicit analysis!
	{
		// problem was diverging after stiffness reformation
		feLogError("Problem diverging uncontrollably.");
		return false;
	}
	catch (FEMultiScaleException)
	{
		// the RVE problem didn't solve
		feLogError("The RVE problem has failed. Aborting macro run.");
		return false;
	}

	return bret;
}


//-----------------------------------------------------------------------------
//! Prepares the data for the time step. 
void FEExplicitSolidSolver2::PrepStep()
{
	int i, j;

	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs  = 0;	// nr of RHS evaluations
	m_nref  = 0;	// nr of stiffness reformations
	m_ntotref = 0;
	m_naug  = 0;	// nr of augmentations

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_ap = ni.m_at;
		ni.UpdateValues();
	}

	const FETimeInfo& tp = fem.GetTime();

	// apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
	vector<double> dummy(m_neq, 0.0);
	zero(m_Fn);
	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
	NodalLoads(Fn, tp); // warning, will ignore nodal loads applied to rigid body nodes

	// apply prescribed displacements
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	int neq = m_neq;
	int nbc = fem.BoundaryConditions();
	for (i=0; i<nbc; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.PrepStep(ui);
	}

	// initialize rigid bodies
	int NO = fem.RigidBodies();
	for (i=0; i<NO; ++i) fem.GetRigidBody(i)->Init();

	// calculate local rigid displacements
	for (i=0; i<(int) fem.RigidPrescribedBCs(); ++i)
	{
		FERigidBodyDisplacement& DC = *fem.GetRigidPrescribedBC(i);
		FERigidBody& RB = *fem.GetRigidBody(DC.GetID());
		if (DC.IsActive())
		{
			int I = DC.GetBC();
			RB.m_dul[I] = DC.Value() - RB.m_Ut[I];
		}
	}

	// calculate global rigid displacements
	for (i=0; i<NO; ++i)
	{
		FERigidBody* prb = fem.GetRigidBody(i);
		if (prb)
		{
			FERigidBody& RB = *prb;
			if (RB.m_prb == 0)
			{
				for (j=0; j<6; ++j) RB.m_du[j] = RB.m_dul[j];
			}
			else
			{
				double* dul = RB.m_dul;
				vec3d dr = vec3d(dul[0], dul[1], dul[2]);
				
				vec3d v = vec3d(dul[3], dul[4], dul[5]);
				double w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
				quatd dq = quatd(w, v);

				FERigidBody* pprb = RB.m_prb;

				vec3d r0 = RB.m_rt;
				quatd Q0 = RB.GetRotation();

				dr = Q0*dr;
				dq = Q0*dq*Q0.Inverse();

				while (pprb)
				{
					vec3d r1 = pprb->m_rt;
					dul = pprb->m_dul;

					quatd Q1 = pprb->GetRotation();
					
					dr = r0 + dr - r1;

					// grab the parent's local displacements
					vec3d dR = vec3d(dul[0], dul[1], dul[2]);
					v = vec3d(dul[3], dul[4], dul[5]);
					w = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
					quatd dQ = quatd(w, v);

					dQ = Q1*dQ*Q1.Inverse();

					// update global displacements
					quatd Qi = Q1.Inverse();
					dr = dR + r1 + dQ*dr - r0;
					dq = dQ*dq;

					// move up in the chain
					pprb = pprb->m_prb;
					Q0 = Q1;
				}

				// set global displacements
				double* du = RB.m_du;

				du[0] = dr.x;
				du[1] = dr.y;
				du[2] = dr.z;

				v = dq.GetVector();
				w = dq.GetAngle();
				du[3] = w*v.x;
				du[4] = w*v.y;
				du[5] = w*v.z;
			}
		}
	}

	// store rigid displacements in Ui vector
	for (i=0; i<NO; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		for (j=0; j<6; ++j)
		{
			int I = -RB.m_LM[j]-2;
			if (I >= 0) ui[I] = RB.m_du[j];
		}
	}

	// intialize material point data
	for (i=0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	fem.Update();

}

//-----------------------------------------------------------------------------
bool FEExplicitSolidSolver2::DoSolve()
{
	// Get the current step
	// FEModel& fem = *GetFEModel();
	// FEAnalysis* pstep = fem.GetCurrentStep();
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());

	// prepare for solve
	PrepStep();

	double dt = fem.GetTime().timeIncrement;

	feLog(" %d\n", m_niter+1);

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes(); // this is the total number of nodes in the mesh


	// collect accelerations, velocities, displacements
	vector<double> an(m_neq, 0.0), vn(m_neq, 0.0), un(m_neq, 0.0);
//#pragma omp parallel for
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d vt = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) { un[n] = node.m_rt.x - node.m_r0.x; vn[n] = vt.x; an[n] = node.m_at.x; }
		if ((n = node.m_ID[m_dofU[1]]) >= 0) { un[n] = node.m_rt.y - node.m_r0.y; vn[n] = vt.y; an[n] = node.m_at.y; }
		if ((n = node.m_ID[m_dofU[2]]) >= 0) { un[n] = node.m_rt.z - node.m_r0.z; vn[n] = vt.z; an[n] = node.m_at.z; }

		if ((n = node.m_ID[m_dofSU[0]]) >= 0) { un[n] = node.get(m_dofSU[0]); vn[n] = node.get(m_dofSV[0]); an[n] = node.get(m_dofSA[0]); }
		if ((n = node.m_ID[m_dofSU[1]]) >= 0) { un[n] = node.get(m_dofSU[1]); vn[n] = node.get(m_dofSV[1]); an[n] = node.get(m_dofSA[1]); }
		if ((n = node.m_ID[m_dofSU[2]]) >= 0) { un[n] = node.get(m_dofSU[2]); vn[n] = node.get(m_dofSV[2]); an[n] = node.get(m_dofSA[2]); }
	}
	for (int i=0; i<fem.RigidBodies(); ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		
		int n;
		
		// TODO: does RB.m_du hold the correct values for this?
		// linear acceleration and velocity of center of mass
		// if ((n = RB.m_LM[0]) >= 0) { un[n] = RB.m_rt.x - RB.m_r0.x; vn[n] = RB.m_vt.x; an[n] = RB.m_at.x; }
		// if ((n = RB.m_LM[1]) >= 0) { un[n] = RB.m_rt.y - RB.m_r0.y; vn[n] = RB.m_vt.y; an[n] = RB.m_at.y; }
		// if ((n = RB.m_LM[2]) >= 0) { un[n] = RB.m_rt.z - RB.m_r0.z; vn[n] = RB.m_vt.z; an[n] = RB.m_at.z; }		
		if ((n = RB.m_LM[0]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_vt.x; an[n] = RB.m_at.x; }
		if ((n = RB.m_LM[1]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_vt.y; an[n] = RB.m_at.y; }
		if ((n = RB.m_LM[2]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_vt.z; an[n] = RB.m_at.z; }
		
		// angular acceleration and velocity of center of mass
		if ((n = RB.m_LM[3]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_wt.x; an[n] = RB.m_alt.x; }
		if ((n = RB.m_LM[4]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_wt.y; an[n] = RB.m_alt.y; }
		if ((n = RB.m_LM[5]) >= 0) { un[n] = RB.m_du[n]; vn[n] = RB.m_wt.z; an[n] = RB.m_alt.z; }
	}

	// predictor
	vector<double> v_pred(m_neq, 0.0);
//#pragma omp parallel for
	for (int i = 0; i < m_neq; ++i)
	{
		v_pred[i] = vn[i] + an[i] * dt*0.5; // update velocity
		m_ui[i] = dt * v_pred[i]; // update displacements
	}
	Update(m_ui);

	// evaluate acceleration
	Residual(m_R1);
	vector<double> anp1(m_neq);
//#pragma omp parallel for
	for (int i = 0; i < m_nreq; ++i) // regular equations only
	{
		anp1[i] = m_R1[i] * m_Mi[i];
	}

	// update velocity
	vector<double> vnp1(m_neq, 0.0);
//#pragma omp parallel for
	double ke = 0.0; // kinetic energy
	for (int i = 0; i < m_nreq; ++i)  // regular equations only
	{
		vnp1[i] = v_pred[i] + anp1[i]*dt*0.5;
		ke += 0.5 * vnp1[i]*vnp1[i] / m_Mi[i];
	}

	// TODO hjs: is there other non-regular equations than the rigid ones we are missing??

	// handle rigid bodies, (if not constrained)
	for (int i=0; i<fem.RigidBodies(); ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);

		int n; 

		// mass matrix
		double M = RB.m_mass;
		double M_inv = 1.0 / M; 

		// compute translational acceleration and store in vector, update velocity add contribute to kinetic energy
		if ((n = RB.m_LM[0]) >= 0) {
			anp1[n] = m_R1[n] * M_inv; 
			vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; 
			ke += 0.5 * vnp1[n]*vnp1[n] * M;
		}
		if ((n = RB.m_LM[1]) >= 0) {
			anp1[n] = m_R1[n] * M_inv; 
			vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; 
			ke += 0.5 * vnp1[n]*vnp1[n] * M;
		}
		if ((n = RB.m_LM[2]) >= 0) {
			anp1[n] = m_R1[n] * M_inv; 
			vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; 
			ke += 0.5 * vnp1[n]*vnp1[n] * M;
		}

		// If all rotational dofs are constrained, do not try to compute angular acceleration
		if ( RB.m_LM[3] < 0 &&  RB.m_LM[4] < 0 && RB.m_LM[5] < 0) continue;
		// TODO: what to do if rotation is partly constrained??
		if ( RB.m_LM[3] < 0 ||  RB.m_LM[4] < 0 || RB.m_LM[5] < 0 ){
			feLogError("Explicit solver cannot handle partly constrained rotational dofs of rigid body %d\n", i);
			assert(false);
		}
		
		// evaluate mass moment of inertia at t
		mat3d Rt = RB.GetRotation().RotationMatrix();
		mat3ds Jt = (Rt*RB.m_moi*Rt.transpose()).sym();
		mat3ds Jt_inv = Jt.inverse();

		// get torque
		vec3d RB_torque = vec3d(0.0, 0.0, 0.0);
		n = RB.m_LM[3]; RB_torque.x = m_R1[n];
		n = RB.m_LM[4]; RB_torque.y = m_R1[n];
		n = RB.m_LM[5]; RB_torque.z = m_R1[n];

		// compute angular acceleration and store in vector and update velocity
		vec3d RB_al = Jt_inv * RB_torque;	
		vec3d RB_w = vec3d(0.0, 0.0, 0.0);	// angular velocity
		n = RB.m_LM[3]; anp1[n] = RB_al.x; vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; RB_w.x = vnp1[n];
		n = RB.m_LM[4]; anp1[n] = RB_al.y; vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; RB_w.y = vnp1[n];
		n = RB.m_LM[5]; anp1[n] = RB_al.z; vnp1[n] = v_pred[n] + anp1[n]*dt*0.5; RB_w.z = vnp1[n];

		// kinetic energy
		ke += 0.5 * (RB_w * (Jt * RB_w)); 

	}

	// increase iteration number
	m_niter++;

	// scatter velocity and accelerations
//#pragma omp parallel for
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) { node.set(m_dofV[0], vnp1[n]); node.m_at.x = anp1[n]; }
		if ((n = node.m_ID[m_dofU[1]]) >= 0) { node.set(m_dofV[1], vnp1[n]); node.m_at.y = anp1[n]; }
		if ((n = node.m_ID[m_dofU[2]]) >= 0) { node.set(m_dofV[2], vnp1[n]); node.m_at.z = anp1[n]; }

		if ((n = node.m_ID[m_dofSU[0]]) >= 0) { node.set(m_dofSV[0], vnp1[n]); node.set(m_dofSA[0], anp1[n]); }
		if ((n = node.m_ID[m_dofSU[1]]) >= 0) { node.set(m_dofSV[1], vnp1[n]); node.set(m_dofSA[1], anp1[n]); }
		if ((n = node.m_ID[m_dofSU[2]]) >= 0) { node.set(m_dofSV[2], vnp1[n]); node.set(m_dofSA[2], anp1[n]); }
	}
	for (int i=0; i<fem.RigidBodies(); ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		
		int n;
		
		// linear acceleration and velocity of center of mass
		if ((n = RB.m_LM[0]) >= 0) { RB.m_vt.x = vnp1[n]; RB.m_at.x = anp1[n]; }
		if ((n = RB.m_LM[1]) >= 0) { RB.m_vt.y = vnp1[n]; RB.m_at.y = anp1[n]; }
		if ((n = RB.m_LM[2]) >= 0) { RB.m_vt.z = vnp1[n]; RB.m_at.z = anp1[n]; }

		// angular acceleration and velocity of center of mass
		if ((n = RB.m_LM[3]) >= 0) { RB.m_wt.x = vnp1[n]; RB.m_alt.x = anp1[n]; }
		if ((n = RB.m_LM[4]) >= 0) { RB.m_wt.y = vnp1[n]; RB.m_alt.y = anp1[n]; }
		if ((n = RB.m_LM[5]) >= 0) { RB.m_wt.z = vnp1[n]; RB.m_alt.z = anp1[n]; }

	}

	// do minor iterations callbacks
	fem.DoCallback(CB_MINOR_ITERS);

	// update the total displacements
	m_Ut += m_ui;

	m_R0 = m_R1;

	// print energies
	feLog("Total kinetic energy: %7.2e\n", ke);

	double se = 0.0;  // total strain energy
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{		
		FEExplicitSolidDomain* pbd = dynamic_cast<FEExplicitSolidDomain*>(&mesh.Domain(nd));
		if (pbd) se += pbd->StrainEnergy();
	}
	feLog("Total strain energy: %7.2e\n", se);


	return true;
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEExplicitSolidSolver2::Residual(vector<double>& R)
{
	// get the time information
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	const FETimeInfo& tp = fem.GetTime();

	// initialize residual with concentrated nodal loads
	R = m_Fn;

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	// FEGlobalVector RHS(fem, R, m_Fr);  // hjs: does not sum rigid body forces
	FEResidualVector RHS(fem, R, m_Fr);  // hjs: inherites from FEGlobalVector, but sums rigid body forces during assembly

	// zero rigid body reaction forces
	int NRB = fem.RigidBodies();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = *fem.GetRigidBody(i);
		RB.m_Fr = RB.m_Mr = vec3d(0,0,0);
	}

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// calculate the internal (stress) forces
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
		dom.InternalForces(RHS);
	}

	// calculate the body forces
	for (int j = 0; j<fem.BodyLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
		if (pbf && pbf->IsActive())
		{
			for (int i = 0; i<pbf->Domains(); ++i)
			{
				FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(*pbf->Domain(i));
				if (pbf) dom.BodyForce(RHS, *pbf);
			}
		}
	}

	// calculate forces due to surface loads
	int nsl = fem.SurfaceLoads();
	for (int i=0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->LoadVector(RHS, tp);
	}

	// calculate contact forces
	if (fem.SurfacePairConstraints() > 0)
	{
		ContactForces(RHS);
	}

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// forces due to point constraints
//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->LoadVector(this, R);

	// set the nodal reaction forces 
	// TODO: Is this a good place to do this?
	int count_rnodes = 0;
	vec3d f_sum = vec3d(0.0, 0.0, 0.0);
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.set_load(m_dofU[0], 0);
		node.set_load(m_dofU[1], 0);
		node.set_load(m_dofU[2], 0);
		node.set_load(m_dofSU[0], 0);
		node.set_load(m_dofSU[1], 0);
		node.set_load(m_dofSU[2], 0);

		int n;
		if ((n = -node.m_ID[m_dofU[0]] - 2) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[1]] - 2) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[2]] - 2) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);

		if ((n = -node.m_ID[m_dofSU[0]] - 2) >= 0) node.set_load(m_dofSU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofSU[1]] - 2) >= 0) node.set_load(m_dofSU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofSU[2]] - 2) >= 0) node.set_load(m_dofSU[2], -m_Fr[n]);

	}

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FEExplicitSolidSolver2::ContactForces(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) {
			pci->LoadVector(R, tp);
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FEExplicitSolidSolver2::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->LoadVector(R, tp);
	}
}