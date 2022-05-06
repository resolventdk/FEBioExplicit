#include "FEExplicitData.h"
#include "FEExplicitSolidSolver2.h"
#include "FECore/FESolidElement.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
double FENodeForceX::value(int nnode) 
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FEExplicitSolidSolver2* psolid_solver = dynamic_cast<FEExplicitSolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<double>& Fn = psolid_solver->m_Fn;
		vector<int>& id = mesh.Node(nnode).m_ID;

		double Fx = 0.0;
		if (id[0] >= 0) Fx = Fn[id[0]];
		else if (-id[0] - 2 >= 0) Fx = Fr[-id[0] - 2];
		return Fx;
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceY::value(int nnode) 
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FEExplicitSolidSolver2* psolid_solver = dynamic_cast<FEExplicitSolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<int>& id = mesh.Node(nnode).m_ID;
		return (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceZ::value(int nnode) 
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FEExplicitSolidSolver2* psolid_solver = dynamic_cast<FEExplicitSolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		vector<int>& id = mesh.Node(nnode).m_ID;
		return (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FELogElemVolumetricStrainRate::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();		
		val += pt.RateOfDeformation().tr();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemVolume::value(FEElement& el)
{
	FESolidElement* sel = dynamic_cast<FESolidElement*>(&el);
	if (sel) {
		double val = 0.0;
		int nint = sel->GaussPoints();
		double* w = sel->GaussWeights();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();		
			val += pt.m_J*w[i];
		}
		return val;
	} else {
		return 0.0;
	}
}