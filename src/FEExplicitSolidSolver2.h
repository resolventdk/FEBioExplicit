#pragma once
#include "FEBioMech/FERigidSolver.h"
#include "FECore/FESolver.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/FETimeInfo.h"
#include "FECore/FEDofList.h"

//-----------------------------------------------------------------------------
//! This class implements a nonlinear explicit solver for solid mechanics
//! problems.
class FEExplicitSolidSolver2 : public FESolver
{
	enum MassLumpingMethod
	{
		NO_MASS_LUMPING,	// use consistent mass matrix
		ROW_SUM_LUMPING,	// use simple row-sum lumping
		HRZ_LUMPING			// use Hinton-Rock-Zienkiewicz lumping
	};

public:
	//! constructor
	FEExplicitSolidSolver2(FEModel* pfem);

	//! destructor
	virtual ~FEExplicitSolidSolver2() {}

public:
	//! Data initialization
	bool Init() override;

	bool InitEquations() override;

	//! clean up
	void Clean() override;

	//! Solve an analysis step
	bool SolveStep() override;

	//! Update data
	void Update(vector<double>& ui) override;

	//! Serialize data
	void Serialize(DumpStream& ar) override;

public:
	//! update kinematics
	void UpdateKinematics(vector<double>& ui);

	//! Update rigid bodies 
	void UpdateRigidBodies(vector<double>& ui);

	//! solve the step
	bool DoSolve();

	void PrepStep();

	bool Residual(vector<double>& R);

	void UpdateExplicitMaterialPoints();

	void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

	void ContactForces(FEGlobalVector& R);

private:
	bool CalculateMassMatrix();

	double CalculateStableTimeIncrement();

public:
	int			m_mass_lumping;	//!< specify mass lumping method
	double		m_dyn_damping;	//!< velocity damping for the explicit solver
	int			m_print_stride;	//!< how often to print to screen
	double      m_dtcrit;
	double      m_adjust_step_size_safety;  //!< safety factor to critical est. time step
	int         m_adjust_step_size_stride;  //!< how often to adjust time step size 
	double      m_bv_c1;
	double      m_bv_c2;

public:
	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

	vector<double> m_Mi;	//!< inverse mass vector for explicit analysis
	vector<double> m_Fn;	//!< concentrated nodal force vector
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)

	vector<double> m_ui;	//!< displacement increment vector

	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

protected:
	FEDofList	m_dofU, m_dofV, m_dofSQ, m_dofRQ;
	FEDofList	m_dofSU, m_dofSV, m_dofSA;

protected:
    FERigidSolverNew	m_rigidSolver;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
