#pragma once
#include <FECore/NodeDataRecord.h>

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FENodeForceX: public FENodeLogData
{ 
public: 
	FENodeForceX(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceY: public FENodeLogData
{ 
public: 
	FENodeForceY(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceZ: public FENodeLogData
{ 
public: 
	FENodeForceZ(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};
