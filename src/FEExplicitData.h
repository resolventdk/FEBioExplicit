#pragma once
#include <FECore/NodeDataRecord.h>
#include <FECore/ElementDataRecord.h>

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

//=============================================================================
// E L E M E N T  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogElemVolumetricStrainRate: public FELogElemData
{ 
public: 
	FELogElemVolumetricStrainRate(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemVolume: public FELogElemData
{ 
public: 
	FELogElemVolume(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};