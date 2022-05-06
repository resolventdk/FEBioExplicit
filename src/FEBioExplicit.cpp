#include <FECore/sdk.h>
#include "FEExplicitSolidSolver2.h"
#include "FEExplicitData.h"
#include "FEScaledDensityMapGenerator.h"
#include "FEExplicitElastic.h"

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 0;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// create the module
	const char* info = \
		"{ "
		"   \"title\" : \"Tool kit for explicit dynamic structural simulations in FEBio\","
		"   \"info\"  : \"Reworked solver and features for dynamic structural simulations.\","
		"   \"author\": \"Henrik Spietz\","
		"   \"version\": \"0.0\""
        "}";

	fecore.CreateModule("explicit-solid2");
	fecore.SetModuleDependency("solid");
	
	// Solver classes
	REGISTER_FECORE_CLASS(FEExplicitSolidSolver2, "explicit-solid2");

	// Materials
	REGISTER_FECORE_CLASS(FEExplicitElastic, "explicit elastic");

	// Derived from FENodeLogData
	REGISTER_FECORE_CLASS(FENodeForceX, "Rx2");
	REGISTER_FECORE_CLASS(FENodeForceY, "Ry2");
	REGISTER_FECORE_CLASS(FENodeForceZ, "Rz2");

	// Derived from FENodeLogData
	REGISTER_FECORE_CLASS(FELogElemVolume, "V");
	REGISTER_FECORE_CLASS(FELogElemVolumetricStrainRate, "dVdt");

	// Derived from FEDataGenerator
	REGISTER_FECORE_CLASS(FEScaledDensityMapGenerator, "scaled density");

	fecore.SetActiveModule(0);
}
