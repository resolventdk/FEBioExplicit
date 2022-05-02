#include "FEExplicitUtils.h"
#include "FEScaledDensityMapGenerator.h"

#include "FECore/FEModel.h"
#include "FECore/FEDomainMap.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEScaledDensityMapGenerator, FEDataGenerator)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
    ADD_PARAMETER(m_density0, FE_RANGE_GREATER(0.0), "density");
    ADD_PARAMETER(m_dt, FE_RANGE_GREATER(0.0), "time_step");
END_FECORE_CLASS();

FEScaledDensityMapGenerator::FEScaledDensityMapGenerator(FEModel* fem) : FEDataGenerator(fem)
{	
}

FEScaledDensityMapGenerator::~FEScaledDensityMapGenerator()
{
}

bool FEScaledDensityMapGenerator::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	return FEDataGenerator::Init();
}

// generate the data array for the given element set
bool FEScaledDensityMapGenerator::Generate(FEDomainMap& map)
{
	const FEElementSet& set = *map.GetElementSet();

	FEMesh& mesh = *set.GetMesh();

	FEDataType dataType = map.DataType();
	if (dataType != FE_DOUBLE) return false;

	int storageFormat = map.StorageFormat();
	if (storageFormat != FMT_MATPOINTS) return false;

	int N = set.Elements();
	for (int iel = 0; iel < N; ++iel)
	{
		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID(set[iel]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

        // get minimum element edge lenght (squared)
        double h = get_elem_hmin_alt(mesh, el);

        // loop Gauss points
        int nint = el.GaussPoints();
        for (int j = 0; j < nint; ++j)
        {

            // get material parameters
            FEMaterialPoint* pt = el.GetMaterialPoint(j);
            double mE = m_E(*pt);
            double mv = m_v(*pt);
            double density0 = m_density0(*pt);

            // lame parameters
            double lam = (mv*mE/((1+mv)*(1-2*mv)));
            double mu  = (0.5*mE/(1+mv));

            // critical density at target dt
            double density = 0.25 * m_dt*m_dt * (lam+2*mu) / (h*h); // 0.25 == 1/(2**2) -> scheme facotr
            if (density > density0){
                map.setValue(iel, j, density);  // density needs to be increased to achieve dt
            } else {
                map.setValue(iel, j, density0);  // density is already sufficiently high
            }

        } // for gauss point

	}
	return true;
}