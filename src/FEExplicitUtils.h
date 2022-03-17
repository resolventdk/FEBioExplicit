#pragma once 

#include "FECore/FEElement.h"
#include "FECore/FEMesh.h"

double get_elem_volume(FEMesh& mesh, FESolidElement& el);

double get_elem_hmin2(FEMesh& mesh, FESolidElement& el);

double get_elem_hmin2_alt(FEMesh& mesh, FESolidElement& el);
