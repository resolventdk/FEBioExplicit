#pragma once 

#include "FECore/FEElement.h"
#include "FECore/FEMesh.h"

double get_elem_volume(FEMesh& mesh, FESolidElement& el);

double get_elem_hmin(FEMesh& mesh, FESolidElement& el);

double get_elem_hmin_alt(FEMesh& mesh, FESolidElement& el);
