#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FEIsotropicElastic.h"
#include "FEBioMech/FEMechModel.h"

#include "FECore/FEAnalysis.h"
#include "FECore/FEMesh.h"
#include "FECore/FEModel.h"
#include "FECore/FESurfaceElementShape.h"
#include "FECore/log.h"

double get_elem_volume(FEMesh& mesh, FESolidElement& el){

    vec3d rt[FEElement::MAX_NODES];

    int neln = el.Nodes();
    for (int i = 0; i < neln; ++i) rt[i] = mesh.Node(el.m_node[i]).m_rt;

    int nint = el.GaussPoints();
    double* w = el.GaussWeights();
    double V = 0;
    for (int n = 0; n < nint; ++n)
    {
        // shape function derivatives
        double* Grn = el.Gr(n);
        double* Gsn = el.Gs(n);
        double* Gtn = el.Gt(n);

        // jacobian matrix
        double J[3][3] = { 0 };
        for (int i = 0; i < neln; ++i)
        {
            const double& Gri = Grn[i];
            const double& Gsi = Gsn[i];
            const double& Gti = Gtn[i];

            const double& x = rt[i].x;
            const double& y = rt[i].y;
            const double& z = rt[i].z;

            J[0][0] += Gri * x; J[0][1] += Gsi * x; J[0][2] += Gti * x;
            J[1][0] += Gri * y; J[1][1] += Gsi * y; J[1][2] += Gti * y;
            J[2][0] += Gri * z; J[2][1] += Gsi * z; J[2][2] += Gti * z;
        }

        // calculate the determinant
        double detJ0 = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
            + J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
            + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

        V += detJ0 * w[n];
    }

    return V;
    
}

//-----------------------------------------------------------------------------
double get_elem_hmin(FEMesh& mesh, FESolidElement& el){
        // get minimum edge lenght (squared)
        // TODO: wont work for quadratic elements
        double h2 = 1e99;
        int ne = el.Nodes();
        for (int i = 0; i < ne; ++i)
        {
            for (int j = 0; j < ne; ++j)
            {
                if (i != j)
                {
                    vec3d a = mesh.Node(el.m_node[i]).m_rt;
                    vec3d b = mesh.Node(el.m_node[j]).m_rt;

                    double L2 = (a - b).norm2();
                    if (L2 < h2) h2 = L2;
                }
            }
        }
        return sqrt(h2);
}

//-----------------------------------------------------------------------------
double get_elem_hmin_alt(FEMesh& mesh, FESolidElement& el){
    
    // get element volume
    double V = get_elem_volume(mesh, el);

    // get minimum element edge lenght
    double hmin = 1e99;

    // shape function derivative in center
    double Gr[3], Gs[3];
    int neln;
    double r, s, weight;
    switch (el.Shape())
    {
    case ET_TET4:
        {
            neln   = 3;
            r      = 1.0/3.0;
            s      = 1.0/3.0;
            weight = 0.5;
            FETri3 shape = FETri3();
            shape.shape_deriv(Gr, Gs, r, s);
        }
        break;
    case ET_TET10:
        {
            neln   = 6;
            r      = 1.0/3.0;
            s      = 1.0/3.0;
            weight = 0.5;
            FETri6 shape = FETri6();
            shape.shape_deriv(Gr, Gs, r, s);
        }
        break;
    case ET_HEX8:
        {
            neln   = 4;
            r      = 0.0;
            s      = 0.0;
            weight = 4.0;
            FEQuad4 shape = FEQuad4();
            shape.shape_deriv(Gr, Gs, r, s);
        }
        break;
    case ET_HEX20:
        {
            neln   = 8;
            r      = 0.0;
            s      = 0.0;
            weight = 4.0;
            FEQuad8 shape = FEQuad8();
            shape.shape_deriv(Gr, Gs, r, s);
        }
        break;
    default:
        printf("Error determining area of face of element shape type: %d (nodes=%d)\n", el.Shape(), el.Nodes());
        assert(false);
        break;
    }

    // loop faces
    int nfaces = el.Faces();
    int nf[FEElement::MAX_NODES];
    vec3d dxr, dxs;
    for (int j=0; j<nfaces; ++j)
    {																												
    	// get nodes of this face
    	el.GetFace(j, nf);  // nodelist
                                
    	// calculate area
    	dxr = dxs = vec3d(0,0,0);
    	for (int k=0; k<neln; ++k)
    	{
    		vec3d r0 = mesh.Node(nf[k]).m_r0;
    		dxr.x += Gr[k]*r0.x;
    		dxr.y += Gr[k]*r0.y;
    		dxr.z += Gr[k]*r0.z;

    		dxs.x += Gs[k]*r0.x;
    		dxs.y += Gs[k]*r0.y;
    		dxs.z += Gs[k]*r0.z;
    	}
    	double detJ = (dxr ^ dxs).norm();
    	double A = weight*detJ;
        double hrel = V / A;
    	if (hrel < hmin) hmin = hrel;
    }  // for face

    return hmin;
}
