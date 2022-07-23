#include "element_constraint.hpp"

using namespace std;
using namespace Eigen;

void Element::stiffnessMatrix(int ndim, int voigt, const MatrixXd &D,
                              vector<Triplet<double>> &triplets,
                              vector<Node> &nodes)
{
  // coordinate
  MatrixXd X(ne, ndim);
  for(int i = 0; i < ne; i++)
    for(int j = 0; j < ndim; j++)
      X(i, j) = nodes[nodeID[i]].x[j];

  MatrixXd K(numdof, numdof);
  K.setZero();
  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    // get Bmatrix & jacobian & gauss weight
    Eigen::MatrixXd B(voigt, numdof);
    B.setZero();
    double jac = 0.0, weight = 0.0;
    bmatrixTri3(B, X, jac, weight, ip);

    K += B.transpose() * D * B * jac * weight;
  }

  // global dof's array
  VectorXi idof(numdof);
  for(int i = 0; i < ne; i++)
    for(int j = 0; j < ndim; j++)
      idof[ndim * i + j] = ndim * nodeID[i] + j;

  for(int i = 0; i < numdof; i++)
  {
    for(int j = 0; j < numdof; j++)
    {
      Triplet<double> trplt(idof[i], idof[j], K(i, j));
      triplets.push_back(trplt);
    }
  }
}

double Element::misesStress(int ndim, int voigt, const MatrixXd &D,
                            vector<Node> &nodes, VectorXd &displacement)
{
  // coordinate & displacement
  MatrixXd X(ne, ndim);
  VectorXd disp(ndim * ne);
  for(int i = 0; i < ne; i++)
    for(int j = 0; j < ndim; j++)
    {
      X(i, j) = nodes[nodeID[i]].x[j];
      disp[ndim * i + j] = displacement[ndim * nodeID[i] + j];
    }

  Eigen::VectorXd sigma(voigt);
  sigma.setZero();
  for(int ip = 0; ip < ipmax; ip++)
  {
    Eigen::MatrixXd B(voigt, numdof);
    B.setZero();
    double jac = 0.0, weight = 0.0;
    bmatrixTri3(B, X, jac, weight, ip);
    sigma += D * B * disp / (double)(ipmax);
  }

  // Mises stress (2D)
  assert(ndim == 2);
  return sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] +
              3.0 * sigma[2] * sigma[2]);
}

void bmatrixTri3(Eigen::MatrixXd &B, Eigen::MatrixXd &X, double &jac,
                 double &weight, int ip)
{
  // C is the area in this element
  MatrixXd C(3, 3);
  VectorXd ones(3);
  ones.setOnes();
  //    1  x1  y1
  // C= 1  x2  y2
  //    1  x3  y3
  C << ones, X;

  MatrixXd IC = C.inverse();

  for(int i = 0; i < 3; i++)
  {
    B(0, 2 * i + 0) = IC(1, i);
    B(0, 2 * i + 1) = 0.0;
    B(1, 2 * i + 0) = 0.0;
    B(1, 2 * i + 1) = IC(2, i);
    B(2, 2 * i + 0) = IC(2, i);
    B(2, 2 * i + 1) = IC(1, i);
  }
  jac = C.determinant();
  weight = 0.5;
}