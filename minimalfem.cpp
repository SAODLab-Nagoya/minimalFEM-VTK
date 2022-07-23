#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

const int ndim = 2;  ///< dimension
const int voigt = 3; ///< number of rows of B&D matrix in Voigt expression

// Create a Bmatrix of triangular elements.
// see 竹内則雄ら「計算力学」section 6.1.3
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

// constraint class
struct Dirichlet
{
  int flag[3];   ///< 1=constrainted, 0=not constrainted
  double val[3]; ///< constrained value
};

/// node class
struct Node
{
  double x[3];                       ///< coordinate
  std::shared_ptr<Dirichlet> dirich; ///< Dirichlet constraint
};

// element class
class Element
{
public:
  // make element stiffness matrix's triplets
  void stiffnessMatrix(int ndim, int voigt, const Eigen::MatrixXd &D,
                       std::vector<Eigen::Triplet<double>> &triplets,
                       std::vector<Node> &nodes);

  // return element-wise mises stress
  double misesStress(int ndim, int voigt, const Eigen::MatrixXd &D,
                     std::vector<Node> &nodes, Eigen::VectorXd &displacement);

  const int ne = 3;    /// nodes in this element
  const int ipmax = 1; ///< number of gauss point

  int numdof;              ///< dof in this element
  std::vector<int> nodeID; ///< node ID in this element
};

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

int numnp; ///< number of node
vector<Node> node;

// variables about element
int nelx; ///< number of element
vector<Element> element;

VectorXd load;

void applyConstraint(SparseMatrix<double> &K, const vector<Node> &node,
                     VectorXd &load)
{
  vector<Triplet<double>> triplets;
  for(int i = 0; i < numnp; i++)
  {
    if(node[i].dirich == nullptr)
    {
      for(int j = 0; j < ndim; j++)
      {
        Triplet<double> tmp(ndim * i + j, ndim * i + j, 1);
        triplets.push_back(tmp);
      }
    }
    else
    {
      for(int j = 0; j < ndim; j++)
      {
        if(node[i].dirich->flag[j] == 0)
        {
          Triplet<double> tmp(ndim * i + j, ndim * i + j, 1);
          triplets.push_back(tmp);
        }
        // fix load to make it consistent with the constraint conditions.
        else
          load[ndim * i + j] = 0.0;
      }
    }
  }

  SparseMatrix<double> N(ndim * numnp, ndim * numnp);
  N.setFromTriplets(triplets.begin(), triplets.end());

  SparseMatrix<double> I(ndim * numnp, ndim * numnp);
  I.setIdentity();

  K = N.transpose() * K * N + (I - N);
}

void output(char *outputPass, VectorXd &displacements,
            vector<double> &sigma_mises)
{
  ofstream outfile(outputPass);
  // header
  outfile << "# vtk DataFile Version 4.0\n"
          << "output of FEM program\n"
          << "ASCII\n\n"
          << "DATASET UNSTRUCTURED_GRID\n";

  // coordinate
  outfile << "POINTS"
          << " " << numnp << " "
          << "float" << endl;
  for(int i = 0; i < numnp; i++)
  {
    outfile << node[i].x[0] << " " << node[i].x[1] << " " << 0.0 << endl;
  }
  outfile << endl;

  // connectivity
  outfile << "CELLS"
          << " " << nelx << " " << (element[0].ne + 1) * nelx << std ::endl;
  for(int i = 0; i < nelx; i++)
  {
    outfile << element[i].ne << " ";
    for(int j = 0; j < element[i].ne; j++)
    {
      outfile << element[i].nodeID[j] << " ";
    }
    outfile << endl;
  }
  outfile << endl;

  // cell shape(triangle is 5,square is 9)
  outfile << "CELL_TYPES"
          << " " << nelx << endl;
  for(int i = 0; i < nelx; i++)
  {
    outfile << 5 << endl;
  }
  outfile << endl;

  // displacement
  outfile << "POINT_DATA"
          << " " << numnp << endl;
  outfile << "VECTORS displacement float" << endl;
  if(ndim == 2)
    for(int i = 0; i < numnp; i++)
    {
      outfile << displacements[2 * i] << " " << displacements[2 * i + 1] << " "
              << 0.0 << endl;
    }
  else // 3D
    for(int i = 0; i < numnp; i++)
    {
      outfile << displacements[3 * i] << " " << displacements[3 * i + 1] << " "
              << displacements[3 * i + 2] << endl;
    }
  outfile << endl;

  // mises stress
  outfile << "CELL_DATA"
          << " " << nelx << endl
          << "SCALARS mises_stress float" << endl
          << "LOOKUP_TABLE default" << endl;
  for(int i = 0; i < nelx; i++)
  {
    outfile << sigma_mises[i] << endl;
  }

  outfile << endl;
}

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    cout << "usage: " << argv[0] << " <input file> <output file>\n";
    return 1;
  }

  ifstream infile(argv[1]);

  puts("input material");
  double poisson, young;
  infile >> poisson >> young;
  // set D matrix(plane stress)
  MatrixXd De(voigt, voigt);
  De << 1.0, poisson, 0.0, poisson, //
      1.0, 0.0, 0.0,                //
      0.0, (1.0 - poisson) / 2.0;
  De *= young / (1.0 - pow(poisson, 2.0));

  puts("input corrdinate");
  infile >> numnp;
  node.resize(numnp);
  for(int i = 0; i < numnp; ++i)
    infile >> node[i].x[0] >> node[i].x[1];

  puts("input connectivity");
  infile >> nelx;
  for(int i = 0; i < nelx; ++i)
  {
    Element actele;
    actele.numdof = ndim * actele.ne;
    actele.nodeID.resize(actele.ne);
    infile >> actele.nodeID[0] >> actele.nodeID[1] >> actele.nodeID[2];
    element.push_back(actele);
  }

  puts("input constraint");
  int nconst;
  infile >> nconst;
  for(int i = 0; i < nconst; ++i)
  {
    int nodeID, type;
    infile >> nodeID >> type;
    if(node[nodeID].dirich == nullptr)
    {
      shared_ptr<Dirichlet> dirich(new Dirichlet);
      // x
      if(type == 1)
      {
        dirich->flag[0] = 1;
        dirich->flag[1] = 0;
        dirich->flag[2] = 0;
      }
      // y
      else if(type == 2)
      {
        dirich->flag[0] = 0;
        dirich->flag[1] = 1;
        dirich->flag[2] = 0;
      }
      // x & y
      else if(type == 3)
      {
        dirich->flag[0] = 1;
        dirich->flag[1] = 1;
        dirich->flag[2] = 0;
      }
      else
      {
        cerr << "error in reading Dirichlet condition." << endl;
        exit(1);
      }

      node[nodeID].dirich = dirich;
    }
  }

  puts("input load");
  load.resize(ndim * numnp);
  load.setZero();
  int loadCount;
  infile >> loadCount;
  for(int i = 0; i < loadCount; ++i)
  {
    int node;
    infile >> node;
    for(int j = 0; j < ndim; j++)
    {
      double val;
      infile >> val;
      load[ndim * node + j] = val;
    }
  }

  puts("make stiffness matrix");
  vector<Triplet<double>> triplets;
  for(int i = 0; i < nelx; i++)
    element[i].stiffnessMatrix(ndim, voigt, De, triplets, node);

  puts("assembling");
  SparseMatrix<double> globalK(ndim * numnp, ndim * numnp);
  globalK.setFromTriplets(triplets.begin(), triplets.end());

  puts("fix matrix");
  applyConstraint(globalK, node, load);
  puts("solve Ku=f");
  SimplicialLDLT<SparseMatrix<double>> solver;
  solver.compute(globalK);
  VectorXd displacements = solver.solve(load);

  puts("make von-Mises stress");
  vector<double> mises(nelx);
  for(int i = 0; i < nelx; i++)
    mises[i] = element[i].misesStress(ndim, voigt, De, node, displacements);

  puts("make output");
  output(argv[2], displacements, mises);

  puts("finish.");
  return 0;
}
