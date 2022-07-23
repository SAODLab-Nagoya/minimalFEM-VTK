#include "element_constraint.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

void applyConstraint(SparseMatrix<double> &K, const vector<Node> &node,
                     VectorXd &load);
void output(char *outputPass, VectorXd &displacements,
            vector<double> &sigma_mises);

const int ndim = 2;  ///< dimension
const int voigt = 3; ///< number of rows of B&D matrix in Voigt expression

int numnp; ///< number of node
vector<Node> node;

// variables about element
int nelx; ///< number of element
vector<Element> element;

VectorXd load;

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