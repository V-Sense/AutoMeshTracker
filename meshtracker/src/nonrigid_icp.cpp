/*
  Modified by Matt Moynihan, https://github.com/mjkmoynihan
  Original Author: Yizhi Tang
  From: https://github.com/Tonsty/Non-Rigid-Registar
  Accessed September 7th 2018
*/

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkCellLocator.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <vector>

#include <boost/shared_ptr.hpp>

#include "nonrigid_icp.h"
#include "log.h"

void NonRigidICP::init()
{
	edgesInit();
	verticesInit();
	nearestSearchInit();
}

void NonRigidICP::initCompute()
{
	correspondencesInit();
	weightsInit();
}

void NonRigidICP::initCompute
(
	const std::vector<std::array<double,3>>& mt_matches
)
{
	setCorrespondences(mt_matches);
	weightsInit();
}

void NonRigidICP::edgesInit()
{
	if (edges_ == NULL)
	{
		edges_.reset(new std::vector<Edge>);
	} 

	for (int i = 0; i < template_->GetNumberOfCells(); ++i)
	{
		vtkCell* cell = template_->GetCell(i);
		int a = cell->GetPointId(0);
		int b = cell->GetPointId(1);
		int c = cell->GetPointId(2);

		edges_->push_back(Edge(a,b));
		edges_->push_back(Edge(b,c));
		edges_->push_back(Edge(c,a));
	}
}

void NonRigidICP::nearestSearchInit()
{
	cellLocator_ = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator_->SetDataSet(target_);
	cellLocator_->BuildLocator();
}

void NonRigidICP::verticesInit()
{
	vertices_ = template_->GetPoints();
}

void NonRigidICP::correspondencesInit()
{
	 

	// TODO: Check this doesn't screw up the intermediate step processes 
	if (correspondences_ != NULL)
	{
		return; 
	}

	correspondences_ = Vertices::New();

	correspondences_->SetNumberOfPoints(vertices_->GetNumberOfPoints());
	for (int i = 0; i < vertices_->GetNumberOfPoints(); i++)
	{
		double testPoint[3];
		vertices_->GetPoint(i, testPoint);
		double closestPoint[3];
		double closestPointDist2;
		vtkIdType cellId;
		int subId;
		cellLocator_->FindClosestPoint(testPoint, closestPoint, cellId, subId, closestPointDist2);
		correspondences_->SetPoint(i, closestPoint);
	}
}

void NonRigidICP::setCorrespondences
(
	const std::vector<std::array<double,3>>& mt_matches
)
{

	if (mt_matches.size() != vertices_->GetNumberOfPoints())
	{
		std::string e = 
		"In NonRigidICP::setCorrespondences. Matches size different to mesh size.";
		BOOST_LOG_TRIVIAL(error) << e;
		throw std::runtime_error(e);
	}

	correspondences_ = Vertices::New();
	correspondences_->SetNumberOfPoints(vertices_->GetNumberOfPoints());
	for (size_t pt = 0; pt < mt_matches.size(); pt++)
	{
		std::array<double,3> matchPoint = mt_matches[pt];
		//TODO: replace this
		double point[3];
		point[0] = matchPoint[0];
		point[1] = matchPoint[1];
		point[2] = matchPoint[2];
		correspondences_->SetPoint(pt,point);
	}
}

// initialize all weights to 1
// Should modify this: wrt the paper, weights may refer to the correspondence reliability
// in which case a value of 1 assumes perfect 1-to-1 matches. 
// TODO: add a flag function to SetCorrespondences which can be read to set the 
// relevant weight indices. i.e. weight could be some function of the arg between 
// the matched point dist and normal
void NonRigidICP::weightsInit()
{
	if (weights_ == NULL)
	{
		weights_.reset(new std::vector<float>());
	}

	weights_->resize(vertices_->GetNumberOfPoints());

	for (int i = 0; i < weights_->size(); ++i)
	{
		(*weights_)[i] = 1.0f;
	}
}


// alpha: stiffness, gamma: weight diff between rot/skew agains trans
int NonRigidICP::compute(
	float _alpha, 
	float _gamma
)
{
	int n = vertices_->GetNumberOfPoints();
	int m = edges_->size();

	Eigen::SparseMatrix<float> A(4*m + n, 4*n);

	// Global Cost Function : E(X) = || [alpha M x G, WD] X - [0, WU] || ^2
	// 							   = || AX - B || ^ 2 (linearized)


	// Set up alpha M x G
	std::vector< Eigen::Triplet<float> > alpha_M_G;
	for (int i = 0; i < m; ++i){
		Edge edge = (*edges_)[i];
		int a = edge.first;
		int b = edge.second;

		for (int j = 0; j < 3; j++)
		{
 			alpha_M_G.push_back(Eigen::Triplet<float>(i*4 + j, a*4 + j, _alpha));
		}

		// apply gamma
		alpha_M_G.push_back(Eigen::Triplet<float>(i*4 + 3, a*4 + 3, _alpha * _gamma));

		for (int j = 0; j < 3; j++)
		{ 
			alpha_M_G.push_back(Eigen::Triplet<float>(i*4 + j, b*4 + j, -_alpha));
		}

		// apply gamma 
		alpha_M_G.push_back(Eigen::Triplet<float>(i*4 + 3, b*4 + 3, -_alpha * _gamma));
	}

	// set up W D 
	std::vector< Eigen::Triplet<float> > W_D;
	for (int i = 0; i < n; ++i)
	{
		double xyz[3];
		vertices_->GetPoint(i, xyz);

		float weight = (*weights_)[i];

		for (int j = 0; j < 3; ++j)
		{
			W_D.push_back(Eigen::Triplet<float>(4*m + i, i*4 + j, weight * xyz[j]));
		}

		W_D.push_back(Eigen::Triplet<float>(4*m + i, i*4 + 3, weight));
	}


	// We ignore the landmark term, beta D_L as we're not using tailored correspondences
	
	std::vector< Eigen::Triplet<float> > _A = alpha_M_G;
	_A.insert(_A.end(), W_D.begin(), W_D.end());

	A.setFromTriplets(_A.begin(), _A.end());

	// Set up last term of cost function, B = [0,WU] 
	Eigen::MatrixX3f B = Eigen::MatrixX3f::Zero(4*m + n, 3);
	for (int i = 0; i < n; ++i)
	{
		double xyz[3];
		correspondences_->GetPoint(i, xyz);

		float weight = (*weights_)[i];
		
		for (int j = 0; j < 3; j++)
		{
			B(4*m + i, j) = weight * xyz[j];
		}
	}

	// E(X) is minimised for X = (A^TA)^-1 A^TB

	Eigen::SparseMatrix<float> ATA = Eigen::SparseMatrix<float>(A.transpose()) * A;
	Eigen::MatrixX3f ATB = Eigen::SparseMatrix<float>(A.transpose()) * B;

	// Eigen::ConjugateGradient< Eigen::SparseMatrix<float> > solver;
	Eigen::ConjugateGradient< Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper > solver;

	solver.compute(ATA);

	if (solver.info()!=Eigen::Success)
	{
		BOOST_LOG_TRIVIAL(error) << "Decomposition failed";
		return 1;
	}

	Eigen::MatrixX3f X = solver.solve(ATB);

	Eigen::Matrix3Xf XT = X.transpose();
	
	for (int i = 0; i < n; ++i)
	{
		double xyz[3];
		vertices_->GetPoint(i, xyz);
		Eigen::Vector4f point(xyz[0], xyz[1], xyz[2], 1.0f);
		Eigen::Vector3f point_transformed = XT.block<3, 4>(0, 4*i) * point;
		vertices_->SetPoint(i, point_transformed[0], point_transformed[1], point_transformed[2]);

		double txyz[3];
		correspondences_->GetPoint(i, txyz);
	}

	return 0;
}
