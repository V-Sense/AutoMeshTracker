/*
  Modified by Matt Moynihan, https://github.com/mjkmoynihan
  Original Author: Yizhi Tang
  From: https://github.com/Tonsty/Non-Rigid-Registar
  Accessed September 7th 2018
*/

#ifndef NONRIGID_ICP_H
#define NONRIGID_ICP_H

#include <vtkCellLocator.h>

typedef vtkSmartPointer<vtkPolyData> VTKMesh;

typedef VTKMesh Template;
typedef VTKMesh Target;
typedef std::pair<int,int> Edge;
typedef boost::shared_ptr< std::vector<Edge> > Edges;
typedef vtkSmartPointer<vtkPoints> Vertices;
typedef boost::shared_ptr< std::vector<float> > Weights;

class NonRigidICP
{
public:
	NonRigidICP(
		Template _template, 
		Target _target):
		template_(_template),
		target_(_target){}
	
	~NonRigidICP(){}

	void init();

	void initCompute();

	void initCompute(const std::vector<std::array<double,3>>&);

	void edgesInit();

	void nearestSearchInit();

	void verticesInit();

	void correspondencesInit();

	void setCorrespondences(
		const std::vector<std::array<double,3>>& mt_matches);

	void weightsInit();

	int compute(float alpha, float gamma);

	inline Template GetTemplate(){
		return this->template_;
	}

protected:
	Template template_;
	Target target_;

	Edges edges_;
	Vertices vertices_;

	// corr[src_id] = target_id ? 
	Vertices correspondences_;
	Weights weights_;

private:
	vtkSmartPointer<vtkCellLocator> cellLocator_;
};

#endif
