#include <Engine/MeshEdit/ASAP.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;


bool ASAP::Iterate()
{
	iter_count_++;
	cout << iter_count_ << "th iteration" << endl;
	UpdatePara();
	UpdateTriMesh();
	return true;
}


void ASAP::UpdatePara()
{
	MatrixXf B(nV, 2);
	for (size_t i = 0; i < nV; i++)
	{
		B.row(i) = RowVector2f::Zero();
	}

	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		Matrix<float, 3, 2> xt = flat_tri_[t];
		Matrix<float, 3, 2> ut = Matrix<float, 3, 2>();
		ut.row(0) = para_solution_.row(heMesh->Index(edge->Origin()));
		ut.row(1) = para_solution_.row(heMesh->Index(edge->End()));
		ut.row(2) = para_solution_.row(heMesh->Index(edge->Next()->End()));
		Matrix2f J = Jacobian(xt, ut);
		JacobiSVD<Matrix2f> svd(J, ComputeThinU | ComputeThinV);
		Matrix2f U = svd.matrixU();
		Matrix2f V = svd.matrixV();
		Vector2f A = svd.singularValues();
		Matrix2f A2 = Matrix2f::Zero();
		A2(0, 0) = (A(0) + A(1)) / 2;
		A2(1, 1) = (A(0) + A(1)) / 2;
		Matrix2f L = U * V;

		for (int k = 0; k < 3; k++, edge = edge->Next())
		{
			auto vi = edge->Origin();
			auto vj = edge->End();
			auto vk = edge->Next()->End();
			int i = heMesh->Index(vi);
			int j = heMesh->Index(vj);
			vecf3 eki = vi->pos - vk->pos;
			vecf3 ekj = vj->pos - vk->pos;
			float cot_theta = eki.cos_theta(ekj) / eki.sin_theta(ekj);
			//B.row(i) += (para_solution_.row(i) - para_solution_.row(j)) * V.transpose() * V * cot_theta;
			B.row(i) += (xt.row(k) - xt.row((k + 1) % 3)) * L * cot_theta;
			//B.row(j) += (para_solution_.row(j) - para_solution_.row(i)) * V.transpose() * V * cot_theta;
			B.row(j) += (xt.row((k + 1) % 3) - xt.row(k)) * L * cot_theta;
		}
	}

	B.row(start) = RowVector2f::Zero();
	B.row(end) = RowVector2f(1, 1);
	para_solution_ = para_solver_.solve(B);
}

