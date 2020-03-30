#include <Engine/MeshEdit/AXAP.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

AXAP::AXAP(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMesh);
}

void AXAP::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool AXAP::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID)
	{
		printf("ERRIR::Minsurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	nV = triMesh->GetPositions().size();
	nT = triMesh->GetTriangles().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(nT);
	for (auto triangle : triMesh->GetTriangles())
	{
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	}
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary())
	{
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundary\n");
		heMesh->Clear();
		return false;
	}

	// positions of triangle mesh -> positions of half-edge structure
	for (size_t i = 0; i < nV; i++)
	{
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool AXAP::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitPara();
	InitFlatTri();
	InitParaSolver();
	iter_count_ = 0;

	UpdateTriMesh();
	return true;
}

void AXAP::UpdateTriMesh()
{
	// half-edge structure -> triangle mesh
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = vector<normalf>();
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nT);
	normals.reserve(nV);
	texcoords.reserve(nV);
	for (size_t i = 0; i < nV; i++)
	{
		//positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
		positions.push_back({ para_solution_(i, 0), para_solution_(i, 1), 0 });
		texcoords.push_back({ para_solution_(i, 0), para_solution_(i, 1) });
	}
	for (auto f : heMesh->Polygons())
	{
		for (auto v : f->BoundaryVertice())
		{
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
		}
	}

	triMesh->Init(indice, positions, normals, texcoords);
}

void AXAP::InitPara()
{
	SparseMatrix<float> A(nV, nV);
	MatrixXf B(nV, 2);
	// construct matrix A, B

	// Initialize
	for (size_t i = 0; i < nV; i++)
	{
		A.insert(i, i) = 0;
		B.row(i) = RowVector2f::Zero();
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			A.insert(i, heMesh->Index(vj)) = 0;
		}
	}

	// add cotangent weight
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
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
			A.coeffRef(i, i) += cot_theta;
			A.coeffRef(i, j) -= cot_theta;
			A.coeffRef(j, i) -= cot_theta;
			A.coeffRef(j, j) += cot_theta;
		}
	}

	// set boundary to given value
	size_t nB = heMesh->Boundaries()[0].size();
	float length_count = 0;
	for (size_t k = 0; k < nB; k++)
	{
		auto edge = heMesh->Boundaries()[0][k];
		double length = (edge->End()->pos - edge->Origin()->pos).norm();
		length_count += length;
	}

	bool flag = false;
	float length_pos = 0;
	for (size_t k = 0; k < nB; k++)
	{
		auto edge = heMesh->Boundaries()[0][k];
		auto vi = edge->Origin();
		size_t i = heMesh->Index(vi);
		float ratio = length_pos / length_count;
		if (k == 0)
		{
			start = i;
		}
		else if ((ratio >= 0.5) && !flag)
		{
			end = i;
			flag = true;
		}

		for (auto vj : vi->AdjVertices())
		{
			A.coeffRef(i, heMesh->Index(vj)) = 0;
		}
		A.coeffRef(i, i) = 1;
		if (ratio < 0.25)
		{
			B(i, 0) = ratio * 4;
			B(i, 1) = 0;
		}
		else if (ratio >= 0.25 && ratio < 0.5)
		{
			B(i, 0) = 1;
			B(i, 1) = (ratio - 0.25) * 4;
		}
		else if (ratio >= 0.5 && ratio < 0.75)
		{
			B(i, 0) = 1 - (ratio - 0.5) * 4;
			B(i, 1) = 1;
		}
		else
		{
			B(i, 0) = 0;
			B(i, 1) = 1 - (ratio - 0.75) * 4;
		}
		double length = (edge->End()->pos - edge->Origin()->pos).norm();
		length_pos += length;
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<float>> solver;
	solver.compute(A);
	para_solution_ = solver.solve(B);

	cout << "Success::Parametrize::Parametrize:" << endl
		<< "\t" << "parametrization successfully constructed" << endl;
}

void AXAP::InitFlatTri()
{
	flat_tri_ = (Matrix<float, 3, 2>*)malloc(nT * sizeof(Matrix<float, 3, 2>));
	if (flat_tri_ == nullptr)
	{
		cout << "No enough memory!" << endl;
	}

	for (size_t i = 0; i < nT; i++)
	{
		auto ti = heMesh->Polygons().at(i);
		auto vi0 = ti->HalfEdge()->Origin();
		auto vi1 = ti->HalfEdge()->End();
		auto vi2 = ti->HalfEdge()->Next()->End();
		vecf3 ei1 = vi1->pos - vi0->pos;
		vecf3 ei2 = vi2->pos - vi0->pos;
		double cos_theta = ei1.cos_theta(ei2);
		flat_tri_[i] = Matrix<float, 3, 2>();
		flat_tri_[i].row(0) = RowVector2f::Zero();
		flat_tri_[i].row(1) = RowVector2f(ei1.norm(), 0);
		flat_tri_[i].row(2) = RowVector2f(ei2.norm() * cos_theta, ei2.norm() * sqrt(1 - pow(cos_theta, 2)));
	}
}

void AXAP::InitParaSolver()
{
	SparseMatrix<float> A(nV, nV);
	for (size_t i = 0; i < nV; i++)
	{
		if (i == start || i == end)
		{
			A.insert(i, i) = 1;
		}
		else
		{
			auto vi = heMesh->Vertices().at(i);
			A.insert(i, i) = 0;
			for (auto vj : vi->AdjVertices())
			{
				A.insert(i, heMesh->Index(vj)) = 0;
			}
		}
	}

	// add cotangent weight
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
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
			if (i != start && i != end)
			{
				A.coeffRef(i, i) += cot_theta;
				A.coeffRef(i, j) -= cot_theta;
			}
			if (j != start && j != end)
			{
				A.coeffRef(j, i) -= cot_theta;
				A.coeffRef(j, j) += cot_theta;
			}
		}
	}
	para_solver_.compute(A);
}

Matrix2f AXAP::Jacobian(Matrix<float, 3, 2>& x, Matrix<float, 3, 2>& u)
{
	Matrix2f A = Matrix2f();
	Matrix2f B = Matrix2f();
	A.row(0) = x.row(1) - x.row(0);
	A.row(1) = x.row(2) - x.row(0);
	B.row(0) = u.row(1) - u.row(0);
	B.row(1) = u.row(2) - u.row(0);
	Matrix2f J =  A.lu().solve(B);
	return J;
}

bool AXAP::Iterate()
{
	return true;
}

void AXAP::UpdatePara()
{

}