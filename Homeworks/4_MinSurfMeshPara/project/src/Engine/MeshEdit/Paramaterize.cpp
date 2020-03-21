#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>()), weight_type_(kEqual), bound_shape_(kCircle)
{
	// TODO
	Init(triMesh);
}

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh, WeightType weight_type, BoundShape bound_shape)
	: heMesh(make_shared<HEMesh<V>>()), weight_type_(weight_type), bound_shape_(bound_shape)
{
	// TODO
	Init(triMesh);
}

void Paramaterize::Clear() {
	// TODO
	heMesh->Clear();
	triMesh = nullptr;
}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
	// TODO
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
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
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

bool Paramaterize::Run() {
	// TODO
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Parametrize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = vector<normalf>();
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	normals.reserve(nV);
	texcoords.reserve(nV);
	for (size_t i = 0; i < nV; i++)
	{
		positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
		texcoords.push_back({ tex_solution_(i, 0), tex_solution_(i, 1) });
	}
	for (auto f : heMesh->Polygons())
	{
		for (auto v : f->BoundaryVertice())
		{
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
		}
	}

	triMesh->Init(indice, positions, normals, texcoords);
	return true;
}

void Paramaterize::Parametrize()
{
	size_t nV = heMesh->Vertices().size();
	SparseMatrix<float> A(nV, nV);
	MatrixXf B(nV, 2);
	// construct matrix A, B
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		A.insert(i, i) = 0;
		B.row(i) = Vector3f::Zero();
		for (auto vj : vi->AdjVertices())
		{
			float weight;
			switch (weight_type_)
			{
			case kEqual:
				weight = 1;
				break;
			case kCotangent:
				weight = CotangentWeight(vi, vj);
				break;
			default:
				break;
			}
			A.coeffRef(i, i) += weight;
			A.insert(i, heMesh->Index(vj)) = -weight;
		}
	}

	// set boundary to given value
	size_t nB = heMesh->Boundaries()[0].size();
	for (size_t k = 0; k < nB; k++)
	{
		auto vi = heMesh->Boundaries()[0][k]->Origin();
		size_t i = heMesh->Index(vi);
		for (auto vj : vi->AdjVertices())
		{
			A.coeffRef(i, heMesh->Index(vj)) = 0;
		}
		A.coeffRef(i, i) = 1;
		float ratio;
		switch (bound_shape_)
		{
		case kCircle:
			B(i, 0) = 0.5 + 0.5 * cos(k * 2 * 3.1415926 / nB);
			B(i, 1) = 0.5 + 0.5 * sin(k * 2 * 3.1415926 / nB);
			break;
		case kSquare:
			ratio = (float)k / nB;
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
			break;
		default:
			break;
		}
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<float>> solver;
	solver.compute(A);
	tex_solution_ = solver.solve(B);

	cout << "Success::Parametrize::Parametrize:" << endl
		<< "\t" << "parametrization successfully constructed" << endl;
}

float Paramaterize::CotangentWeight(V* vi, V* vj)
{
	V* v_pre = vi->HalfEdgeTo(vj)->RotatePre()->End();
	V* v_next = vi->HalfEdgeTo(vj)->RotateNext()->End();
	vecf3 v_pi = vi->pos - v_pre->pos;
	vecf3 v_pj = vj->pos - v_pre->pos;
	vecf3 v_ni = vi->pos - v_next->pos;
	vecf3 v_nj = vj->pos - v_next->pos;
	float cot_alpha = v_pi.cos_theta(v_pj) / v_pi.sin_theta(v_pj);
	float cot_beta = v_ni.cos_theta(v_nj) / v_ni.sin_theta(v_nj);
	float weight = cot_alpha + cot_beta;
	return weight > 0 ? weight : 0;
}
