#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class MinSurf;
	class Paramaterize;

	// mesh boundary == 1
	class AXAP : public HeapObj {
	public:
		AXAP(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<AXAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<AXAP>(triMesh);
		}
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();
		virtual bool Iterate();

	protected:
		void InitPara();
		void UpdateTriMesh();
		void InitFlatTri();
		void InitParaSolver();
		virtual void UpdatePara();
		Matrix2f Jacobian(Matrix<float, 3, 2>& x, Matrix<float, 3, 2>& u);

	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P : public TPolygon<V, E, P> { };
	protected:
		size_t							nV;
		size_t							nT;
		int								start;
		int								end;
		int								iter_count_;
		Ptr<TriMesh>					triMesh;
		const Ptr<HEMesh<V>>			heMesh;
		Matrix<float, 3, 2>				*flat_tri_;
		SparseLU<SparseMatrix<float>>	para_solver_;
		MatrixXf						para_solution_;
	};
}

