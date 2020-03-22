#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include <UI/Attribute.h>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class MinSurf;
	class AXAP;

	// mesh boundary == 1
	class Paramaterize : public HeapObj {
	public:
		Paramaterize(Ptr<TriMesh> triMesh);
		Paramaterize(Ptr<TriMesh> triMesh, WeightType weight_type, BoundShape bound_shape);
	public:
		static const Ptr<Paramaterize> New(Ptr<TriMesh> triMesh, WeightType weight_type, BoundShape bound_shape) {
			return Ubpa::New<Paramaterize>(triMesh, weight_type, bound_shape);
		}
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);

		bool Run();

	private:
		void Parametrize();

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
	private:
		float CotangentWeight(V* vi, V* vj);
	private:
		friend class Minsurf;
		friend class AXAP;

		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh;
		WeightType weight_type_;
		BoundShape bound_shape_;
		MatrixXf tex_solution_;
	};
}
