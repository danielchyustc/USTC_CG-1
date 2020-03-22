#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <Engine/MeshEdit/AXAP.h>
#include <UGM/UGM>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class MinSurf;
	class Paramaterize;
	class AXAP;

	// mesh boundary == 1
	class ASAP : public AXAP {
	public:
		ASAP(Ptr<TriMesh> triMesh) : AXAP(triMesh) {}
	public:
		static const Ptr<ASAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<ASAP>(triMesh);
		}
	public:
		bool Iterate();

	protected:
		void UpdatePara();

	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P : public TPolygon<V, E, P> { };\
	};
}

