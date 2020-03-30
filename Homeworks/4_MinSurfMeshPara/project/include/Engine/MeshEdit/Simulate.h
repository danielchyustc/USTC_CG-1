#pragma once

#include <Basic/HeapObj.h>
//#include <Engine/Primitive/MassSpring.h>
#include <UGM/UGM>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

namespace Ubpa {
	class Simulate : public HeapObj {
	public:
		Simulate(const std::vector<pointf3>& plist,
			const std::vector<unsigned>& elist) {
			edgelist = elist;
			this->positions.resize(plist.size());
			for (int i = 0; i < plist.size(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					this->positions[i][j] = plist[i][j];
				}
			}
		};
	public:
		static const Ptr<Simulate> New(const std::vector<pointf3>& plist,
			const std::vector<unsigned> &elist) {
			return Ubpa::New<Simulate>(plist, elist);
		}
	public:
		// clear cache data
		void Clear();

		// init cache data (eg. half-edge structure) for Run()
		bool Init();
		//bool Init();

		// call it after Init()
		bool Run();
		
		const std::vector<pointf3>& GetPositions() const { return positions; };

		const float GetStiff() { return stiff; };
		void SetStiff(float k) { stiff = k; Init();};
		const float GetTimeStep() { return h; };
		void SetTimeStep(float k) { h = k; Init();};
		std::vector<unsigned>& GetFix() { return this->fixed_id; };
		void SetFix(const std::vector<unsigned>& f) { this->fixed_id = f; Init();};
		const std::vector<pointf3>& GetVelocity() { return velocity; };
		//void SetVelocity(const std::vector<pointf3>& v) { velocity = v; };

		void SetLeftFix();


	private:
		// kernel part of the algorithm
		void SimulateOnce();
		VectorXf FunctionG(VectorXf x);
		void UpdateGradient(VectorXf x);

	private:
		float					h = 0.03f;  //²½³¤
		float					stiff = 1e5;	//spring stiffness
		float					G = 9.8;	//acceleration of gravity
		std::vector<unsigned>	fixed_id = std::vector<unsigned>();  //fixed point id


		//mesh data
		int						nV;				//number of vertices
		int						nE;				//number of edges
		int						mV;				//number of free vertices
		std::vector<unsigned>	edgelist;
		VectorXf				original_length;


		//simulation data
		std::vector<pointf3>	positions;
		std::vector<pointf3>	velocity;
		VectorXf				y;
		SparseMatrix<float>		K;
		SparseMatrix<float>		gradient;
		
	};
}
