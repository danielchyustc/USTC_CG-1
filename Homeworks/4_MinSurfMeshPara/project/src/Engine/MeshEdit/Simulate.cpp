#include <Engine/MeshEdit/Simulate.h>
#include <iostream>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;


void Simulate::Clear() {
	this->positions.clear();
	this->velocity.clear();
}

bool Simulate::Init() {
	cout << "fixed id: " << fixed_id.size() << endl;
	nV = positions.size();
	mV = nV - fixed_id.size();
	nE = edgelist.size() / 2;
	cout << "nV: " << nV << ",	mV: " << mV << endl;
	this->velocity .resize(nV);
	for (int i = 0; i < nV; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->velocity[i][j] = 0;
		}
	}

	original_length.resize(nE);
	for (int k = 0; k < nE; k++)
	{
		int i = edgelist[2 * k];
		int j = edgelist[2 * k + 1];
		original_length[k] = (positions[i] - positions[j]).norm();
	}

	K.resize(3 * mV, 3 * nV);
	for (int i = 0, fix_index = 0, row_index = 0; i < nV; i++)
	{
		if (fix_index == fixed_id.size() || i != fixed_id[fix_index])
		{
			K.insert(3 * row_index, 3 * i) = 1;
			K.insert(3 * row_index + 1, 3 * i + 1) = 1;
			K.insert(3 * row_index + 2, 3 * i + 2) = 1;
			row_index++;
		}
		else
		{
			fix_index++;
		}
	}

	y.resize(3 * nV);
	gradient.resize(3 * nV, 3 * nV);

	return true;
}

bool Simulate::Run() {
	SimulateOnce();

	// half-edge structure -> triangle mesh

	return true;
}

void Ubpa::Simulate::SetLeftFix()
{
	//固定网格x坐标最小点
	fixed_id.clear();
	double x = 100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][0] < x)
		{
			x = positions[i][0];
		}
	}

	for (int i = 0; i < positions.size(); i++)
	{
		if (abs(positions[i][0] - x) < 1e-5)
		{
			fixed_id.push_back(i);
		}
	}

	Init();
}

void Simulate::SimulateOnce() {
	for (int i = 0, fix_index = 0; i < nV; i++)
	{
		if (fix_index == fixed_id.size() || i != fixed_id[fix_index])
		{
			y[3 * i] = positions[i][0] + h * velocity[i][0];
			y[3 * i + 1] = positions[i][1] + h * velocity[i][1] - h * h * G;
			y[3 * i + 2] = positions[i][2] + h * velocity[i][2];
		}
		else
		{
			y[3 * i] = positions[i][0];
			y[3 * i + 1] = positions[i][1];
			y[3 * i + 2] = positions[i][2];
			fix_index++;
		}
	}

	VectorXf x = y;
	VectorXf xf = K * x;
	VectorXf b = x - K.transpose() * xf;
	VectorXf g = FunctionG(x);
	for (int i = 0; g.norm() > mV * 1e-3 && i < 3; i++)
	{
		UpdateGradient(x);
		SparseLU<SparseMatrix<float>> solver;
		solver.compute(K * gradient * K.transpose());
		MatrixXf delta_x = solver.solve(g);
		xf = xf - delta_x;
		x = K.transpose() * xf + b;
		g = FunctionG(x);
	}
	
	for (int i = 0; i < nV; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			velocity[i][j] = (x[3 * i + j] - positions[i][j]) / h;
			positions[i][j] = x[3 * i + j];
		}
	}

}

VectorXf Simulate::FunctionG(VectorXf x)
{
	VectorXf g = x - y;
	for (int k = 0; k < edgelist.size() / 2; k++)
	{
		int i = edgelist[2 * k];
		int j = edgelist[2 * k + 1];
		Vector3f xi = Vector3f(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
		Vector3f xj = Vector3f(x[3 * j], x[3 * j + 1], x[3 * j + 2]);
		Vector3f r = xi - xj;
		Vector3f f = stiff * (r.norm() - original_length[k]) / r.norm() * r;
		for (int l = 0; l < 3; l++)
		{
			g[3 * i + l] += h * h * f[l];
			g[3 * j + l] -= h * h * f[l];
		}
	}
	return K * g;
}

void Simulate::UpdateGradient(VectorXf x)
{
	gradient.setIdentity();
	for (int k = 0; k < edgelist.size() / 2; k++)
	{
		int i = edgelist[2 * k];
		int j = edgelist[2 * k + 1];
		Vector3f xi = Vector3f(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
		Vector3f xj = Vector3f(x[3 * j], x[3 * j + 1], x[3 * j + 2]);
		Vector3f r = xi - xj;
		Matrix3f df_dx = stiff * (original_length[k] / r.norm() - 1) * Matrix3f::Identity()
			- stiff * original_length[k] * pow(r.norm(), -3) * r * r.transpose();
		for (int m = 0; m < 3; m++)
		{
			for (int n = 0; n < 3; n++)
			{
				gradient.coeffRef(3 * i + m, 3 * i + n) -= h * h * df_dx(m, n);
				gradient.coeffRef(3 * i + m, 3 * j + n) += h * h * df_dx(m, n);
				gradient.coeffRef(3 * j + m, 3 * i + n) += h * h * df_dx(m, n);
				gradient.coeffRef(3 * j + m, 3 * j + n) -= h * h * df_dx(m, n);
			}
		}
	}
}
