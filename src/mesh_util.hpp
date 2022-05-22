#pragma once

#include "Common/hexahedral_mesh.hpp"
#include "Common/tetrahedral_mesh.hpp"

#include <utility>
#include <map>
#include <fstream>

namespace materials {
	template<typename T>
	HexahedralMesh<T> vox2hex(const std::string & filename) {
		// read the voxel file
		std::ifstream fin(filename);
		Eigen::Matrix<T,3,1> _pmin;
		T _dx;
		int _nx, _ny, _nz;
		std::vector<std::vector<std::string>> _corner_flags;
		fin >> _pmin.x() >> _pmin.y() >> _pmin.z();
		fin >> _dx;
		fin >> _nx >> _ny >> _nz;
		_corner_flags.resize(_nx + 1);
		for (int i = 0; i < _nx + 1; ++i) {
			_corner_flags[i].resize(_ny + 1);
			for (int j = 0; j < _ny + 1; ++j) {
				fin >> _corner_flags[i][j];
			}
		}

		// Convert voxels into vertices and elements
		Matrix3X<T> V;
		Matrix8Xi<T> E;
		V.resize(Eigen::NoChange, _nx*_ny);
		E.resize(Eigen::NoChange, _nx*_ny);
		std::map<std::tuple<int, int, int>, int> index;
		int nv = 0; // Number of vertices
		auto v_index = [&](std::tuple<int, int, int> v)->int {
			if (index[v] == 0) {
				int i, j, k;
				std::tie(i, j, k) = v;
				Eigen::Matrix<T,3,1> vec;
				vec(0,0) = (T)i; vec(1,0) = T(j); vec(2,0) = T(k);
				if (V.cols() <= nv) {
					V.conservativeResize(Eigen::NoChange, V.cols() * 2);
				}
				V.col(nv) = _pmin + _dx * vec;
				index[v] = nv + 1;
				++nv;
			}
			return index[v] - 1;
		};
		std::cout << "Converting Voxels to Hex Mesh" << std::endl;
		int c = 0; // Number of columns
		for (int i = 0; i < _nx; ++i) {
			for (int j = 0; j < _ny; ++j) {
				for (int k = 0; k < _nz; ++k) {
					if (_corner_flags[i][j][k] == '1') {
						if (E.cols() <= c) {
							E.conservativeResize(Eigen::NoChange, 2 * E.cols());
						}
						E(0, c) = v_index(std::make_tuple(i, j, k));
						E(1, c) = v_index(std::make_tuple(i, j, k + 1));
						E(2, c) = v_index(std::make_tuple(i, j + 1, k));
						E(3, c) = v_index(std::make_tuple(i, j + 1, k + 1));
						E(4, c) = v_index(std::make_tuple(i + 1, j, k));
						E(5, c) = v_index(std::make_tuple(i + 1, j, k + 1));
						E(6, c) = v_index(std::make_tuple(i + 1, j + 1, k));
						E(7, c) = v_index(std::make_tuple(i + 1, j + 1, k + 1));
						++c;
					}
				}
			}
		}
		std::cout << std::endl;
		V.conservativeResize(Eigen::NoChange, nv);
		E.conservativeResize(Eigen::NoChange, c);

		return HexahedralMesh<T>(V, E);
	}

	template<typename T>
	TetrahedralMesh<T> hex2tet(const HexahedralMesh<T> & hex_mesh) {
		Matrix3X<T> V = hex_mesh.vertex();
		const auto E = hex_mesh.element();
		Eigen::Matrix4Xi tet_el;
		tet_el.resize(4, 5 * E.cols());
		for (int e = 0; e < E.cols(); ++e) {
			tet_el.col(5 * e + 0) = Eigen::Vector4i(E(0, e), E(4, e), E(2, e), E(1, e));
			tet_el.col(5 * e + 1) = Eigen::Vector4i(E(2, e), E(1, e), E(7, e), E(3, e));
			tet_el.col(5 * e + 2) = Eigen::Vector4i(E(1, e), E(4, e), E(7, e), E(5, e));
			tet_el.col(5 * e + 3) = Eigen::Vector4i(E(2, e), E(7, e), E(4, e), E(6, e));
			tet_el.col(5 * e + 4) = Eigen::Vector4i(E(2, e), E(4, e), E(7, e), E(1, e));
		}
		std::cout << std::endl;
		return TetrahedralMesh<T>(V, tet_el);
	}

	template <typename T>
	void tet2tri(
		const TetrahedralMesh<T> & tet_mesh,
		Eigen::MatrixXd& V,
		Eigen::MatrixXi& F)
	{
		std::cout << "Converting tet mesh to tri mesh" << std::endl;
		V = tet_mesh.vertex().transpose().eval();
		const auto E = tet_mesh.element();
		F.resize(4 * E.cols(), 3);
		for (int c = 0; c < E.cols(); ++c) {
			F.row(4 * c + 0) = Eigen::Vector3i(E(2, c), E(1, c), E(0, c));
			F.row(4 * c + 1) = Eigen::Vector3i(E(3, c), E(2, c), E(0, c));
			F.row(4 * c + 2) = Eigen::Vector3i(E(1, c), E(3, c), E(0, c));
			F.row(4 * c + 3) = Eigen::Vector3i(E(2, c), E(3, c), E(1, c));
		}
		std::cout << std::endl;
	}
}