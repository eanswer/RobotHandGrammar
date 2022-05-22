// Tao Du
// taodu@csail.mit.edu
// Oct 5, 2016
#pragma once
#include "Common/poly_mesh.hpp"
#include "Common/typedefs.hpp"

namespace materials {





    template<typename T>
    class HexahedralMesh : public PolyMesh<3, T> {



    public:


        HexahedralMesh(const Matrix3X<T>& vertex, const Matrix8Xi<T>& element)
        : PolyMesh<3, T>(vertex, element,
                        (Eigen::Matrix<int, 2, 12>() << 0, 0, 1, 2, 4, 4, 5, 6, 0, 1, 2, 3,
                                                        1, 2, 3, 3, 5, 6, 7, 7, 4, 5, 6, 7).finished()) {}

        explicit HexahedralMesh(const std::string& file_name) : PolyMesh<3, T>(file_name) {}
        HexahedralMesh(const HexahedralMesh& hex_mesh) : PolyMesh<3, T>(hex_mesh.vertex_, hex_mesh.element_,
                                                                     hex_mesh.edge_in_element_) {}

        HexahedralMesh& operator=(const HexahedralMesh& hex_mesh) {
            PolyMesh<3, T>::operator=(hex_mesh);
            return *this;
        }
        ~HexahedralMesh() {}
    };

    // Do not use HexahedralMeshCuboid below to build a rigid or soft body, as
// HexRigidBody and HexDeformableBody assumes all hex elements are cubes.
// TODO: fix this inconsistency.
    template<typename T>
    const HexahedralMesh<T> HexahedralMeshCuboid(const Vector3<T>& pmin,
                                                 const Eigen::Vector3i& vertex_counts, const Vector3<T>& size){
        const int nx = vertex_counts.x();
        const int ny = vertex_counts.y();
        const int nz = vertex_counts.z();
        Matrix3X<T> vertex = MatrixX<T>::Zero(3, nx * ny * nz);
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k)
                    vertex.col(i * ny * nz + j * nz + k) = pmin
                                                           + size.cwiseProduct(Vector3<T>(i, j, k));
        Matrix8Xi<T> element = Eigen::MatrixXi::Zero(8,
                                                     (nx - 1) * (ny - 1) * (nz - 1));
        for (int i = 0; i < nx - 1; ++i)
            for (int j = 0; j < ny - 1; ++j)
                for (int k = 0; k < nz - 1; ++k) {
                    const int col_index = i * (ny - 1) * (nz - 1) + j * (nz - 1) + k;
                    element(0, col_index) = i * ny * nz + j * nz + k;
                    element(1, col_index) = i * ny * nz + j * nz + k + 1;
                    element(2, col_index) = i * ny * nz + (j + 1) * nz + k;
                    element(3, col_index) = i * ny * nz + (j + 1) * nz + k + 1;
                    element(4, col_index) = (i + 1) * ny * nz + j * nz + k;
                    element(5, col_index) = (i + 1) * ny * nz + j * nz + k + 1;
                    element(6, col_index) = (i + 1) * ny * nz + (j + 1) * nz + k;
                    element(7, col_index) = (i + 1) * ny * nz + (j + 1) * nz + k + 1;
                }
        return HexahedralMesh<T>(vertex, element);
    }

    template<typename T>
    const HexahedralMesh<T> HexahedralMeshCuboid(const Vector3<T>& pmin,
                                const Eigen::Vector3i& vertex_counts, const T size){
        return HexahedralMeshCuboid(pmin, vertex_counts, Vector3<T>(size, size, size));
    }


}

