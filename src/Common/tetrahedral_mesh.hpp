// Jie Xu
// jiex@csail.mit.edu
// Aug. 21, 2019
#pragma once
#include "Common/poly_mesh.hpp"
#include "Common/typedefs.hpp"

namespace materials {
    template<typename T>
    class TetrahedralMesh : public PolyMesh<3, T> {
    public:
        TetrahedralMesh(const Matrix3X<T>& vertex, const Matrix4Xi<T>& element)
        : PolyMesh<3, T>(vertex, element,
                        (Eigen::Matrix<int, 2, 6>() << 0, 0, 0, 1, 1, 2,
                                                       1, 2, 3, 2, 3, 3).finished()) {}
        TetrahedralMesh(const TetrahedralMesh& tet_mesh) : PolyMesh<3, T>(tet_mesh.vertex_, tet_mesh.element_,
                                                                     tet_mesh.edge_in_element_) {}

        ~TetrahedralMesh() {};
        TetrahedralMesh& operator=(const TetrahedralMesh& tet_mesh) {
            PolyMesh<3, T>::operator=(tet_mesh);
            return *this;
        }
    };

    int get_cuboid_id(int x, int y, int z, const Eigen::Vector3i& vertex_counts) {
        return x * vertex_counts[1] * vertex_counts[2] + y * vertex_counts[2] + z;
    }

    // generate a cuboid tet mesh from hex parameters
    template<typename T>
    const TetrahedralMesh<T> TetrahedralMeshCuboid(const Eigen::Vector3i& vertex_counts, const T size) {
        const int vertex_num = vertex_counts[0] * vertex_counts[1] * vertex_counts[2];
        const int element_num = 5 * (vertex_counts[0] - 1) * (vertex_counts[1] - 1) * (vertex_counts[2] - 1);
        Matrix3X<T> vertex = MatrixX<T>::Zero(3, vertex_num);
        Matrix4Xi<T> element = Eigen::MatrixXi::Zero(4, element_num);
        int id = 0;
        for (int i = 0;i < vertex_counts[0];i ++)
            for (int j = 0;j < vertex_counts[1];j ++)
                for (int k = 0;k < vertex_counts[2];k++) {
                    vertex(0, id) = i * size;
                    vertex(1, id) = j * size;
                    vertex(2, id) = k * size;
                    id ++;
                }
        id = 0;
        for (int i = 0;i < vertex_counts[0] - 1;i++)
            for (int j = 0;j < vertex_counts[1] - 1;j++)
                for (int k = 0;k < vertex_counts[2] - 1;k++) {
                    std::vector<int> hex_id;
                    hex_id.clear();
                    hex_id.push_back(get_cuboid_id(i, j, k, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i + 1, j, k, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i + 1, j + 1, k, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i, j + 1, k, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i, j, k + 1, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i + 1, j, k + 1, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i + 1, j + 1, k + 1, vertex_counts));
                    hex_id.push_back(get_cuboid_id(i, j + 1, k + 1, vertex_counts));
                    // 0, 2, 1, 5
                    element(0, id) = hex_id[0];
                    element(1, id) = hex_id[2];
                    element(2, id) = hex_id[1];
                    element(3, id) = hex_id[5];
                    id ++;
                    // 0, 4, 7, 5
                    element(0, id) = hex_id[0];
                    element(1, id) = hex_id[4];
                    element(2, id) = hex_id[7];
                    element(3, id) = hex_id[5];
                    id ++;
                    // 0, 2, 5, 7
                    element(0, id) = hex_id[0];
                    element(1, id) = hex_id[2];
                    element(2, id) = hex_id[5];
                    element(3, id) = hex_id[7];
                    id ++;
                    // 2, 6, 5, 7
                    element(0, id) = hex_id[2];
                    element(1, id) = hex_id[6];
                    element(2, id) = hex_id[5];
                    element(3, id) = hex_id[7];
                    id ++;
                    // 0, 7, 3, 2
                    element(0, id) = hex_id[0];
                    element(1, id) = hex_id[7];
                    element(2, id) = hex_id[3];
                    element(3, id) = hex_id[2];
                    id ++;
                }

        return TetrahedralMesh<T>(vertex, element);
    }
}