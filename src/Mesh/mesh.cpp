#include "Mesh/mesh.h"
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>

TriMesh::TriMesh(Eigen::Matrix3Xd V, Eigen::Matrix3Xi F) {
    _V = V;
    _F = F;
}

TriMesh::TriMesh(std::string stl_path, bool abs_path) {
    Eigen::MatrixXd V, V_unique;
    Eigen::MatrixXi F;
    Eigen::VectorXi SVI, SVJ;
    std::string full_stl_path;
    if (abs_path) {
        full_stl_path = stl_path;
    } else {
        full_stl_path = PROJECT_SOURCE_DIR + stl_path;
    }
    
    bool flag = igl::read_triangle_mesh(full_stl_path, V, F);
    if (!flag) return;
    igl::remove_duplicate_vertices(V, 0, V_unique, SVI, SVJ);
                std::for_each(F.data(), F.data() + F.size(), [&SVJ](int& f) {f = SVJ(f); });
                V = V_unique;

    _V = V.transpose();
    _F = F.transpose();
}

void TriMesh::merge(TriMesh& mesh) {
    merge(mesh._V, mesh._F);
}

void TriMesh::merge(Eigen::Matrix3Xd V, Eigen::Matrix3Xi F) {
    int num_vertices_old = _V.cols();
    int num_faces_old = _F.cols();
    _V.conservativeResize(3, num_vertices_old + V.cols());
    _F.conservativeResize(3, num_faces_old + F.cols());
    _V.block(0, num_vertices_old, 3, V.cols()) = V;
    _F.block(0, num_faces_old, 3, F.cols()) = F + Matrix3Xi::Ones(3, F.cols()) * num_vertices_old;
}

void TriMesh::rotate(Vector3d axis, dtype angle) {
    Matrix3d R = Eigen::AngleAxis<dtype>(angle, axis).matrix();
    _V = R * _V;
}

void TriMesh::translate(Vector3d translation) {
    _V.colwise() += translation;
}

void TriMesh::export_obj(std::string filename) {
    FILE* fp = fopen(filename.c_str(), "w");
    for (int i = 0;i < _V.cols();i++) {
        fprintf(fp, "v %.5lf %.5lf %.5lf\n", _V(0, i), _V(1, i), _V(2, i));
    }
    for (int i = 0;i < _F.cols();i++) {
        fprintf(fp, "f %d %d %d\n", _F(0, i) + 1, _F(1, i) + 1, _F(2, i) + 1);
    }
    fclose(fp);
}

/********************** Quad Mesh ************************/
QuadMesh::QuadMesh(Eigen::Matrix3Xd V, Eigen::Matrix4Xi F) {
    _V = V;
    _F = F;
}

void QuadMesh::export_obj(std::string filename) {
    FILE* fp = fopen(filename.c_str(), "w");
    for (int i = 0;i < _V.cols();i++) {
        fprintf(fp, "v %.5lf %.5lf %.5lf\n", _V(0, i), _V(1, i), _V(2, i));
    }
    for (int i = 0;i < _F.cols();i++) {
        fprintf(fp, "f %d %d %d %d\n", _F(0, i) + 1, _F(1, i) + 1, _F(2, i) + 1, _F(3, i) + 1);
    }
    fclose(fp);
}