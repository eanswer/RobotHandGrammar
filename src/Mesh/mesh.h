/**
 * @file mesh.h
 * @author Jie Xu 
 */ 
#pragma once

#include "Common/Common.h"

/**
 * Implement an abstract class for mesh. 
 */
class Mesh {
public:
    Eigen::Matrix3Xd _V;
};

/**
 * Implement a Triangle Mesh class. 
 */
class TriMesh : public Mesh {
public:
    TriMesh() {}
    TriMesh(Eigen::Matrix3Xd V, Eigen::Matrix3Xi F);
    TriMesh(std::string stl_path, bool abs_path = false);
    
    void merge(TriMesh& mesh);
    void merge(Eigen::Matrix3Xd V, Eigen::Matrix3Xi F);
    void rotate(Vector3d axis, dtype angle);
    void translate(Vector3d translation);

    void export_obj(std::string filename);
    
    Eigen::Matrix3Xd _V;
    Eigen::Matrix3Xi _F;
};

/**
 * Implement a Quad Mesh class.
 */
class QuadMesh : public Mesh {
public:
    QuadMesh() {}
    QuadMesh(Eigen::Matrix3Xd V, Eigen::Matrix4Xi F);

    void export_obj(std::string filename);

    Eigen::Matrix3Xd _V;
    Eigen::Matrix4Xi _F;
};