/**
 * @file node.h
 * @author Jie Xu 
 */ 
#pragma once

#include "Common/Common.h"
#include "Mesh/mesh.h"
#include "handles.hpp"

/**
 * Implement a face connection constraints.
 * 
 * Each connection face is a face supposed to be aligned by two objects.
 * the contains rotation matrix and the translation matrix represented in
 * the body frame.
 */
class ConnectionFace {
public:
    ConnectionFace() {}
    ConnectionFace(Eigen::Vector3d origin, Eigen::Matrix3d R);

    Eigen::Vector3d _origin;
    Eigen::Matrix3d _R;
};

enum DeformationType {
    FFD = 0,
    SCALE = 1
};

/**
 * Implement a node in design graph (tree).
 * 
 * Each node can be terminal or nonterminal, associated with a mesh, 
 * a transformation matrix relative to parent frame, and a deformation cage.
 */
class Node {
public:
    enum NodeType {NONTERMINAL, TERMINAL, NONE};
    enum PartType {BODY, JOINT};
    enum MeshType {SINGLE, SEPARATE};
    
    Node();
    Node(NodeType type, std::string symbol, 
            PartType part_type, std::pair<Vector3d, Vector3d> joint_axis,
            MeshType mesh_type, std::vector<TriMesh> meshes, 
            Node* parent, std::vector<Node*> children,
            ConnectionFace connection_parent,
            std::vector<ConnectionFace> connection_children,
            DeformationType deformation_type,
            std::vector<Eigen::Vector3d> scale_axis,
            WeightHandles handles);
    Node(NodeType type, std::string symbol, 
            PartType part_type, std::pair<Vector3d, Vector3d> joint_axis,
            MeshType mesh_type, std::vector<TriMesh> meshes, 
            Node* parent,
            ConnectionFace connection_parent,
            std::vector<ConnectionFace> connection_children,
            DeformationType deformation_type,
            std::vector<Eigen::Vector3d> scale_axis,
            WeightHandles handles);
    // Node(const Node &node);

    void copy_from(const Node &node);

    WeightHandles get_transformed_handles();
    void get_transformed_meshes(std::vector<Eigen::Matrix3Xd> &V, 
                                std::vector<Eigen::Matrix3Xi> &F);
    std::pair<Vector3d, Vector3d> get_transformed_joint_axis();
    std::vector<Eigen::RowVector3d> get_transformed_scale_axis();
    Matrix3X get_transformed_tactile_quad_mesh();

    NodeType _type;
    std::string _symbol;
    PartType _part_type;
    std::pair<Vector3d, Vector3d> _joint_axis;
    MeshType _mesh_type;
    std::vector<TriMesh> _meshes;
    ConnectionFace _connection_parent;
    std::vector<ConnectionFace> _connection_children;
    DeformationType _deformation_type;
    std::vector<Eigen::Vector3d> _scale_axis;
    WeightHandles _handles;
    std::vector<Vector3> _tactile_quad_mesh_verts;

    Eigen::Matrix4d _T;

    Node* _parent;          // parent in the graph
    std::vector<Node*> _children;

    int _id_in_objects;
};