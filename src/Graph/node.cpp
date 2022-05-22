#include "Graph/node.h"

// ConnectionFace class
ConnectionFace::ConnectionFace(Eigen::Vector3d origin, Eigen::Matrix3d R) {
    _origin = origin;
    _R = R;
}


// Node class
Node::Node() {
    _parent = nullptr;
    _children.clear();
    _T.setIdentity();
    _deformation_type = DeformationType::FFD;
    _type = NodeType::NONE;
    _symbol = "empty";
    _scale_axis.clear();
}

Node::Node(NodeType type, std::string symbol, 
            PartType part_type, std::pair<Vector3d, Vector3d> joint_axis,
            MeshType mesh_type, std::vector<TriMesh> meshes, 
            Node* parent, std::vector<Node*> children,
            ConnectionFace connection_parent,
            std::vector<ConnectionFace> connection_children,
            DeformationType deformation_type,
            std::vector<Eigen::Vector3d> scale_axis,
            WeightHandles handles) {
    _type = type;
    _symbol = symbol;
    _part_type = part_type;
    _joint_axis = joint_axis;
    _mesh_type = mesh_type;
    _meshes = meshes;
    _parent = parent;
    _children.clear();
    for (int i = 0;i < children.size();i++)
        _children.push_back(children[i]);
    _connection_parent = connection_parent;
    _connection_children = connection_children;
    _deformation_type = deformation_type;
    _scale_axis = scale_axis;
    _handles = handles;
    _T.setIdentity();
}

Node::Node(NodeType type, std::string symbol, 
            PartType part_type, std::pair<Vector3d, Vector3d> joint_axis,
            MeshType mesh_type, std::vector<TriMesh> meshes, 
            Node* parent,
            ConnectionFace connection_parent,
            std::vector<ConnectionFace> connection_children,
            DeformationType deformation_type,
            std::vector<Eigen::Vector3d> scale_axis,
            WeightHandles handles) {
    _type = type;
    _symbol = symbol;
    _part_type = part_type;
    _joint_axis = joint_axis;
    _mesh_type = mesh_type;
    _meshes = meshes;
    _parent = parent;
    _children.clear();
    _connection_parent = connection_parent;
    _connection_children = connection_children;
    _deformation_type = deformation_type;
    _scale_axis = scale_axis;
    _handles = handles;
    _T.setIdentity();
}

// Node::Node(const Node &node) {
//     _type = node._type;
//     _symbol = node._symbol;
//     _mesh = node._mesh;
//     _connection_parent = node._connection_parent;
//     _connection_child = node._connection_child;
//     _parent = nullptr;
//     _children.clear();
// }

void Node::copy_from(const Node &node) {
    _type = node._type;
    _symbol = node._symbol;
    _part_type = node._part_type;
    _joint_axis = node._joint_axis;
    _mesh_type = node._mesh_type;
    _meshes = node._meshes;
    _connection_parent = node._connection_parent;
    _connection_children = node._connection_children;
    _parent = nullptr;
    _children.clear();
    _deformation_type = node._deformation_type;
    _scale_axis = node._scale_axis;
    _handles = node._handles;
    _T = node._T;
    _tactile_quad_mesh_verts = node._tactile_quad_mesh_verts;
}

WeightHandles Node::get_transformed_handles() {
    Eigen::Matrix3d R = _T.topLeftCorner(3, 3);
    Eigen::Vector3d p = _T.topRightCorner(3, 1);

    WeightHandles handles = _handles;
    Eigen::Matrix3Xd handle_positions = handles.handle_positions.transpose();
    handle_positions = R * handle_positions;
    
    handle_positions.colwise() += p;
    handles.set_handle_positions(handle_positions.transpose());

    return handles;
}

void Node::get_transformed_meshes(std::vector<Eigen::Matrix3Xd> &V, 
                                    std::vector<Eigen::Matrix3Xi> &F) {
    Eigen::Matrix3d R = _T.topLeftCorner(3, 3);
    Eigen::Vector3d p = _T.topRightCorner(3, 1);
    V.clear(); F.clear();
    for (int i = 0;i < _meshes.size();i++) {
        V.push_back((R * _meshes[i]._V).colwise() + p);
        F.push_back(_meshes[i]._F);
    }
}

std::pair<Vector3d, Vector3d> Node::get_transformed_joint_axis() {
    Eigen::Matrix3d R = _T.topLeftCorner(3, 3);
    Eigen::Vector3d p = _T.topRightCorner(3, 1);
    std::pair<Vector3d, Vector3d> joint_axis;
    joint_axis.first = R * _joint_axis.first + p;
    joint_axis.second = R * _joint_axis.second + p;
    
    return joint_axis;
}

std::vector<Eigen::RowVector3d> Node::get_transformed_scale_axis() {
    Eigen::Matrix3d R = _T.topLeftCorner(3, 3);
    std::vector<Eigen::RowVector3d> scale_axis;
    for (int i = 0;i < _scale_axis.size();i++)
        scale_axis.push_back((R * _scale_axis[i]).transpose());
    return scale_axis;
}

Matrix3X Node::get_transformed_tactile_quad_mesh() {
    Matrix3X transformed_mesh(3, _tactile_quad_mesh_verts.size());
    Eigen::Matrix3d R = _T.topLeftCorner(3, 3);
    Eigen::Vector3d p = _T.topRightCorner(3, 1);
    for (int i = 0; i < _tactile_quad_mesh_verts.size();i++)
        transformed_mesh.col(i) = R * _tactile_quad_mesh_verts[i] + p;
    return transformed_mesh;
}