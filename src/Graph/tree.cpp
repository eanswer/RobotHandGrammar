#include "Graph/tree.h"
#include "Mesh/mesh.h"
#include "linear_blend_skinning_object.h"
#include <igl/opengl/glfw/Viewer.h>
#include "Grammar/finger_grammar.h"
#include "Grammar/palm_grammar.h"

Tree::Tree() {
    _root = nullptr;
}

Tree::Tree(Node* root) {
    _root = root;
}

void Tree::copy_from(const Tree &tree) {
    _root = copy(tree._root);
}

Node* Tree::copy(const Node* root) {
    Node *node = new Node();
    node->copy_from(*root);
    for (int i = 0;i < root->_children.size();i++) {
        node->_children.push_back(copy(root->_children[i]));
        node->_children[i]->_parent = node;
    }
    return node;
}

Node* Tree::find_node(const FingerGrammarRule rule) const {
    if (_root != nullptr)
        return find_node(_root, rule);
    else
        return nullptr;
}

Node* Tree::find_node(Node* node, const FingerGrammarRule rule) const {
    if (node->_symbol == rule._LHS) {
        int num_children_after = node->_children.size() + (rule._RHS.size() > 1 ? 1 : 0);
        if (num_children_after <= rule._RHS_nodes[0]._connection_children.size()) {
            return node;
        }
    }

    for (int i = 0;i < node->_children.size();i++) {
        Node* res = find_node(node->_children[i], rule);
        if (res != nullptr) {
            return res;
        }
    }

    return nullptr;
}

// // type 1 replace: if original tree is like A->B, after replacing A by a->c, it becomes to a->c->B
// void Tree::replace(int node_id, std::vector<Node*> chain) {

//     Node* deleted_node = _nodes[node_id];
//     Node* parent = deleted_node->_parent;
//     std::vector<Node*> &children = deleted_node->_children;

//     if (chain.size() == 0) {
//         _nodes.erase(_nodes.begin() + node_id);
//         if (parent != nullptr) {
//             for (int i = 0;i < parent->_children.size();i++)
//                 if (parent->_children[i] == deleted_node) {
//                     parent->_children.erase(parent->_children.begin() + i);
//                     break;
//                 }
//             for (int i = 0;i < children.size();i++)
//                 parent->_children.push_back(children[i]);
//         }
//         for (int i = 0;i < children.size();i++)
//             children[i]->_parent = parent;
//     } else {
//         _nodes.resize(_nodes.size() + chain.size() - 1);
//         for (int i = _nodes.size() - 1;i >= node_id + chain.size();i--) {
//             _nodes[i] = _nodes[i - chain.size() + 1];
//         }
//         for (int i = node_id;i < node_id + chain.size();i++)
//             _nodes[i] = chain[i - node_id];

//         // connect head and parent
//         chain[0]->_parent = parent;
//         if (parent != nullptr) {
//             for (int i = 0;i < parent->_children.size();i++)
//                 if (parent->_children[i] == deleted_node) {
//                     parent->_children[i] = chain[0];
//                 }
//         }

//         // connect tail and children
//         Node* tail = chain[chain.size() - 1];
//         for (int i = 0;i < children.size();i++) {
//             children[i]->_parent = tail;
//             tail->_children.push_back(children[i]);
//         }
//     }

//     if (deleted_node == _root) {
//         _root = chain[0];
//     }

//     deleted_node->_parent = nullptr;
//     deleted_node->_children.clear();

//     delete deleted_node;
// }

// type 2 replace: if original tree is like A->B, after replacing A by a->c, it becomes to a->B
//                                                                                         |
//                                                                                         c
void Tree::replace(Node* node, std::vector<Node*> chain) {
    Node* deleted_node = node;
    Node* parent = deleted_node->_parent;
    std::vector<Node*> &children = deleted_node->_children;

    if (chain.size() == 0) { // very special case, be careful
        if (parent != nullptr) {
            for (int i = 0;i < parent->_children.size();i++)
                if (parent->_children[i] == deleted_node) {
                    parent->_children.erase(parent->_children.begin() + i);
                    break;
                }
            for (int i = 0;i < children.size();i++)
                parent->_children.push_back(children[i]);
        }
        for (int i = 0;i < children.size();i++)
            children[i]->_parent = parent;
    } else {
        // connect to parent
        chain[0]->_parent = parent;
        if (parent != nullptr) {
            for (int i = 0;i < parent->_children.size();i++)
                if (parent->_children[i] == deleted_node) {
                    parent->_children[i] = chain[0];
                }
        }

        // connect to children
        chain[0]->_children.clear();
        for (int i = 0;i < children.size();i++) {
            chain[0]->_children.push_back(children[i]);
            children[i]->_parent = chain[0];
        }
        
        if (chain.size() > 1) {
            chain[0]->_children.push_back(chain[1]);
        }
    }

    if (deleted_node == _root) {
        _root = chain[0];
    }

    deleted_node->_parent = nullptr;
    deleted_node->_children.clear();
    
    delete deleted_node;
}

void Tree::get_deformation_objects(std::vector<LinearBlendSkinningObject*> &objects,
                        std::vector<Node::NodeType> &types,
                        PalmGrammar* palm_grammar) {
    objects.clear();
    types.clear();
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topRightCorner(3, 1) = _root->_connection_parent._origin;
    get_deformation_objects(_root, -1, T, objects, types, palm_grammar);
}


void Tree::get_deformation_objects(
        Node* root, 
        int parent_id,
        Eigen::Matrix4d T,
        std::vector<LinearBlendSkinningObject*> &objects,
        std::vector<Node::NodeType> &types,
        PalmGrammar* palm_grammar) {

    Eigen::Matrix3d R = T.topLeftCorner(3, 3);
    Eigen::Vector3d p = T.topRightCorner(3, 1);

    root->_T = T;
    
    // special process for palm, supposed to be root
    if (root->_symbol == "P" || root->_symbol == "palm") {
        std::vector<TriMesh> palm_meshes;
        std::vector<Vector3> colors;
        std::vector<TriMesh> parent_meshes;
        std::vector<TriMesh> child_meshes;
        std::vector<bool> connector_masks;
        std::vector<WeightHandles> handles_list;
        std::vector<std::pair<int, int> > cell_coord;
        palm_grammar->construct_palm_meshes(palm_meshes, colors, parent_meshes, child_meshes, connector_masks, handles_list, cell_coord);

        for (int i = 0;i < palm_meshes.size();i++) {
            LinearBlendSkinningObject* object = new LinearBlendSkinningObject();
            object->set_mesh_type(Node::MeshType::SINGLE);
            std::vector<Eigen::Matrix3Xd> V;
            std::vector<Eigen::Matrix3Xi> F;
            V.push_back(palm_meshes[i]._V);
            F.push_back(palm_meshes[i]._F);
            object->set_mesh(V, F);
            object->set_part_type(Node::PartType::BODY);
            object->set_joint_axis(std::make_pair(Vector3::Zero(), Vector3::Zero()));
            
            object->_parent_mesh = parent_meshes[i];
            object->_child_mesh = child_meshes[i];

            object->_parent = nullptr;

            object->_symbol = root->_symbol;

            object->_cell_row = cell_coord[i].first;
            object->_cell_col = cell_coord[i].second;

            object->_handles = handles_list[i];
            // object->compute_weights();
            
            objects.push_back(object);
            types.push_back(root->_type);
        }

        int connector_id = 0;
        for (int i = 0;i < root->_children.size();i++) {
            while (!connector_masks[connector_id]) connector_id ++;

            Node* child = root->_children[i];
            if (child->_type == Node::NodeType::NONE) {
                continue;
            }
            
            Eigen::Matrix3d R_pc;
            Eigen::Vector3d p_pc;
            R_pc = root->_connection_children[i]._R * child->_connection_parent._R.transpose();
            p_pc = root->_connection_children[i]._origin - 
                        R_pc * child->_connection_parent._origin;

            Eigen::Matrix4d T_pc = Eigen::Matrix4d::Identity();
            T_pc.topLeftCorner(3, 3) = R_pc;
            T_pc.topRightCorner(3, 1) = p_pc;

            Eigen::Matrix4d T_child = T * T_pc;

            get_deformation_objects(root->_children[i], connector_id, T_child, objects, types, palm_grammar);

            connector_id += 1;
        }
    } else {

        LinearBlendSkinningObject* object = new LinearBlendSkinningObject();

        // get transformed mesh
        std::vector<Eigen::Matrix3Xd> V;
        std::vector<Eigen::Matrix3Xi> F;
        root->get_transformed_meshes(V, F);

        object->set_mesh_type(root->_mesh_type);
        object->set_mesh(V, F);
        
        object->set_part_type(root->_part_type);
        object->set_joint_axis(root->get_transformed_joint_axis());

        // get transformed handle positions
        WeightHandles handles = root->get_transformed_handles();
        object->_handles = handles;

        object->_tactile_quad_mesh_V = root->get_transformed_tactile_quad_mesh().transpose();
        // object->_tactile_quad_mesh_V = MatrixX3d(root->_tactile_quad_mesh_verts.size(), 3);
        // for (int i = 0;i < root->_tactile_quad_mesh_verts.size();i++) {
        //     object->_tactile_quad_mesh_V.row(i) = root->_tactile_quad_mesh_verts[i].transpose();
        // }

        // object->compute_weights();

        object->_symbol = root->_symbol;

        // get interface constraints
        if (parent_id != -1) {
            object->_parent = objects[parent_id];
            object->_parent->_children.push_back(object);
            object->_parent->_interface_vertices_children.push_back(std::vector<std::pair<int, int> >());
            std::vector<std::pair<int, int> > &parent_interface_vertices_children
                = object->_parent->_interface_vertices_children[object->_parent->_interface_vertices_children.size() - 1]; 

            WeightHandles parent_handles = object->_parent->_handles;
            for (int i = 0;i < parent_handles.handle_positions.rows();i++)
                for (int j = 0;j < handles.handle_positions.rows();j++)
                    if ((parent_handles.handle_positions.row(i) - 
                            handles.handle_positions.row(j)).norm() < 1e-5) {
                        object->_interface_vertices_parent.push_back(std::make_pair(j, i));
                        parent_interface_vertices_children.push_back(std::make_pair(i, j));
                        break;
                    }
        } else {
            object->_parent = nullptr;
        }

        // set deformation type
        object->_deformation_type = root->_deformation_type;

        // set scale origin
        object->_scale_origin = (R * root->_connection_parent._origin + p).transpose();
        object->_scale_origin_old = object->_scale_origin;
        
        // set scale axis
        object->_scale_axis = root->get_transformed_scale_axis();

        // set origin and directions
        object->_parent_face_origin = object->_scale_origin;
        object->_direction = -R * root->_connection_parent._R.col(2);
        object->_direction_x = R * root->_connection_parent._R.col(0);
        object->_direction_y = R * root->_connection_parent._R.col(1);

        root->_id_in_objects = objects.size();
        
        int object_id = objects.size();

        // TODO: remove, fix for palm
        // find handle correspondence

        objects.push_back(object);
        types.push_back(root->_type);
        
        for (int i = 0;i < root->_children.size();i++) {
            Node* child = root->_children[i];
            if (child->_type == Node::NodeType::NONE) {
                continue;
            }

            Eigen::Matrix3d R_pc;
            Eigen::Vector3d p_pc;
            // R_pc = child->_connection_parent._R.transpose() * 
            //             root->_connection_children[i]._R;
            R_pc = root->_connection_children[i]._R * child->_connection_parent._R.transpose();
            p_pc = root->_connection_children[i]._origin - 
                        R_pc * child->_connection_parent._origin;

            Eigen::Matrix4d T_pc = Eigen::Matrix4d::Identity();
            T_pc.topLeftCorner(3, 3) = R_pc;
            T_pc.topRightCorner(3, 1) = p_pc;

            Eigen::Matrix4d T_child = T * T_pc;

            get_deformation_objects(child, object_id, T_child, objects, types, palm_grammar);
        }
    }
}

void Tree::get_objects_and_constraints(std::vector<LinearBlendSkinningObject*> &objects,
                        std::vector<Node::NodeType> &types,
                        std::map<std::pair<int, int>, std::pair<int, int> > &constraints) {
    objects.clear();
    types.clear();
    constraints.clear();
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topRightCorner(3, 1) = _root->_connection_parent._origin;
    get_objects_and_constraints(_root, -1, T, objects, types, constraints);
}

void Tree::get_objects_and_constraints(
        Node* root, 
        int parent_id,
        Eigen::Matrix4d T,
        std::vector<LinearBlendSkinningObject*> &objects,
        std::vector<Node::NodeType> &types,
        std::map<std::pair<int, int>, std::pair<int, int> > &constraints) {

    Eigen::Matrix3d R = T.topLeftCorner(3, 3);
    Eigen::Vector3d p = T.topRightCorner(3, 1);

    root->_T = T;
    
    LinearBlendSkinningObject* object = new LinearBlendSkinningObject();

    // get transformed mesh
    std::vector<Eigen::Matrix3Xd> V;
    std::vector<Eigen::Matrix3Xi> F;
    root->get_transformed_meshes(V, F);

    object->set_mesh_type(root->_mesh_type);
    object->set_mesh(V, F);
    
    object->set_part_type(root->_part_type);
    object->set_joint_axis(root->get_transformed_joint_axis());

    // get transformed handle positions
    WeightHandles handles = root->get_transformed_handles();
    object->_handles = handles;
    object->compute_weights();

    // get interface constraints
    if (parent_id != -1) {
        object->_parent = objects[parent_id];
        object->_parent->_children.push_back(object);
        object->_parent->_interface_vertices_children.push_back(std::vector<std::pair<int, int> >());
        std::vector<std::pair<int, int> > &parent_interface_vertices_children
            = object->_parent->_interface_vertices_children[object->_parent->_interface_vertices_children.size() - 1]; 

        WeightHandles parent_handles = object->_parent->_handles;
        for (int i = 0;i < parent_handles.handle_positions.rows();i++)
            for (int j = 0;j < handles.handle_positions.rows();j++)
                if ((parent_handles.handle_positions.row(i) - 
                        handles.handle_positions.row(j)).norm() < 1e-5) {
                    object->_interface_vertices_parent.push_back(std::make_pair(j, i));
                    parent_interface_vertices_children.push_back(std::make_pair(i, j));
                    break;
                }
    } else {
        object->_parent = nullptr;
    }

    // set deformation type
    object->_deformation_type = root->_deformation_type;

    // set scale origin
    object->_scale_origin = (R * root->_connection_parent._origin + p).transpose();
    object->_scale_origin_old = object->_scale_origin;
    
    // set scale axis
    object->_scale_axis = root->get_transformed_scale_axis();

    int object_id = objects.size();

    // TODO: remove, fix for palm
    // find handle correspondence
    if (parent_id != -1) {
        // TODO: accelerate by specifying in file or using data structure.
        WeightHandles parent_handles = root->_parent->get_transformed_handles();
        for (int i = 0;i < parent_handles.handle_positions.rows();i++)
            for (int j = 0;j < handles.handle_positions.rows();j++) {
                if ((parent_handles.handle_positions.row(i) - 
                        handles.handle_positions.row(j)).norm() < 1e-5) {
                    constraints[std::make_pair(parent_id, i)] = std::make_pair(object_id, j);
                    constraints[std::make_pair(object_id, j)] = std::make_pair(parent_id, i);
                    break;
                }
            }
    }

    objects.push_back(object);
    types.push_back(root->_type);
    
    for (int i = 0;i < root->_children.size();i++) {
        Node* child = root->_children[i];
        if (child->_type == Node::NodeType::NONE) {
            continue;
        }

        Eigen::Matrix3d R_pc;
        Eigen::Vector3d p_pc;
        // R_pc = child->_connection_parent._R.transpose() * 
        //             root->_connection_children[i]._R;
        R_pc = root->_connection_children[i]._R * child->_connection_parent._R.transpose();
        p_pc = root->_connection_children[i]._origin - 
                    R_pc * child->_connection_parent._origin;

        Eigen::Matrix4d T_pc = Eigen::Matrix4d::Identity();
        T_pc.topLeftCorner(3, 3) = R_pc;
        T_pc.topRightCorner(3, 1) = p_pc;

        Eigen::Matrix4d T_child = T * T_pc;

        get_objects_and_constraints(child, object_id, T_child, objects, types, constraints);
    }
}

void Tree::print_graph() {
    std::cerr << "Graph: " << std::endl;
    std::cerr << "root = " << _root << std::endl;
    print_graph(_root);
    std::cerr << std::endl;
}

void Tree::print_graph(Node* node) {
    if (node->_parent == nullptr) {
        std::cerr << node->_symbol << "(" << node << "), root" << std::endl;
    } else {
        std::cerr << node->_symbol << "(" << node << "), parent = " << node->_parent << std::endl;
    }
    std::cerr << "children: ";
    for (int i = 0;i < node->_children.size();i++)
        std::cerr << node->_children[i] << " ";
    std::cerr << std::endl;
    for (int i = 0;i < node->_children.size();i++)
        print_graph(node->_children[i]);
}