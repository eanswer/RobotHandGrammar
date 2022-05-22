#include "Grammar/finger_grammar.h"
#include "Mesh/mesh.h"
#include "Graph/node.h"
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include "handles.hpp"
#include "pugixml.hpp"

FingerGrammarRule::FingerGrammarRule(std::string LHS, std::vector<std::string> RHS) {
    _LHS = LHS;
    _RHS = RHS;
}

Eigen::VectorXd FingerGrammar::str_to_eigen(std::string str) {
    std::stringstream iss(str);

    dtype value;
    std::vector<dtype> values;
    while (iss >> value) {
        values.push_back(value);
    }

    VectorXd vec(values.size());
    for (int i = 0;i < values.size();i++)
        vec(i) = values[i];
    
    return vec;
}

FingerGrammar::FingerGrammar() {
    _rules.clear();
}

FingerGrammar::FingerGrammar(std::string grammar_file) {
    ifstream fin(grammar_file);

    std::string folder_dir = grammar_file.substr(0, grammar_file.find_last_of('/') + 1);

    // grammar name
    fin >> _name;

    // load components
    std::string component_xml_path;
    fin >> component_xml_path;
    std::string full_component_xml_path = folder_dir + component_xml_path;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(full_component_xml_path.c_str());

    if (!result) {
        std::cerr << grammar_file << std::endl;
        std::cerr << full_component_xml_path.c_str() << std::endl;
        throw_error("Component xml file does not exist.");
    }

    // std::map<std::string, Node> template_nodes;
    _template_nodes.clear();
    _template_nodes["empty"] = Node();

    for (auto component : doc.children()) {
        if ((std::string)(component.name()) != "component")
            continue;

        std::string symbol = component.attribute("symbol").value();
        
        std::string symbol_type_str = component.attribute("symbol_type").value();
        Node::NodeType symbol_type;
        if (symbol_type_str == "NONTERMINAL" || symbol_type_str == "nonterminal") {
            symbol_type = Node::NodeType::NONTERMINAL;
        } else if (symbol_type_str == "TERMINAL" || symbol_type_str == "terminal") {
            symbol_type = Node::NodeType::TERMINAL;
        } else {
            throw_error("Undefined symbol type.");
        }

        std::string part_type_str = component.attribute("part_type").value();
        Node::PartType part_type;
        if (part_type_str == "Joint" || part_type_str == "joint") {
            part_type = Node::PartType::JOINT;
        } else if (part_type_str == "Body" || part_type_str == "body") {
            part_type = Node::PartType::BODY;
        } else {
            throw_error("Undefined part type.");
        }

        // load deformation type
        Eigen::VectorXd deformation_flags = str_to_eigen(component.child("deformation").attribute("flags").value());
        DeformationType deformation_type;
        std::vector<Eigen::Vector3d> scale_axis;
        if (deformation_flags[0] > 0.5) {
            deformation_type = DeformationType::FFD;
        } else {
            deformation_type = DeformationType::SCALE;
            if (deformation_flags[1] > 0.5) scale_axis.push_back(Eigen::Vector3d::UnitX());
            if (deformation_flags[2] > 0.5) scale_axis.push_back(Eigen::Vector3d::UnitY());
            if (deformation_flags[3] > 0.5) scale_axis.push_back(Eigen::Vector3d::UnitZ());
        }

        // load joint axes
        std::pair<Eigen::Vector3d, Eigen::Vector3d> joint_axis;
        if (part_type == Node::PartType::JOINT) {
            joint_axis.first = str_to_eigen(component.child("joint_axis").attribute("start").value());
            joint_axis.second = str_to_eigen(component.child("joint_axis").attribute("end").value());
        }

        // load meshes
        std::string mesh_type_str = component.child("mesh").attribute("type").value();
        Node::MeshType mesh_type;
        std::vector<TriMesh> meshes;
        if (mesh_type_str == "single" || mesh_type_str == "Single") {
            mesh_type = Node::MeshType::SINGLE;
            std::string stl_path = folder_dir + component.child("mesh").attribute("path").value();

            TriMesh tri_mesh(stl_path, true);
            
            meshes.push_back(tri_mesh);
        } else {
            mesh_type = Node::MeshType::SEPARATE;

            meshes.push_back(TriMesh());
            for (auto parent_mesh : component.child("mesh").children()) {
                if ((std::string)parent_mesh.name() == "parent_mesh") {
                    std::string stl_path = folder_dir + parent_mesh.attribute("path").value();
                    TriMesh tri_mesh(stl_path, true);
                    meshes[0].merge(tri_mesh._V, tri_mesh._F);
                }
            }

            meshes.push_back(TriMesh());
            for (auto parent_mesh : component.child("mesh").children()) {
                if ((std::string)parent_mesh.name() == "child_mesh") {
                    std::string stl_path = folder_dir + parent_mesh.attribute("path").value();
                    TriMesh tri_mesh(stl_path, true);
                    meshes[1].merge(tri_mesh._V, tri_mesh._F);
                }
            }

            for (auto parent_mesh : component.child("mesh").children()) {
                if ((std::string)parent_mesh.name() == "aux_mesh") {
                    std::string stl_path = folder_dir + parent_mesh.attribute("path").value();
                    TriMesh tri_mesh(stl_path, true);
                    meshes.push_back(tri_mesh);
                }
            }
        }

        // load cage handles
        std::string handle_path = component.child("cage").attribute("path").value();
        WeightHandles handles;
        std::string full_handle_path = folder_dir + handle_path;
        handles.load_handle_file(full_handle_path);
        
        // read parent connection faces
        Eigen::Vector3d p0 = str_to_eigen(component.child("connection").child("parent_connection").attribute("origin").value());
        Eigen::Matrix3d R0;
        R0.col(0) = str_to_eigen(component.child("connection").child("parent_connection").attribute("x_axis").value());
        R0.col(1) = str_to_eigen(component.child("connection").child("parent_connection").attribute("y_axis").value());
        R0.col(2) = str_to_eigen(component.child("connection").child("parent_connection").attribute("z_axis").value());
        ConnectionFace connection_parent(p0, R0);

        int num_children = component.child("connection").child("child_connection").attribute("count").as_int();
        std::vector<ConnectionFace> connection_children = std::vector<ConnectionFace>(num_children, ConnectionFace());
        for (auto connection_face : component.child("connection").child("child_connection").children()) {
            int id = connection_face.attribute("id").as_int();
            Eigen::Vector3d p1 = str_to_eigen(connection_face.attribute("origin").value());
            Eigen::Matrix3d R1;
            R1.col(0) = str_to_eigen(connection_face.attribute("x_axis").value()).normalized();
            R1.col(1) = str_to_eigen(connection_face.attribute("y_axis").value()).normalized();
            R1.col(2) = str_to_eigen(connection_face.attribute("z_axis").value()).normalized();
            connection_children[id] = ConnectionFace(p1, R1);
        }

        _template_nodes[symbol] = Node(symbol_type, symbol,
                                        part_type, joint_axis,
                                        mesh_type, meshes,
                                        nullptr, 
                                        connection_parent, connection_children,
                                        deformation_type, scale_axis,
                                        handles);

        // load tactile quad mesh points
        _template_nodes[symbol]._tactile_quad_mesh_verts.clear();
        if (component.child("tactile_quad_mesh")) {
            std::string tactile_quad_mesh_path = component.child("tactile_quad_mesh").attribute("path").value();
            std::string full_tactile_quad_mesh_path = folder_dir + tactile_quad_mesh_path;
            FILE* fp = fopen(full_tactile_quad_mesh_path.c_str(), "r");
            int N_verts;
            fscanf(fp, "%d", &N_verts);
            if (N_verts != 8) {
                throw_error("Only support tactile quad mesh with 8 vertices now.");
            }
            for (int i = 0;i < N_verts;i++) {
                double x, y, z;
                fscanf(fp, "%lf %lf %lf", &x, &y, &z);
                _template_nodes[symbol]._tactile_quad_mesh_verts.push_back(Vector3(x, y, z));
            }
            fclose(fp);
            std::vector<Vector3> project;
            project.clear();
            for (int i = 0;i < N_verts;i++) {
                project.push_back(Vector3::Zero());
                project[i](0) = (_template_nodes[symbol]._tactile_quad_mesh_verts[i] - connection_parent._origin).dot(-connection_parent._R.col(2));
                project[i](1) = (_template_nodes[symbol]._tactile_quad_mesh_verts[i] - connection_parent._origin).dot(connection_parent._R.col(0));
                project[i](2) = (_template_nodes[symbol]._tactile_quad_mesh_verts[i] - connection_parent._origin).dot(connection_parent._R.col(1));
            }
            for (int i = 0;i < N_verts;i++)
                for (int j = i;j < N_verts;j++)
                    if ((project[j](0) < project[i](0) - 1e-7) ||
                        (fabs(project[j](0) - project[i](0)) < 1e-7 && project[j](1) < project[i](1) - 1e-7) ||
                        (fabs(project[j](0) - project[i](0)) < 1e-7 && fabs(project[j](1) - project[i](1)) < 1e-7 && project[j](2) < project[i](2) - 1e-7)) {
                            swap(project[i], project[j]);
                            swap(_template_nodes[symbol]._tactile_quad_mesh_verts[i], _template_nodes[symbol]._tactile_quad_mesh_verts[j]);
                        }
        }
    }

    _rules.clear();

    int m; fin >> m;   
    for (int i = 0;i < m;i++) {
        std::string str;
        fin >> str; // read '{'
        std::string LHS;
        fin >> LHS;
        fin >> str; // read "->"
        std::vector<std::string> RHS;
        for (;;) {
            fin >> str;
            if (str != "}") {
                // Node rhs_node;
                // rhs_node.copy_from(template_nodes[str]);
                RHS.push_back(str);
            } else {
                break;
            }
        }
        
        _rules.push_back(FingerGrammarRule(LHS, RHS));
    }
}

void FingerGrammar::add_rule(FingerGrammarRule rule) {
    _rules.push_back(rule);
}

void FingerGrammar::reset() {
    _graph = get_init_graph();
    _graph_his.clear();
    _graph_his.push_back(_graph);
}

void FingerGrammar::undo() {
    if (_graph_his.size() > 1) {
        _graph_his.pop_back();
        _graph = _graph_his[_graph_his.size() - 1];
    }
}

Tree FingerGrammar::get_init_graph() {
    Node* node = new Node();
    node->copy_from(_template_nodes[_starting_symbol]);

    Tree graph(node);
    return graph;
}

std::vector<int> FingerGrammar::get_available_rules() {
    std::vector<int> available_rules;
    available_rules.clear();
    for (int i = 0;i < _rules.size();i++) {
        Node* res = _graph.find_node(_rules[i]);
        if (res != nullptr) {
            available_rules.push_back(i);
        }
    }
    return available_rules;
}

bool FingerGrammar::apply_rule(int rule_id) {   
    Tree new_tree;
    new_tree.copy_from(_graph);
    Node* application_node = new_tree.find_node(_rules[rule_id]);
    if (application_node != nullptr) {
        // create the chain of the rhs
        std::vector<Node*> chain;
        chain.clear();
        for (int k = 0;k < _rules[rule_id]._RHS.size();k++) {
            Node* node = new Node();
            node->copy_from(_rules[rule_id]._RHS_nodes[k]);

            if (k > 0) {
                node->_parent = chain[k - 1];
                chain[k - 1]->_children.push_back(node);
            }
            chain.push_back(node);
        }

        // replace the symbol in the tree by the chain
        new_tree.replace(application_node, chain);

        _graph = new_tree;
        _graph_his.push_back(_graph);
        return true;
    }
    else {
        return false;
    }
}

void FingerGrammar::finalize_rule_nodes() {
    for (int i = 0;i < _rules.size();i++) {
        _rules[i]._RHS_nodes.clear();
        for (int j = 0;j < _rules[i]._RHS.size();j++) {
            Node rhs_node;
            rhs_node.copy_from(_template_nodes[_rules[i]._RHS[j]]);
            _rules[i]._RHS_nodes.push_back(rhs_node);
        }
    }
}

void FingerGrammar::print_grammar() {
    std::cerr << "------------- Grammar ------------------" << std::endl;
    for (int i = 0;i < _rules.size();i++) {
        std::cerr << _rules[i]._LHS << " ->";
        for (int j = 0;j < _rules[i]._RHS.size();j++)
            std::cerr << ' ' << _rules[i]._RHS[j];
        std::cerr << std::endl;
    }
    std::cerr << "----------------------------------------" << std::endl;
}
