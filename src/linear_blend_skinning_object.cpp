#include "linear_blend_skinning_object.h"
#include "nearest_neighbor_weights.hpp"
#include "bounded_biharmonic_weights.hpp"
#include "mean_value_weights.h"
#include "linear_weights.hpp"
#include "Grammar/palm_grammar.h"

LinearBlendSkinningObject::LinearBlendSkinningObject() {
    _weight_type = 0;
    _stale = true;
    _parent = nullptr;
    _children.clear();
    _interface_vertices_parent.clear();
    _interface_vertices_children.clear();
    _lbs_mat.clear();
}

LinearBlendSkinningObject::~LinearBlendSkinningObject() {
    _parent = nullptr;
    _children.clear();
    _interface_vertices_parent.clear();
    _interface_vertices_children.clear();
}

void LinearBlendSkinningObject::set_mesh_type(Node::MeshType mesh_type) {
    _mesh_type = mesh_type;
}

void LinearBlendSkinningObject::set_mesh(
    std::vector<Matrix3Xd> V, 
    std::vector<Matrix3Xi> F) {
    _V.clear(); _F.clear();
    for (int i = 0;i < V.size();i++) {
        _V.push_back(V[i].transpose());
        _F.push_back(F[i].transpose());
    }
    _stale = true;
}

void LinearBlendSkinningObject::set_mesh(
    std::vector<MatrixX3d> V, 
    std::vector<MatrixX3i> F) {
    _V = V;
    _F = F;
    _stale = true;
}

void LinearBlendSkinningObject::set_part_type(Node::PartType part_type) {
    _part_type = part_type;
}

void LinearBlendSkinningObject::set_joint_axis(std::pair<Vector3d, Vector3d> joint_axis) {
    _joint_axis = joint_axis;
}

void LinearBlendSkinningObject::set_weight_type(int weight_type) {
    _weight_type = weight_type;
    _stale = true;
}

std::vector<MatrixX3d> LinearBlendSkinningObject::transform_mesh(unsigned int naive_transformation_mode) {
    if (_stale || _handles._stale) {
        _transformed_V.clear();
        if (naive_transformation_mode == 1) {
            // infer the face center point from corner points
            for (int i = 0;i < _handles.cage_faces.size();i++) {
                int center_idx = 8 + i;
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for (int j = 0;j < _handles.cage_faces[i].size();j++) {
                    Eigen::Vector3d pos = _handles.handle_positions.row(_handles.cage_faces[i][j]) 
                                        + _handles.handle_transforms.block(4 * _handles.cage_faces[i][j] + 3, 0, 1, 3);
                    centroid += pos;
                }
                centroid = centroid / _handles.cage_faces[i].size();
                _handles.move_handle(center_idx, centroid);
            }
        }
        for (int i = 0;i < _lbs_mat.size();i++) {
            if (_handles.positions().rows() > 0) {
                _transformed_V.push_back(_lbs_mat[i] * _handles.transform());
            } else {
                _transformed_V.push_back(_V[i]);
            }
        }
        _stale = false;
        _handles._stale = false;
    }
    return _transformed_V;
}

std::pair<RowVectorXd, RowVectorXd> LinearBlendSkinningObject::transform_joint_axis() {
    std::pair<RowVectorXd, RowVectorXd> transformed_joint_axis;
    if (_handles.positions().rows() > 0) {
        transformed_joint_axis.first = _lbs_mat_joint_axis.first * _handles.transform();
        transformed_joint_axis.second = _lbs_mat_joint_axis.second * _handles.transform();
    } else {
        transformed_joint_axis = _joint_axis;
    }
    return transformed_joint_axis;
}

MatrixX3d LinearBlendSkinningObject::transform_tactile_mesh(unsigned int naive_transformation_mode) {
    if (naive_transformation_mode == 1) {
        // infer the face center point from corner points
        for (int i = 0;i < _handles.cage_faces.size();i++) {
            int center_idx = 8 + i;
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (int j = 0;j < _handles.cage_faces[i].size();j++) {
                Eigen::Vector3d pos = _handles.handle_positions.row(_handles.cage_faces[i][j]) 
                                    + _handles.handle_transforms.block(4 * _handles.cage_faces[i][j] + 3, 0, 1, 3);
                centroid += pos;
            }
            centroid = centroid / _handles.cage_faces[i].size();
            _handles.move_handle(center_idx, centroid);
        }
    }
    MatrixX3d transformed_tactile_mesh = _lbs_mat_tactile_mesh * _handles.transform();

    return transformed_tactile_mesh;
}

MatrixX3d LinearBlendSkinningObject::transform_parent_mesh(unsigned int naive_transformation_mode) {
    if (naive_transformation_mode == 1) {
        // infer the face center point from corner points
        for (int i = 0;i < _handles.cage_faces.size();i++) {
            int center_idx = 8 + i;
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (int j = 0;j < _handles.cage_faces[i].size();j++) {
                Eigen::Vector3d pos = _handles.handle_positions.row(_handles.cage_faces[i][j]) 
                                    + _handles.handle_transforms.block(4 * _handles.cage_faces[i][j] + 3, 0, 1, 3);
                centroid += pos;
            }
            centroid = centroid / _handles.cage_faces[i].size();
            _handles.move_handle(center_idx, centroid);
        }
    }
    MatrixX3d transformed_parent_mesh = _lbs_mat_parent_mesh * _handles.transform();

    return transformed_parent_mesh;
}

MatrixX3d LinearBlendSkinningObject::transform_child_mesh(unsigned int naive_transformation_mode) {
    if (naive_transformation_mode == 1) {
        // infer the face center point from corner points
        for (int i = 0;i < _handles.cage_faces.size();i++) {
            int center_idx = 8 + i;
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (int j = 0;j < _handles.cage_faces[i].size();j++) {
                Eigen::Vector3d pos = _handles.handle_positions.row(_handles.cage_faces[i][j]) 
                                    + _handles.handle_transforms.block(4 * _handles.cage_faces[i][j] + 3, 0, 1, 3);
                centroid += pos;
            }
            centroid = centroid / _handles.cage_faces[i].size();
            _handles.move_handle(center_idx, centroid);
        }
    }
    MatrixX3d transformed_child_mesh = _lbs_mat_child_mesh * _handles.transform();

    return transformed_child_mesh;
}

void LinearBlendSkinningObject::compute_weights() {
    // Only compute weights if there are weights to compute!
    _handles.reset_handle_transformations();
    _lbs_mat.clear();
    if (_handles.positions().rows() > 0) {
        // compute weights for mesh
        for (int i = 0;i < _V.size();i++) {
            Eigen::MatrixXd W;
            ComputeMeanValueWeights(_V[i], _handles.positions(), _handles.cage_tri_mesh(), W);

            igl::normalize_row_sums(W, W);
            MatrixX lbs_mat;
            igl::lbs_matrix(_V[i], W, lbs_mat);
            _lbs_mat.push_back(lbs_mat);
        }

        // compute weights for tactile quad mesh
        if (_tactile_quad_mesh_V.rows() > 0) {
            Eigen::MatrixXd W;
            ComputeMeanValueWeights(_tactile_quad_mesh_V, _handles.positions(), _handles.cage_tri_mesh(), W);

            igl::normalize_row_sums(W, W);
            igl::lbs_matrix(_tactile_quad_mesh_V, W, _lbs_mat_tactile_mesh);
        }

        // compute weights for child mesh
        if (_child_mesh._V.cols() > 0) {
            Eigen::MatrixXd W;
            ComputeMeanValueWeights(_child_mesh._V.transpose(), _handles.positions(), _handles.cage_tri_mesh(), W);

            igl::normalize_row_sums(W, W);
            igl::lbs_matrix(_child_mesh._V.transpose(), W, _lbs_mat_child_mesh);
        }

        // compute weights for parent mesh
        if (_parent_mesh._V.cols() > 0) {
            Eigen::MatrixXd W;
            ComputeMeanValueWeights(_parent_mesh._V.transpose(), _handles.positions(), _handles.cage_tri_mesh(), W);

            igl::normalize_row_sums(W, W);
            igl::lbs_matrix(_parent_mesh._V.transpose(), W, _lbs_mat_parent_mesh);
        }

        _stale = true;
    } else {
        for (int i = 0;i < _V.size();i++) {
            _lbs_mat.push_back(MatrixXd(_V[i].cols(), 0));
        }
    }
}

void LinearBlendSkinningObject::move_handle(int handle_id, const Eigen::RowVector3d pos) {
    // skip operation if try to move parent interface handles.
    for (int i = 0;i < _interface_vertices_parent.size();i++) 
        if (_interface_vertices_parent[i].first == handle_id) {
            // std::cerr << "parent handle: " << i << std::endl;
            return;
        }
    
    // child interface handle
    bool is_interface_child_handle = false;
    for (int i = 0;i < _interface_vertices_children.size();i++)
        for (int j = 0;j < _interface_vertices_children[i].size();j++)
            if (_interface_vertices_children[i][j].first == handle_id) {
                is_interface_child_handle = true;
                // std::cerr << "child handle: " << i << " " << j << std::endl;

                LinearBlendSkinningObject* child = _children[i];
                if (_deformation_type == 0) { // FFD object
                    if (child->_deformation_type == 0) { // if child is also FFD object
                        // std::cerr << "branch FFD -> FFD" << std::endl;
                        _handles.move_handle(handle_id, pos);
                        child->FFD_from_parent();
                    } else {
                        // infer offset and scale by child
                        Eigen::RowVector3d old_pos = _handles.handle_positions.row(handle_id)
                                    + _handles._handle_translation.row(handle_id);
                        Eigen::RowVector3d translation = pos - old_pos;

                        // project onto child's scale plane
                        Eigen::RowVector3d projected_translation = Eigen::RowVector3d::Zero();
                        for (int k = 0;k < child->_scale_axis.size();k++)
                            projected_translation += translation.dot(child->_scale_axis[k]) * child->_scale_axis[k];

                        Eigen::RowVector3d offset = translation - projected_translation;
                        
                        Eigen::RowVector3d real_pos = old_pos + projected_translation;

                        // infer the scale
                        std::vector<double> scales;
                        for (int k = 0;k < child->_scale_axis.size();k++) {
                            if (fabs((old_pos - child->_scale_origin).dot(child->_scale_axis[k])) > 1e-7)
                                scales.push_back((real_pos - child->_scale_origin).dot(child->_scale_axis[k]) 
                                                / (old_pos - child->_scale_origin).dot(child->_scale_axis[k]));
                            else
                                scales.push_back(-10000.0);
                        }
                        
                        // apply offset and scale to all handles on child interface
                        for (int k = 0;k < _interface_vertices_children[i].size();k++) {
                            int affected_handle_id = _interface_vertices_children[i][k].first;
                            Eigen::RowVector3d last_handle_pos = _handles.handle_positions.row(affected_handle_id) 
                                                                + _handles._handle_translation.row(affected_handle_id);
                            Eigen::RowVector3d new_handle_pos = last_handle_pos;
                            for (int l = 0;l < child->_scale_axis.size();l++) {
                                if (scales[l] > -1000.0)
                                    new_handle_pos += ((scales[l] - 1.) * 
                                        (last_handle_pos - _scale_origin).dot(child->_scale_axis[l])) * child->_scale_axis[l];
                            }
                            new_handle_pos += offset;
                            _handles.move_handle(affected_handle_id, new_handle_pos);
                        }
                        child->offset_and_scale_from_parent(offset, scales, child->_scale_axis);
                    }
                } else {
                    move_single_handle(handle_id, pos);
                }

                return;
            }

    // internal handle
    if (!is_interface_child_handle) {
        // std::cerr << "internal handle" << std::endl;
        move_single_handle(handle_id, pos);
    }
}

void LinearBlendSkinningObject::finalize_operation() {
    _scale_origin_old = _scale_origin;
    _handles.finalize_operation();
}

void LinearBlendSkinningObject::visualize_handles(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& points, 
    Eigen::MatrixXd& point_colors, 
    Eigen::MatrixXi& lines, 
    Eigen::MatrixXd& line_colors,
    bool transform,
    unsigned int naive_deformation_mode)
{   
    int num_points = points.rows();
    int num_lines = lines.rows();
    
    Eigen::Vector3d point_handle_color(1.0, 0.7, 0.3); // Red
    Eigen::Vector3d bone_handle_color(0.0, 0.0, 1.0); // Green
    Eigen::Vector3d cage_handle_color(0.0, 1.0, 1.0); // Blue
    Eigen::Vector3d inactive_point_handle_color(0.5, 0.5, 0.5); // Grey

    if (naive_deformation_mode == 0) {
        points.conservativeResize(num_points + _handles.handle_positions.rows(), 3);
        lines.conservativeResize(num_lines + _handles.bone_edge_indices.rows() 
                                        + _handles.cage_edge_points.rows() 
                                        + _handles.cage_tri_mesh_elements.rows() * 3, 2);

        if (transform) {
            points.block(num_points, 0, _handles.handle_positions.rows(), 3) = 
                _handles.handle_lbs_matrix * _handles.handle_transforms;
        }
        else {
            points.block(num_points, 0, _handles.handle_positions.rows(), 3) = 
                _handles.handle_positions;
        }
        
        point_colors.conservativeResizeLike(points);
        lines.block(num_lines, 0, _handles.bone_edge_indices.rows(), 2) = _handles.bone_edge_indices;
        line_colors.conservativeResize(lines.rows(), 3);

        for (int i = 0; i < _handles.point_handle_indices.size(); ++i) {
            point_colors.row(num_points + _handles.point_handle_indices(i)) = point_handle_color;
        }

        for (int i = 0; i < _handles.bone_edge_indices.rows(); ++i) {
            point_colors.row(num_points + _handles.bone_edge_indices(i, 0)) = bone_handle_color;
            point_colors.row(num_points + _handles.bone_edge_indices(i, 1)) = bone_handle_color;
            line_colors.row(num_points + i) = bone_handle_color;
        }
        for (int i = 0; i < _handles.cage_edge_points.rows(); ++i) {
            int u = _handles.point_handle_indices(_handles.cage_edge_points(i, 0));
            int v = _handles.point_handle_indices(_handles.cage_edge_points(i, 1));
            point_colors.row(num_points + u) = cage_handle_color;
            point_colors.row(num_points + v) = cage_handle_color;
            int line_row = num_lines + _handles.bone_edge_indices.rows() + i;
            lines.row(line_row) = Eigen::Vector2i(num_points + u, num_points + v);
            line_colors.row(line_row) = cage_handle_color;
        }
        for (int i = 0;i < _handles.cage_tri_mesh_elements.rows();i++) {
            for (int j = 0;j < 3;j++) {
                int u = _handles.cage_tri_mesh_elements(i, j);
                int v = _handles.cage_tri_mesh_elements(i, (j + 1) % 3);
                point_colors.row(num_points + u) = cage_handle_color;
                point_colors.row(num_points + v) = cage_handle_color;
                int line_row = num_lines + _handles.bone_edge_indices.rows() + _handles.cage_edge_points.rows() + i * 3 + j;
                lines.row(line_row) = Eigen::Vector2i(num_points + u, num_points + v);
                line_colors.row(line_row) = cage_handle_color;
            }
        }

        for (int i = 0; i < _interface_vertices_parent.size(); ++i) {
            point_colors.row(num_points + _interface_vertices_parent[i].first) = inactive_point_handle_color;
        }
    } else { // only show the cuboid
        assert(_handles.handle_positions.rows() >= 8);
        points.conservativeResize(num_points + 8, 3);
        lines.conservativeResize(num_lines + 24, 2);

        if (transform) {
            points.block(num_points, 0, 8, 3) = 
                (_handles.handle_lbs_matrix * _handles.handle_transforms).topRows(8);
        }
        else {
            points.block(num_points, 0, 8, 3) = 
                _handles.handle_positions.topRows(8);
        }
        
        point_colors.conservativeResizeLike(points);
        
        line_colors.conservativeResize(lines.rows(), 3);

        for (int i = 0; i < 8; ++i) {
            point_colors.row(num_points + _handles.point_handle_indices(i)) = cage_handle_color;
        }

        int line_idx = num_lines;
        for (int i = 0;i < _handles.cage_faces.size();i++) {
            for (int j = 0;j < 4;j++) {
                int u = _handles.cage_faces[i][j];
                int v = _handles.cage_faces[i][(j + 1) % 4];
                lines.row(line_idx) = Eigen::Vector2i(num_points + u, num_points + v);
                line_colors.row(line_idx) = cage_handle_color;
                line_idx += 1;
            }
        }

        // for (int i = 0; i < _interface_vertices_parent.size(); ++i) {
        //     if (_interface_vertices_parent[i].first < 8)
        //         point_colors.row(num_points + _interface_vertices_parent[i].first) = inactive_point_handle_color;
        // }
    }
}

void LinearBlendSkinningObject::export_design_to_simulation(
        std::string name,
        std::string tab,
        std::string indent,
        std::string& robot_xml, 
        std::string& actuator_xml,
        std::map<LinearBlendSkinningObject*, TriMesh>& meshes,
        std::map<LinearBlendSkinningObject*, std::string>& mesh_paths,
        LinearBlendSkinningObject* current_joint_object,
        int &link_count,
        int depth) {
    
    Matrix3 coord_transform = (Matrix3() << 0., 0., -1.,
                                            -1., 0., 0.,
                                            0., 1., 0.).finished();

    std::vector<MatrixX3d> transformed_vertices = transform_mesh(1);
    
    // change the unit from mm into cm
    for (int i = 0;i < transformed_vertices.size();i++) {
        transformed_vertices[i] = transformed_vertices[i] * coord_transform.transpose(); // change coordinate
        transformed_vertices[i] /= 10.;
    }

    std::string link_name = "link" + to_string(link_count);
    std::string joint_name = "joint" + to_string(link_count);
    std::string body_name = "body" + to_string(link_count);
    std::string mesh_path = name + '/' + body_name + ".obj";

    bool new_link = true;
    // if (_part_type == Node::PartType::BODY && current_joint_object != nullptr && current_joint_object->_part_type == Node::PartType::JOINT) {
    //     new_link = false;
    // } else {
    //     new_link = true;
    //     robot_xml += indent + "<link name=\"" + link_name + "\">\n";
    //     link_count += 1;
    // }
    if (_part_type == Node::PartType::BODY) {
        // new_link = false;
        // attach to last joint
        // std::cerr << "current_joint_object = " << current_joint_object << std::endl;
        // meshes[current_joint_object].merge(transformed_vertices[0].transpose(), _F[0].transpose());
        // if (current_joint_object != nullptr && current_joint_object->_part_type == Node::PartType::JOINT) {
        if (current_joint_object != nullptr) {
            new_link = false;
            // attach to last joint
            meshes[current_joint_object].merge(transformed_vertices[0].transpose(), _F[0].transpose());
        } else {
            new_link = true;
            robot_xml += indent + "<link name=\"" + link_name + "\">\n";
            link_count += 1;
            Vector3d color(max(0., 0.7 - depth * 0.1), max(0., 0.7 - depth * 0.1), max(0., 0.7 - depth * 0.1));
            depth += 1;

            robot_xml += indent + tab + "<joint name=\"" + joint_name + "\"" + 
                        " type=\"fixed\"" +
                        " pos=\"0 0 0\""+
                        " quat=\"1 0 0 0\"" + 
                        " frame=\"WORLD\"/>\n";
            robot_xml += indent + tab + "<body name=\"" + body_name + "\"" + 
                            " type=\"mesh\" " +
                            " filename=\"" + mesh_path + "\" " +
                            " pos=\"0 0 0\" " +
                            " quat=\"1 0 0 0\" " +
                            " transform_type=\"OBJ_TO_WORLD\" " +
                            " density=\"1\" " +
                            " mu=\"0\"" + 
                            " rgba=\"" + to_string(color(0)) + " " + to_string(color(1)) + " " + to_string(color(2)) + " 1\"/>\n";
            TriMesh mesh;
            mesh._V = transformed_vertices[0].transpose();
            mesh._F = _F[0].transpose();
            meshes[this] = mesh;
            mesh_paths[this] = mesh_path;
            current_joint_object = this;
        }
    } else if (_part_type == Node::PartType::JOINT) {
        new_link = true;
        robot_xml += indent + "<link name=\"" + link_name + "\">\n";
        link_count += 1;

        std::pair<RowVector3d, RowVector3d> joint_axis = transform_joint_axis();
        // change the unit from mm into cm for joint axis
        joint_axis.first = joint_axis.first * coord_transform.transpose();
        joint_axis.first /= 10.;
        joint_axis.second = joint_axis.second * coord_transform.transpose();
        joint_axis.second /= 10.;
        
        RowVector3d axis = joint_axis.second - joint_axis.first;
        axis /= axis.norm();
        RowVector3d position = (joint_axis.first + joint_axis.second) / 2.;
        
        Vector3d color(max(0., 0.7 - depth * 0.1), max(0., 0.7 - depth * 0.1), max(0., 0.7 - depth * 0.1));
        depth += 1;

        robot_xml += indent + tab + "<joint name=\"" + joint_name + "\"" + 
                                " type=\"revolute\"" +
                                " axis=\"" + to_string(axis(0)) + " " + to_string(axis(1)) + " " + to_string(axis(2)) + "\"" + 
                                " pos=\"" + to_string(position(0)) + " " + to_string(position(1)) + " " + to_string(position(2)) + "\""+
                                " quat=\"1 0 0 0\"" + 
                                " frame=\"WORLD\"" +
                                " damping=\"1e4\"/>\n";

        robot_xml += indent + tab + "<body name=\"" + body_name + "\"" + 
                            " type=\"mesh\" " +
                            " filename=\"" + mesh_path + "\" " +
                            " pos=\"0 0 0\" " +
                            " quat=\"1 0 0 0\" " +
                            " transform_type=\"OBJ_TO_WORLD\" " +
                            " density=\"1\" " +
                            " mu=\"0\"" + 
                            " rgba=\"" + to_string(color(0)) + " " + to_string(color(1)) + " " + to_string(color(2)) + " 1\"/>\n";
        
        actuator_xml += tab + tab + "<motor joint=\"" + joint_name + "\" " + 
                            " ctrl=\"force\" " +
                            " ctrl_range=\"-3e5 3e5\"/>\n";

        if (_mesh_type == Node::MeshType::SINGLE) {
            TriMesh mesh;
            mesh._V = transformed_vertices[0].transpose();
            mesh._F = _F[0].transpose();
            meshes[this] = mesh;
            mesh_paths[this] = mesh_path;
        } else if (_mesh_type == Node::MeshType::SEPARATE) {
            // merge parent part to parent mesh
            if (current_joint_object != nullptr) {
                meshes[current_joint_object].merge(transformed_vertices[0].transpose(), _F[0].transpose());
            }
            // its own mesh
            TriMesh mesh;
            mesh._V = transformed_vertices[1].transpose();
            mesh._F = _F[1].transpose();
            for (int i = 2;i < _V.size();i++)
                mesh.merge(transformed_vertices[i].transpose(), _F[i].transpose());
            meshes[this] = mesh;
            mesh_paths[this] = mesh_path;
        }
        current_joint_object = this;
    }
    
    for (int i = 0;i < _children.size();i++)
        _children[i]->export_design_to_simulation(
            name, tab, indent + tab,
            robot_xml, actuator_xml, 
            meshes, mesh_paths, current_joint_object, 
            link_count,
            depth);
    
    if (new_link)
        robot_xml += indent + "</link>\n";
}

TriMesh LinearBlendSkinningObject::export_printing_meshes(
    std::vector<TriMesh>& meshes) {
    
    TriMesh mesh;
    for (int i = 0;i < _children.size();i++) {
        TriMesh child_mesh = _children[i]->export_printing_meshes(meshes);
        mesh.merge(child_mesh);
    }

    std::vector<MatrixX3d> transformed_vertices = transform_mesh(1);
    if (_part_type == Node::PartType::BODY) {
        mesh.merge(transformed_vertices[0].transpose(), _F[0].transpose());
    } else {
        mesh.merge(transformed_vertices[1].transpose(), _F[1].transpose());
        meshes.push_back(mesh);
        mesh = TriMesh(transformed_vertices[0].transpose(), _F[0].transpose());
    }
    // add additional meshes
    for (int i = 2;i < transformed_vertices.size();i++)
        meshes.push_back(TriMesh(transformed_vertices[i].transpose(), _F[i].transpose()));

    // merge parent's mesh of palm
    if (_parent->_symbol == "palm") {
        MatrixX3d transformed_vertices_parent = _parent->transform_child_mesh(1);
        mesh.merge(transformed_vertices_parent.transpose(), _parent->_child_mesh._F);
        meshes.push_back(mesh);
    }

    return mesh;
}

void LinearBlendSkinningObject::export_tactile_mesh(
        std::vector<QuadMesh>& meshes,
        std::vector<Vector3> tactile_mesh_vertices, 
        std::vector<Vector4i> tactile_mesh_faces) {
    Eigen::MatrixX3d transformed_mesh = transform_tactile_mesh(1);
    if (_symbol != "fork") {
        if (_parent->_symbol == "palm" || _parent->_symbol == "P" || _parent->_symbol == "fork") {
            for (int i = 0;i < 4;i++) {
                _tactile_quad_verts_index.push_back(tactile_mesh_vertices.size());
                tactile_mesh_vertices.push_back(transformed_mesh.row(i).transpose());
            }
        } else {
            for (int i = 4;i < 8;i++)
                _tactile_quad_verts_index.push_back(_parent->_tactile_quad_verts_index[i]);
        }
        for (int i = 4;i < 8;i++) {
            _tactile_quad_verts_index.push_back(tactile_mesh_vertices.size());
            tactile_mesh_vertices.push_back(transformed_mesh.row(i).transpose());
        }
        int faces[4][4] = {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}};
        for (int i = 0;i < 4;i++) {
            Vector4i face;
            for (int j = 0;j < 4;j++)
                face[j] = _tactile_quad_verts_index[faces[i][j]];
            tactile_mesh_faces.push_back(face);
        }
    }

    if ((_symbol == "fork" || _children.size() == 0) && tactile_mesh_vertices.size() > 0 && tactile_mesh_faces.size() > 0) {
        Matrix3Xd V(3, tactile_mesh_vertices.size());
        Matrix4Xi F(4, tactile_mesh_faces.size());
        for (int i = 0;i < tactile_mesh_vertices.size();i++)
            V.col(i) = tactile_mesh_vertices[i];
        for (int j = 0;j < tactile_mesh_faces.size();j++)
            F.col(j) = tactile_mesh_faces[j];
        QuadMesh mesh(V, F);
        meshes.push_back(mesh);

        tactile_mesh_vertices.clear();
        tactile_mesh_faces.clear();
    } 
    
    for (int i = 0;i < _children.size();i++)
        _children[i]->export_tactile_mesh(meshes, tactile_mesh_vertices, tactile_mesh_faces);
}

void LinearBlendSkinningObject::move_single_handle(
        int handle_id, 
        const Eigen::RowVector3d pos) {

    if (_deformation_type == 0) {   // FFD object
        _handles.move_handle(handle_id, pos);
    } else {                        // Scale object, must be a component only have one child (e.g. joint)
        // project onto scale plane
        Eigen::RowVector3d old_pos = _handles.handle_positions.row(handle_id)
                                    + _handles._handle_translation.row(handle_id);
        Eigen::RowVector3d translation = pos - old_pos;

        Eigen::RowVector3d projected_translation = Eigen::RowVector3d::Zero();
        for (int i = 0;i < _scale_axis.size();i++)
            projected_translation += translation.dot(_scale_axis[i]) * _scale_axis[i];

        Eigen::RowVector3d real_pos = old_pos + projected_translation;

        // infer the scale
        std::vector<double> scales;
        for (int i = 0;i < _scale_axis.size();i++) {
            if (fabs((old_pos - _scale_origin).dot(_scale_axis[i])) > 1e-7)
                scales.push_back((real_pos - _scale_origin).dot(_scale_axis[i]) 
                                / (old_pos - _scale_origin).dot(_scale_axis[i]));
            else
                scales.push_back(-10000.0);
        }
        
        // apply scale to all handles
        for (int i = 0;i < _handles.handle_positions.rows();i++) {
            Eigen::RowVector3d last_handle_pos = _handles.handle_positions.row(i) 
                                                + _handles._handle_translation.row(i);
            Eigen::RowVector3d new_handle_pos = last_handle_pos;
            for (int j = 0;j < _scale_axis.size();j++) {
                if (scales[j] > -1000.0)
                    new_handle_pos += ((scales[j] - 1.) * 
                        (last_handle_pos - _scale_origin).dot(_scale_axis[j])) * _scale_axis[j];
            }
            _handles.move_handle(i, new_handle_pos);
        }

        if (_parent != nullptr) {   // upward propagation
            if (_parent->_deformation_type == 0) {
                _parent->FFD_from_child(this);
            } else { // parent must be a component only have one child (e.g. joint)
                _parent->scale_from_child(scales, _scale_axis);
            }
        }

        if (_children.size() > 0 && _children[0] != nullptr) {
            if (_children[0]->_deformation_type == 0) {
                _children[0]->FFD_from_parent();
            } else {
                _children[0]->scale_from_parent(scales, _scale_axis);
            }
        }
    }
}

void LinearBlendSkinningObject::offset_and_scale_from_parent(
        Eigen::RowVector3d offset, 
        std::vector<double> scales, 
        std::vector<Eigen::RowVector3d> scale_axis) {

    // check if valid for scale
    if (_deformation_type == 1) {
        for (int i = 0;i < scale_axis.size();i++) {
            bool flag = false;
            for (int j = 0;j < _scale_axis.size();j++)
                if ((scale_axis[i] - _scale_axis[j]).norm() < 1e-7) {
                    flag = true;
                    break;
                }
            if (!flag) {
                scales.clear();
                scale_axis.clear();
            }
        }
    } else {
        scales.clear();
        scale_axis.clear();
    }

    // scale and offset
    for (int i = 0;i < _handles.handle_positions.rows();i++) {
        Eigen::RowVector3d last_handle_pos = _handles.handle_positions.row(i) 
                                            + _handles._handle_translation.row(i);
        Eigen::RowVector3d new_handle_pos = last_handle_pos;
        for (int j = 0;j < scale_axis.size();j++) {
            if (scales[j] > -1000.0)
                new_handle_pos += ((scales[j] - 1.) * 
                    (last_handle_pos - _scale_origin).dot(scale_axis[j])) * scale_axis[j];
        }
        _handles.move_handle(i, new_handle_pos + offset);
    }

    _scale_origin += offset;

    // if is a FFD object, then align the interface points
    FFD_from_parent();

    // must be a component only have one child (e.g. joint)
    for (int i = 0;i < _children.size();i++) {
        _children[i]->offset_and_scale_from_parent(offset, scales, _scale_axis);
    }
}

void LinearBlendSkinningObject::scale_from_child(
        std::vector<double> scales, 
        std::vector<Eigen::RowVector3d> scale_axis) {
    
    // check if valid
    for (int i = 0;i < scale_axis.size();i++) {
        bool flag = false;
        for (int j = 0;j < _scale_axis.size();j++)
            if ((scale_axis[i] - _scale_axis[j]).norm() < 1e-7) {
                flag = true;
                break;
            }
        assert(flag);
        if (!flag)
            return;
    }

    // scale
    for (int i = 0;i < _handles.handle_positions.rows();i++) {
        Eigen::RowVector3d last_handle_pos = _handles.handle_positions.row(i) 
                                            + _handles._handle_translation.row(i);
        Eigen::RowVector3d new_handle_pos = last_handle_pos;
        for (int j = 0;j < scale_axis.size();j++) {
            if (scales[j] > -1000.0)
                new_handle_pos += ((scales[j] - 1.) * 
                    (last_handle_pos - _scale_origin).dot(scale_axis[j])) * scale_axis[j];
        }
        _handles.move_handle(i, new_handle_pos);
    }

    if (_parent != nullptr) {   // upward propagation
        if (_parent->_deformation_type == 0) {
            _parent->FFD_from_child(this);
        } else { // parent must be a component only have one child (e.g. joint)
            _parent->scale_from_child(scales, _scale_axis);
        }
    }
}

void LinearBlendSkinningObject::scale_from_parent(
        std::vector<double> scales, 
        std::vector<Eigen::RowVector3d> scale_axis) {
    
    // check if valid
    for (int i = 0;i < scale_axis.size();i++) {
        bool flag = false;
        for (int j = 0;j < _scale_axis.size();j++)
            if ((scale_axis[i] - _scale_axis[j]).norm() < 1e-7) {
                flag = true;
                break;
            }
        assert(flag);
        if (!flag)
            return;
    }

    // scale
    for (int i = 0;i < _handles.handle_positions.rows();i++) {
        Eigen::RowVector3d last_handle_pos = _handles.handle_positions.row(i) 
                                            + _handles._handle_translation.row(i);
        Eigen::RowVector3d new_handle_pos = last_handle_pos;
        for (int j = 0;j < scale_axis.size();j++) {
            if (scales[j] > -1000.0)
                new_handle_pos += ((scales[j] - 1.) * 
                    (last_handle_pos - _scale_origin).dot(scale_axis[j])) * scale_axis[j];
        }
        _handles.move_handle(i, new_handle_pos);
    }

    // must be a component only have one child (e.g. joint)
    if (_children[0] != nullptr) {
        if (_children[0]->_deformation_type == 0) {
            _children[0]->FFD_from_parent();
        } else {
            _children[0]->scale_from_parent(scales, _scale_axis);
        }
    }
}

void LinearBlendSkinningObject::FFD_from_child(LinearBlendSkinningObject* child) {
    for (int i = 0;i < _children.size();i++)
        if (_children[i] == child) {
            for (int j = 0;j < _interface_vertices_children[i].size();j++) {
                Eigen::RowVector3d new_handle_pos = 
                    child->_handles.handle_positions.row(_interface_vertices_children[i][j].second)
                    + child->_handles.handle_transforms.block(4 * _interface_vertices_children[i][j].second + 3, 0, 1, 3);
                _handles.move_handle(_interface_vertices_children[i][j].first, new_handle_pos);                    
            }
        }
}

void LinearBlendSkinningObject::FFD_from_parent() {
    for (int j = 0;j < _interface_vertices_parent.size();j++) {
        Eigen::RowVector3d new_handle_pos = 
            _parent->_handles.handle_positions.row(_interface_vertices_parent[j].second)
            + _parent->_handles.handle_transforms.block(4 * _interface_vertices_parent[j].second + 3, 0, 1, 3);
        _handles.move_handle(_interface_vertices_parent[j].first, new_handle_pos);                    
    }
}

void LinearBlendSkinningObject::reset_naive_mode_parameters() {
    if (_symbol != "palm" && _symbol != "P") {
        _length = 0.;
        for (int i = 0;i < 8;i++) {
            double project = (_handles.handle_positions.row(i).transpose() - _parent_face_origin).dot(_direction);
            if (project < 1e-7) {
                _parent_face_vertices.push_back(i);
            } else {
                _child_face_vertices.push_back(i);
            }
            _length = max(_length, (float)project);
        }
        _original_length = _length;
        for (int i = 0;i < 4;i++)
            for (int j = i;j < 4;j++) {
                double project_x_i = (_handles.handle_positions.row(_parent_face_vertices[i]).transpose() - _parent_face_origin).dot(_direction_x);
                double project_x_j = (_handles.handle_positions.row(_parent_face_vertices[j]).transpose() - _parent_face_origin).dot(_direction_x);
                double project_y_i = (_handles.handle_positions.row(_parent_face_vertices[i]).transpose() - _parent_face_origin).dot(_direction_y);
                double project_y_j = (_handles.handle_positions.row(_parent_face_vertices[j]).transpose() - _parent_face_origin).dot(_direction_y);
                if (project_x_j < project_x_i - 1e-7 || (fabs(project_x_j - project_x_i) < 1e-7 && project_y_j < project_y_i)) {
                    swap(_parent_face_vertices[i], _parent_face_vertices[j]);
                }
            }
        for (int i = 0;i < 4;i++)
            for (int j = i;j < 4;j++) {
                double project_x_i = (_handles.handle_positions.row(_child_face_vertices[i]).transpose() - _parent_face_origin).dot(_direction_x);
                double project_x_j = (_handles.handle_positions.row(_child_face_vertices[j]).transpose() - _parent_face_origin).dot(_direction_x);
                double project_y_i = (_handles.handle_positions.row(_child_face_vertices[i]).transpose() - _parent_face_origin).dot(_direction_y);
                double project_y_j = (_handles.handle_positions.row(_child_face_vertices[j]).transpose() - _parent_face_origin).dot(_direction_y);
                if (project_x_j < project_x_i - 1e-7 || (fabs(project_x_j - project_x_i) < 1e-7 && project_y_j < project_y_i)) {
                    swap(_child_face_vertices[i], _child_face_vertices[j]);
                }
            }
        _child_face_width = (_handles.handle_positions.row(_child_face_vertices[3]).transpose() - _parent_face_origin).dot(_direction_x)
                            - (_handles.handle_positions.row(_child_face_vertices[0]).transpose() - _parent_face_origin).dot(_direction_x);
        _original_child_face_height = _child_face_height 
                            = (_handles.handle_positions.row(_child_face_vertices[1]).transpose() - _parent_face_origin).dot(_direction_y)
                            - (_handles.handle_positions.row(_child_face_vertices[0]).transpose() - _parent_face_origin).dot(_direction_y);
        _shear_x = 0.;
        _shear_y = 0.;

        // if parent is fork, source back to the first non-fork ancester and compute the translation to the parent connection face of first fork component.
        if (_parent != nullptr && _parent->_symbol == "fork") {
            auto now = this;
            while (now->_parent != nullptr && now->_parent->_symbol == "fork") {
                now = now->_parent;
            }
            _fork_translation.clear();
            assert(_parent_face_vertices.size() == now->_parent_face_vertices.size());
            for (int i = 0;i < _parent_face_vertices.size();i++)
                _fork_translation.push_back(_handles.handle_positions.row(_parent_face_vertices[i]).transpose() - now->_handles.handle_positions.row(now->_parent_face_vertices[i]).transpose());
        }
    } else {
        _handle_vertices.clear();
        for (int i = 0;i < 8;i++)
            _handle_vertices.push_back(i);
        for (int i = 0;i < 8;i++)
            for (int j = i;j < 8;j++) {
                double project_x_i = _handles.handle_positions.row(_handle_vertices[i])[0];
                double project_x_j = _handles.handle_positions.row(_handle_vertices[j])[0];
                double project_y_i = _handles.handle_positions.row(_handle_vertices[i])[1];
                double project_y_j = _handles.handle_positions.row(_handle_vertices[j])[1];
                double project_z_i = _handles.handle_positions.row(_handle_vertices[i])[2];
                double project_z_j = _handles.handle_positions.row(_handle_vertices[j])[2];
                if ((project_x_j < project_x_i - 1e-7) ||
                    (fabs(project_x_j - project_x_i) < 1e-7 && project_y_j < project_y_i - 1e-7) ||
                    (fabs(project_x_j - project_x_i) < 1e-7 && fabs(project_y_j - project_y_i) < 1e-7 && project_z_j < project_z_i - 1e-7)) {
                    swap(_handle_vertices[i], _handle_vertices[j]);
                }
            } 
        _cell_width = _handles.handle_positions.row(_handle_vertices[2])[1]
                        - _handles.handle_positions.row(_handle_vertices[0])[1];
        _cell_height = _handles.handle_positions.row(_handle_vertices[1])[0]
                        - _handles.handle_positions.row(_handle_vertices[0])[0];
    }
}

void LinearBlendSkinningObject::update_naive_deformation(PalmGrammar& palm_grammar, std::vector<LinearBlendSkinningObject*>& objects) {
    if (_symbol == "palm") {
        // update width
        int center_row = 0, center_col = palm_grammar._num_cols / 2;
        double y = palm_grammar._cell_widths[center_col] / 2.;
        for (int i = 0;i < palm_grammar._num_rows;i++)
            if (palm_grammar._cell_index(i, center_col) != -1) {
                int idx = palm_grammar._cell_index(i, center_col);
                int handle_ids_left[4] = {0, 1, 4, 5};
                int handle_ids_right[4] = {2, 3, 6, 7};
                for (int k = 0;k < 4;k++) {
                    RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_left[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_left[k]] + 3, 0, 1, 3);
                    new_pos[1] = -y;
                    objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_left[k]], new_pos);
                }
                for (int k = 0;k < 4;k++) {
                    RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_right[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_right[k]] + 3, 0, 1, 3);
                    new_pos[1] = y;
                    objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_right[k]], new_pos);
                }
            }
        for (int j = center_col + 1;j < palm_grammar._num_cols;j++)  {
            double next_y = y + palm_grammar._cell_widths[j];
            for (int i = 0;i < palm_grammar._num_rows;i++)
                if (palm_grammar._cell_index(i, j) != -1) {
                    int idx = palm_grammar._cell_index(i, j);
                    int handle_ids_left[4] = {0, 1, 4, 5};
                    int handle_ids_right[4] = {2, 3, 6, 7};
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_left[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_left[k]] + 3, 0, 1, 3);
                        new_pos[1] = y;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_left[k]], new_pos);
                    }
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_right[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_right[k]] + 3, 0, 1, 3);
                        new_pos[1] = next_y;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_right[k]], new_pos);
                    }
                }
            y = next_y;
        }
        y = -palm_grammar._cell_widths[center_col] / 2.;
        for (int j = center_col - 1;j >= 0;j--)  {
            double next_y = y - palm_grammar._cell_widths[j];
            for (int i = 0;i < palm_grammar._num_rows;i++)
                if (palm_grammar._cell_index(i, j) != -1) {
                    int idx = palm_grammar._cell_index(i, j);
                    int handle_ids_left[4] = {0, 1, 4, 5};
                    int handle_ids_right[4] = {2, 3, 6, 7};
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_left[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_left[k]] + 3, 0, 1, 3);
                        new_pos[1] = next_y;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_left[k]], new_pos);
                    }
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_ids_right[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_ids_right[k]] + 3, 0, 1, 3);
                        new_pos[1] = y;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_ids_right[k]], new_pos);
                    }
                }
            y = next_y;
        }

        // update height
        double x = -palm_grammar._cell_heights[0] / 2.;
        for (int i = 0;i < palm_grammar._num_rows;i++) {
            double next_x = x + palm_grammar._cell_heights[i];
            for (int j = 0;j < palm_grammar._num_cols;j++)
                if (palm_grammar._cell_index(i, j) != -1) {
                    int idx = palm_grammar._cell_index(i, j);
                    int handle_idx_bottom[4] = {0, 1, 2, 3};
                    int handle_idx_top[4] = {4, 5, 6, 7};
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_idx_bottom[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_idx_bottom[k]] + 3, 0, 1, 3);
                        new_pos[0] = x;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_idx_bottom[k]], new_pos);
                    }
                    for (int k = 0;k < 4;k++) {
                        RowVector3d new_pos = objects[idx]->_handles.handle_positions.row(objects[idx]->_handle_vertices[handle_idx_top[k]])
                                    + objects[idx]->_handles.handle_transforms.block(4 * objects[idx]->_handle_vertices[handle_idx_top[k]] + 3, 0, 1, 3);
                        new_pos[0] = next_x;
                        objects[idx]->_handles.move_handle(objects[idx]->_handle_vertices[handle_idx_top[k]], new_pos);
                    }
                }
            x = next_x;
        }

        // update finger
        for (int i = 0;i < palm_grammar._num_rows;i++) {
            for (int j = 0;j < palm_grammar._num_cols;j++) {
                if (palm_grammar._connector_symbols.find(palm_grammar._palm_grid(i, j)) != palm_grammar._connector_symbols.end()) {
                    // compute the translation of the connection face
                    Vector3 translation = objects[palm_grammar._cell_index(i, j)]->_handles.get_handle_position(1) - objects[palm_grammar._cell_index(i, j)]->_handles.handle_positions.row(1).transpose();
                    // propagate the translation to child finger component
                    auto child = objects[palm_grammar._cell_index(i, j)]->_children[0];
                    if (child != nullptr) {
                        for (int k = 0;k < child->_parent_face_vertices.size();k++) {
                            int child_idx = child->_parent_face_vertices[k];
                            Vector3 pos = child->_handles.handle_positions.row(child_idx).transpose() + translation;
                            child->_handles.move_handle(child_idx, pos);
                            child->propagate_downwards_naive_deformation();
                        }
                    }
                }
            }
        }

    } else {
        if (_symbol == "joint") {
            _child_face_height = _length * _original_child_face_height / _original_length;
        }
        auto now = this;
        while (now->_parent != nullptr && now->_symbol == "joint") {
            if (now->_parent->_symbol == "palm") {
                break;
            }
            now = now->_parent;
        }
        now->_child_face_height = _child_face_height;
        now->_child_face_width = _child_face_width;
        now->propagate_downwards_naive_deformation();
    }
}

void LinearBlendSkinningObject::propagate_downwards_naive_deformation() {
    Vector3 parent_face_origin = Vector3::Zero();

    // ------------ special case for fork-connected child component ------------
    if (_parent != nullptr && _parent->_symbol == "fork") {
        auto now = this;
        while (now->_parent != nullptr && now->_parent->_symbol == "fork") {
            now = now->_parent;
        }
        for (int i = 0;i < _parent_face_vertices.size();i++) {
            Vector3 new_pos = now->_handles.get_handle_position(now->_parent_face_vertices[i]) + _fork_translation[i];
            _handles.move_handle(_parent_face_vertices[i], new_pos);
        }
    }
    // -------------------------------------------------------------------------

    for (int i = 0;i < _parent_face_vertices.size();i++) {
        int idx = _parent_face_vertices[i];
        parent_face_origin += _handles.handle_positions.row(idx) + _handles.handle_transforms.block(4 * idx + 3, 0, 1, 3);
    }
    parent_face_origin = parent_face_origin / _parent_face_vertices.size();

    float sign_x[4] = {-1, -1, 1, 1};
    float sign_y[4] = {-1, 1, -1, 1};

    if (_symbol == "palm") {
    } else {
        if (_symbol == "joint") {
            _child_face_height = _length * _original_child_face_height / _original_length;
        }
        for (int i = 0;i < 4;i++) {
            int idx = _child_face_vertices[i];
            Vector3 new_position = parent_face_origin;
            new_position += _length * _direction + sign_x[i] * _child_face_width / 2. * _direction_x + sign_y[i] * _child_face_height / 2. * _direction_y;
            new_position += _direction_x * _shear_x + _direction_y * _shear_y;
            _handles.move_handle(idx, new_position);
        }
        if (_children.size() > 0) {
            for (int child_id = 0;child_id < _children.size();child_id++) {
                for (int i = 0;i < 4;i++) {
                    Vector3 pos = _handles.handle_positions.row(_child_face_vertices[i]) + _handles.handle_transforms.block(4 * _child_face_vertices[i] + 3, 0, 1, 3);
                    _children[child_id]->_handles.move_handle(_children[child_id]->_parent_face_vertices[i], pos);
                }
                if (_symbol != "fork" && _children[child_id]->_symbol == "joint") {
                    _children[child_id]->_child_face_width = _child_face_width;
                    _children[child_id]->_length = _children[child_id]->_original_length / _children[child_id]->_original_child_face_height * _child_face_height;
                }
                _children[child_id]->propagate_downwards_naive_deformation();
            }
        }
    }

}

void LinearBlendSkinningObject::scale_design(dtype scale_ratio) {
    for (int i = 0;i < _V.size();i++)
        _V[i] = _V[i] * scale_ratio;
    _tactile_quad_mesh_V = _tactile_quad_mesh_V * scale_ratio;
    _parent_mesh._V = _parent_mesh._V * scale_ratio;
    _child_mesh._V = _child_mesh._V * scale_ratio;

    // scale lbs matrix
    int m = _lbs_mat[0].cols();
    for (int j = 0;j < m;j++) 
        if (j % 4 < 3) {
            for (int i = 0;i < _lbs_mat.size();i++) {
                _lbs_mat[i].col(j) *= scale_ratio;
            }
            _lbs_mat_tactile_mesh.col(j) *= scale_ratio;
            _lbs_mat_parent_mesh.col(j) *= scale_ratio;
            _lbs_mat_child_mesh.col(j) *= scale_ratio;
            _handles.handle_lbs_matrix.col(j) *= scale_ratio;
        }

    _handles.handle_positions *= scale_ratio;
    _handles._handle_translation *= scale_ratio;
    for (int i = 0;i < _handles.handle_positions.rows();i++)
        _handles.handle_transforms.block(4 * i + 3, 0, 1, 3) *= scale_ratio;

    for (int i = 0;i < _transformed_V.size();i++)
        _transformed_V[i] *= scale_ratio;
    
    _scale_origin *= scale_ratio;
    _scale_origin_old *= scale_ratio;

    _parent_face_origin *= scale_ratio;
    
    _length *= scale_ratio;
    _shear_x *= scale_ratio;
    _shear_y *= scale_ratio;
    _child_face_width *= scale_ratio;
    _child_face_height *= scale_ratio;
    _original_length *= scale_ratio;
    _original_child_face_height *= scale_ratio;

    _cell_width *= scale_ratio;
    _cell_height *= scale_ratio;
    
    for (int i = 0;i < _fork_translation.size();i++)
        _fork_translation[i] *= scale_ratio;
}