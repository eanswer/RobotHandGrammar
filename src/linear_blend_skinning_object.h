#pragma once

#include "Common/Common.h"
#include "handles.hpp"
#include "Graph/node.h"

class PalmGrammar;

class LinearBlendSkinningObject {
public:
    LinearBlendSkinningObject();
    ~LinearBlendSkinningObject();

    void set_mesh_type(Node::MeshType mesh_type);

    void set_mesh(std::vector<Matrix3Xd> V, 
                    std::vector<Matrix3Xi> F);
    
    void set_mesh(std::vector<MatrixX3d> V, 
                    std::vector<MatrixX3i> F);

    void set_part_type(Node::PartType part_type);

    void set_joint_axis(std::pair<Vector3d, Vector3d> joint_axis);

    void set_weight_type(int weight_type);

    std::vector<MatrixX3d> transform_mesh(unsigned int naive_transformation_mode);
    std::pair<RowVectorXd, RowVectorXd> transform_joint_axis();
    MatrixX3d transform_tactile_mesh(unsigned int naive_transformation_mode);
    MatrixX3d transform_parent_mesh(unsigned int naive_transformation_mode);
    MatrixX3d transform_child_mesh(unsigned int naive_transformation_mode);

    void compute_weights();

    void move_handle(int handle_id, const Eigen::RowVector3d pos);

    void finalize_operation();

    void visualize_handles(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
		Eigen::MatrixXd& points, 
		Eigen::MatrixXd& point_colors, 
		Eigen::MatrixXi& lines, 
		Eigen::MatrixXd& line_colors,
		bool transform = false,
        unsigned int naive_deformation_mode = 1);
    
    void export_design_to_simulation(
        std::string name,
        std::string tab,
        std::string indent,
        std::string& robot_xml, 
        std::string& actuator_xml,
        std::map<LinearBlendSkinningObject*, TriMesh>& meshes,
        std::map<LinearBlendSkinningObject*, std::string>& mesh_paths,
        LinearBlendSkinningObject* current_joint_object,
        int &link_count,
        int depth);

    TriMesh export_printing_meshes(std::vector<TriMesh>& meshes);

    void export_tactile_mesh(
        std::vector<QuadMesh>& meshes,
        std::vector<Vector3> tactile_mesh_vertices, 
        std::vector<Vector4i> tactile_mesh_faces);

    void reset_naive_mode_parameters();

    void update_naive_deformation(PalmGrammar& palm_grammar, std::vector<LinearBlendSkinningObject*>& objects);
    void propagate_downwards_naive_deformation();

    void scale_design(dtype scale_ratio);

    Node::MeshType _mesh_type;

    std::string _symbol;

    std::vector<MatrixX3d> _V;    // Vertex Positions in rest configuration
    std::vector<MatrixX3i> _F;   // Triangular Faces of exterior surface
    MatrixX3d _tactile_quad_mesh_V;
    std::vector<int> _tactile_quad_verts_index;
    TriMesh _parent_mesh;
    TriMesh _child_mesh;

    Node::PartType _part_type;
    std::pair<Vector3d, Vector3d> _joint_axis;

    int _weight_type; // Type of weights to compute
	                     // 0 = Nearest Neighbor
	                     // 1 = Linear
	                     // 2 = Bounded Biharmonic
                         // 3 = Mean Value Weights
    
    unsigned _deformation_type; 

    std::vector<MatrixXd> _lbs_mat;
    std::pair<RowVectorXd, RowVectorXd> _lbs_mat_joint_axis;
    MatrixXd _lbs_mat_tactile_mesh;
    MatrixXd _lbs_mat_parent_mesh;
    MatrixXd _lbs_mat_child_mesh;
    WeightHandles _handles;

    std::vector<MatrixX3d> _transformed_V;
    bool _stale;

    LinearBlendSkinningObject* _parent;
    std::vector<LinearBlendSkinningObject*> _children;
    // the vertices on the parent interface, a pair of index in itself and parent
    std::vector<std::pair<int, int> > _interface_vertices_parent; 
    // the vertices on the parent interface, a pair of index in itself and parent
    std::vector<std::vector<std::pair<int, int> > > _interface_vertices_children;

    // original points used for scale operation
    Eigen::RowVector3d _scale_origin;
    Eigen::RowVector3d _scale_origin_old;

    Eigen::Vector3d _parent_face_origin;
    Eigen::Vector3d _direction;
    Eigen::Vector3d _direction_x;
    Eigen::Vector3d _direction_y;


    // scale axis
    std::vector<Eigen::RowVector3d> _scale_axis;

    // parameters for naive deformation mode
    float _length, _shear_x, _shear_y, _child_face_width, _child_face_height; // for finger
    float _original_length, _original_child_face_height;

    float _cell_width, _cell_height; // for palm
    int _cell_row, _cell_col;

    std::vector<int> _parent_face_vertices; // for finger
    std::vector<int> _child_face_vertices;
    std::vector<int> _handle_vertices; // for palm 

    // fork translation
    std::vector<Vector3> _fork_translation;

private:
    void move_single_handle(int handle_id, const Eigen::RowVector3d pos);
    void offset_and_scale_from_parent(Eigen::RowVector3d offset, std::vector<double> scales, std::vector<Eigen::RowVector3d> scale_axis);
    void scale_from_child(std::vector<double> scales, std::vector<Eigen::RowVector3d> scale_axis);
    void scale_from_parent(std::vector<double> scales, std::vector<Eigen::RowVector3d> scale_axis);
    void FFD_from_child(LinearBlendSkinningObject* child);
    void FFD_from_parent();
};