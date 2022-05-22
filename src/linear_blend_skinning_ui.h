#pragma once

#define IGL_VIEWER_VIEWER_QUIET 1
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/lbs_matrix.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject_in_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/decimate.h>
#include <GL/freeglut.h>

#include "handles.hpp"
#include "linear_blend_skinning_object.h"
#include "Grammar/finger_grammar.h"
#include "Grammar/palm_grammar.h"
#include "Mesh/mesh.h"

class CameraConfig {
public:
    float _camera_base_translation[3];
    float _camera_base_zoom;
    float _camera_translation[3];
    float _camera_zoom;
    float _trackball_angle[4];

    void get(igl::opengl::glfw::Viewer* viewer);
    void set(igl::opengl::glfw::Viewer* viewer);

    void save_to_file(std::string path);
    void load_from_file(std::string path);
};

class LinearBlendSkinningUI : public igl::opengl::glfw::ViewerPlugin {
public:
	LinearBlendSkinningUI();

	void init(igl::opengl::glfw::Viewer* _viewer);

    void reset_palm();
    void reset_finger();

    // void reset_design_graph();

    void construct_meshes();

    void construct_finger_objects(Tree &graph);

    void update_object_mesh(std::vector<MatrixX3d>& V, 
                                std::vector<MatrixX4i>& T, 
                                std::vector<MatrixX3i>& F);

    void update_viewer_mesh();

    void update_colors();
    
    void update_viewer(bool update_mesh = true);

    void start_finger_grammar_stage();
    void start_deformation_stage();

    void scale_design();

    // deform the mesh of each object
    MatrixX3d transform_mesh();

    void draw_mesh();

	void draw_handles();

	void draw_menu();

    int get_closest_mesh_vertex();

    int get_closest_handle_id(int object_id);

    int query_object_by_vertex(int vertex_id);

	bool mouse_down(int button, int modifier);
	bool mouse_up(int button, int modifier);
	bool mouse_move(int mouse_x, int mouse_y);
    bool key_down(int key, int modifiers);
    bool key_up(int key, int modifiers);

    void export_design(std::string name = "design");

	static const int MAX_FACES = 10000;

    enum STAGE {
        PALM = 0, // 0
        FINGER, // 1
        DEFORMATION // 2
    };

	igl::opengl::glfw::imgui::ImGuiMenu menu;

    PalmGrammar _palm_grammar;

    FingerGrammar _finger_grammar;
    
    // meshes and colors
    std::vector<TriMesh> _meshes;
    std::vector<Vector3> _colors; // color for each mesh
    
    // merged meshed for visualization
	MatrixX3d _V; // Vertex Positions Union in rest configuration (objects)
	MatrixX4i _T; // Tetrahedral Elements Union (objects)
	MatrixX3i _F; // Triangular Faces of exterior surface Union (objects)
    MatrixX3d _vert_colors; // colors of each mesh vertex
    
    // hovering id
    int _hover_rule_id;

	// Handle Manipulation State
	int _moving_handle_id = -1;
	Eigen::RowVector3d _sel_pos; // Position of handle when it was selected

    // list of deformable objects
    std::vector<LinearBlendSkinningObject*> _objects;
    std::vector<Node::NodeType> _object_types;

    // handle constraints is a map from <object_id, handle_id> to <object_id, handle_id>
    // the movement of two paired handles should be identical
    std::map<std::pair<int, int>, std::pair<int, int> > _handle_constraints; 

    std::string _design_name;
    
    int _weight_type; // Type of weights to compute
                        // 0 = Nearest Neighbor
                        // 1 = Linear
                        // 2 = Bounded Biharmonic
                        // 3 = Mean Value Weights

    // selected object id
    int _selected_object_id;

    // handle mode (bit mask)
    int _handle_mode; // 1 for creating mode, 2 for delete mode.
    
    unsigned int _draw_handle;

    // stage of the design
    STAGE _stage;

    // size of the scene
    dtype _scale;

    // size of design scale;
    dtype _design_scale;
    dtype _old_design_scale;

    unsigned int _naive_deformation_mode;

    // colors
    float _background_color[3];
    float _terminal_color[3];
    float _nonterminal_color[3];
    float _select_color[3];

    // camera
    float _camera_base_translation[3];
    float _camera_base_zoom;
    float _camera_translation[3];
    float _camera_zoom;
    float _trackball_angle[4];
    CameraConfig _camera_config;

    // saved camera configs
    std::vector<CameraConfig> _saved_camera_configs;
};