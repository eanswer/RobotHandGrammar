/**
 * @file palm_grammar.h
 * @author Jie Xu 
 */ 
#pragma once

#include "Common/Common.h"
#include "Mesh/mesh.h"
#include "Graph/node.h"

class PalmGrammarRule {
public:
    PalmGrammarRule() {}
    PalmGrammarRule(MatrixXi LHS, MatrixXi RHS, int h, int w);

    MatrixXi _LHS;
    MatrixXi _RHS;
    int _h, _w;
    Vector2i _match_location;
};

/**
* Palm Grammar is a grid grammar.
*/
class PalmGrammar {
public:
    PalmGrammar();
    PalmGrammar(std::string grammar_file);

    std::vector<std::string> line_process(std::string data);
    MatrixXi rotate_90(MatrixXi mat);

    void reset();
    void undo();
    void construct_palm_meshes(std::vector<TriMesh>& meshes, std::vector<Vector3>& colors);
    void construct_palm_meshes(std::vector<TriMesh>& meshes, std::vector<Vector3>& colors,
                                std::vector<TriMesh>& parent_meshes, 
                                std::vector<TriMesh>& child_meshes, 
                                std::vector<bool>& connector_masks,
                                std::vector<WeightHandles>& handles_list,
                                std::vector<std::pair<int, int> >& cell_coord);
    void get_colors_when_hover(std::vector<Vector3>& colors, int hover_rule_idx);

    std::vector<int> get_available_rules();
    bool apply_rule(int rule_idx);

    void get_palm_info(TriMesh& palm_mesh, ConnectionFace& connection_parent, std::vector<ConnectionFace>& connection_children);

    std::vector<PalmGrammarRule> _rules;
    int _num_rows, _num_cols;
    MatrixXi _palm_grid;
    MatrixXi _symbol_dir;
    MatrixXi _hover_mask;
    MatrixXi _cell_index;
    std::vector<MatrixXi> _palm_grid_his;
    std::vector<MatrixXi> _symbol_dir_his;
    std::set<char> _connector_symbols;

    char _starting_symbol;
    std::set<char> _terminal_symbols;
    std::map<char, TriMesh> _symbol_visual_meshes;
    std::map<char, TriMesh> _child_meshes;
    std::map<char, TriMesh> _parent_meshes;
    dtype _component_width, _component_height;

    WeightHandles _handles_template;
    WeightHandles _mount_handles_template;
    std::vector<std::pair<int, int> > _mount_connection_pairs;

    std::vector<float> _cell_widths;
    std::vector<float> _cell_heights;
    std::vector<bool> _fixed_cell_widths;
    std::vector<bool> _fixed_cell_heights;

    Vector3 _connection_origin;
    Matrix3 _connection_R;

    float _terminal_color[3];
    float _nonterminal_color[3];
    float _select_color[3];
};