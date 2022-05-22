#include "Grammar/palm_grammar.h"

PalmGrammarRule::PalmGrammarRule(MatrixXi LHS, MatrixXi RHS, int h, int w) {
    _LHS = LHS;
    _RHS = RHS;
    _h = h;
    _w = w;
    _match_location = Vector2i(-1, -1);
}

PalmGrammar::PalmGrammar() {
    _num_rows = _num_cols = 11;
    _rules.clear();
    _terminal_color[0] = 1., _terminal_color[1] = 228.0 / 255.0, _terminal_color[2] = 58.0 / 255.0;
    _nonterminal_color[0] = 0., _nonterminal_color[1] = 1., _nonterminal_color[2] = 0.;
    _select_color[0] = 1., _select_color[1] = 0., _select_color[2] = 0.;
}

std::vector<std::string> PalmGrammar::line_process(std::string data) {
    _num_rows = _num_cols = 11;

    // remove comments
    std::size_t pos = data.find_first_of('#', 0);
    if (pos != data.npos) {
        data = data.substr(0, pos);
    }
    
    // parse line
    std::vector<std::string> parsed_data;
    parsed_data.clear();
    std::string current_item = "";
    for (int i = 0;i < data.size();i++) {
        if (data[i] == ' ' || data[i] == '\n' || data[i] == '\r') {
            if (current_item.size() > 0) {
                parsed_data.push_back(current_item);
                current_item = "";
            }
        } else {
            current_item = current_item + data[i];
        }
    }

    if (current_item.size() > 0) {
        parsed_data.push_back(current_item);
    }

    return parsed_data;
}

MatrixXi PalmGrammar::rotate_90(MatrixXi mat) {
    MatrixXi mat_res(mat.cols(), mat.rows());
    for (int i = 0;i < mat_res.rows();i++)
        for (int j = 0;j < mat_res.cols();j++)
            mat_res(i, j) = mat(j, mat.cols() - 1 - i);
    return mat_res;
}

PalmGrammar::PalmGrammar(std::string grammar_file) {
    ifstream fin(grammar_file);

    std::string folder_dir = grammar_file.substr(0, grammar_file.find_last_of('/') + 1);

    _component_height = 31.;
    _component_width = 31.;

    // parse starting symbol
    std::string data;
    _starting_symbol = '#';
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 0) {
            continue;
        } else if (line_data.size() == 1 && line_data[0].size() == 1) {
            _starting_symbol = line_data[0][0];
            break;
        } else {
            throw_error("starting symbol should be of length 1.");
        }
    }

    if (_starting_symbol == '#') {
        throw_error("starting symbol not found!");
    }

    // parse number of terminal symbols
    int num_terminals = -1;
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 1) {
            num_terminals = std::stoi(line_data[0]);
            break;
        }
    }

    if (num_terminals == -1) {
        throw_error("number of terminals not found");
    }

    // parse terminal symbols
    _terminal_symbols.clear();
    for (int i = 0;i < num_terminals;i++) {
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && line_data[0].size() == 1) {
                _terminal_symbols.insert(line_data[0][0]);
                break;
            }
        }
    }

    if (_terminal_symbols.size() != num_terminals) {
        throw_error("number of terminals does not match");
    }

    // parse number of visual symbols
    int num_visual_symbols = -1;
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 1) {
            num_visual_symbols = std::stoi(line_data[0]);
            break;
        }
    }

    if (num_visual_symbols == -1) {
        throw_error("number of visual symboles not found");
    }

    // connection surface for knuckles
    _connection_origin = Vector3(0., 41.975, 0.);
    _connection_R.col(0) = Vector3::UnitX();
    _connection_R.col(1) = Vector3::UnitZ();
    _connection_R.col(2) = -Vector3::UnitY();
    _connection_origin = Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitZ()).matrix() * Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitX()).matrix() * _connection_origin;
    _connection_R = Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitZ()).matrix() * Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitX()).matrix() * _connection_R;

    _connector_symbols.clear();
    _connector_symbols.insert('k');
    _connector_symbols.insert('n');

    // parse visual symbols' meshes
    _symbol_visual_meshes.clear();
    _parent_meshes.clear();
    _child_meshes.clear();
    for (int i = 0;i < num_visual_symbols;i++) { // terminal e has no mesh
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (((line_data.size() == 2 && _connector_symbols.find(line_data[0][0]) == _connector_symbols.end()) 
                    || (line_data.size() == 3 && _connector_symbols.find(line_data[0][0]) != _connector_symbols.end())) && line_data[0].size() == 1) {
                _symbol_visual_meshes[line_data[0][0]] = TriMesh(folder_dir + line_data[1], true);
                _parent_meshes[line_data[0][0]] = TriMesh(folder_dir + line_data[1], true);
                for (int j = 2;j < line_data.size();j++) {
                    TriMesh mesh(folder_dir + line_data[j], true);
                    _symbol_visual_meshes[line_data[0][0]].merge(mesh);
                    _child_meshes[line_data[0][0]] = TriMesh(mesh);
                    _child_meshes[line_data[0][0]].rotate(-Vector3::UnitX(), pi / 2.);
                    _child_meshes[line_data[0][0]].rotate(-Vector3::UnitZ(), pi / 2.);    
                }
                _symbol_visual_meshes[line_data[0][0]].rotate(-Vector3::UnitX(), pi / 2.);
                _symbol_visual_meshes[line_data[0][0]].rotate(-Vector3::UnitZ(), pi / 2.);
                _parent_meshes[line_data[0][0]].rotate(-Vector3::UnitX(), pi / 2.);
                _parent_meshes[line_data[0][0]].rotate(-Vector3::UnitZ(), pi / 2.);
                break;
            }
        }
    }

    // parse cage description file
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 1) {
            _handles_template.load_handle_file(folder_dir + line_data[0]);
            MatrixX new_handle_positions(_handles_template.handle_positions.rows(), 3);
            for (int i = 0;i < _handles_template.handle_positions.rows();i++) {
                Vector3 pos = _handles_template.handle_positions.row(i).transpose();
                pos = Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitZ()).matrix() * Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitX()).matrix() * pos;
                new_handle_positions.row(i) = pos.transpose();
            }
            _handles_template.set_handle_positions(new_handle_positions);
            break;
        }
    }

    // parse cage description file for mount
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 1) {
            _mount_handles_template.load_handle_file(folder_dir + line_data[0]);
            MatrixX new_handle_positions(_mount_handles_template.handle_positions.rows(), 3);
            for (int i = 0;i < _mount_handles_template.handle_positions.rows();i++) {
                Vector3 pos = _mount_handles_template.handle_positions.row(i).transpose();
                pos = Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitZ()).matrix() * Eigen::AngleAxis<dtype>(pi / 2., -Vector3::UnitX()).matrix() * pos;
                new_handle_positions.row(i) = pos.transpose();
            }
            _mount_handles_template.set_handle_positions(new_handle_positions);
            break;
        }
    }

    _mount_connection_pairs.push_back(std::make_pair(0, 0));
    _mount_connection_pairs.push_back(std::make_pair(1, 1));
    _mount_connection_pairs.push_back(std::make_pair(4, 2));
    _mount_connection_pairs.push_back(std::make_pair(5, 3));

    // parse number of rules
    int num_rules = -1;
    while (std::getline(fin, data)) {
        std::vector<std::string> line_data = line_process(data);
        
        if (line_data.size() == 1) {
            num_rules = std::stoi(line_data[0]);
            break;
        }
    }

    if (num_rules == -1) {
        throw_error("number of rules not found");
    }

    // parse rules
    for (int i = 0;i < num_rules;i++) {
        // parse {RULE} str
        bool success = false;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && line_data[0] == "{RULE}") {
                success = true;
                break;
            }
        }
        if (!success) {
            throw_error("not find {RULE} for rule " + std::to_string(i));
        }

        // parse num_rows, num_cols
        int num_rows = -1, num_cols = -1;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 2) {
                num_rows = std::stoi(line_data[0]);
                num_cols = std::stoi(line_data[1]);
                break;
            }
        }

        if (num_rows == -1 || num_cols == -1) {
            throw_error("num_rows and num_cols not found for rule " + std::to_string(i));
        }

        // parse {LHS}
        success = false;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && line_data[0] == "{LHS}") {
                success = true;
                break;
            }
        }
        if (!success) {
            throw_error("not find {LHS} for rule " + std::to_string(i));
        }

        // parse lhs
        MatrixXi lhs = MatrixXi(num_rows, num_cols);
        for (int j = 0;j < num_rows;j++) {
            std::getline(fin, data);
            std::vector<std::string> line_data = line_process(data);
            if (line_data.size() != 1 || line_data[0].size() != num_cols) {
                throw_error("lhs of rule " + std::to_string(i) + " does not match the size");
            }
            for (int k = 0;k < num_cols;k++) {
                lhs(j, k) = line_data[0][k];
            }
        }

        // parse {RHS}
        success = false;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && line_data[0] == "{RHS}") {
                success = true;
                break;
            }
        }
        if (!success) {
            throw_error("not find {RHS} for rule " + std::to_string(i));
        }

        // parse rhs
        MatrixXi rhs = MatrixXi(num_rows, num_cols);
        for (int j = 0;j < num_rows;j++) {
            std::getline(fin, data);
            std::vector<std::string> line_data = line_process(data);
            if (line_data.size() != 1 || line_data[0].size() != num_cols) {
                throw_error("rhs of rule " + std::to_string(i) + " does not match the size");
            }
            for (int k = 0;k < num_cols;k++) {
                rhs(j, k) = line_data[0][k];
            }
        }

        _rules.push_back(PalmGrammarRule(lhs, rhs, num_rows, num_cols));

        // parse {ROTATION}
        success = false;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && line_data[0] == "{ROTATION}") {
                success = true;
                break;
            }
        }
        if (!success) {
            throw_error("not find {ROTATION} for rule " + std::to_string(i));
        }

        // parse rotation
        success = false;
        bool rotation = false;
        while (std::getline(fin, data)) {
            std::vector<std::string> line_data = line_process(data);
            
            if (line_data.size() == 1 && (line_data[0] == "True" || line_data[0] == "False")) {
                success = true;
                rotation = (line_data[0] == "True");
                break;
            }
        }

        if (!success) {
            throw_error("not find ROTATION for rule " + std::to_string(i));
        }

        if (rotation) {        
            // create rules with different rotations
            if (num_rows > 1 || num_cols > 1) {
                for (int i = 0;i < 3;i++) {
                    lhs = rotate_90(lhs);
                    rhs = rotate_90(rhs);
                    _rules.push_back(PalmGrammarRule(lhs, rhs, lhs.rows(), lhs.cols()));
                }
            }
        }
    }

    std::cerr << "Load grammar from " << grammar_file << " successfully!" << std::endl;

    fin.close();

    reset();

    _terminal_color[0] = 1., _terminal_color[1] = 228.0 / 255.0, _terminal_color[2] = 58.0 / 255.0;
    _nonterminal_color[0] = 0., _nonterminal_color[1] = 1., _nonterminal_color[2] = 0.;
    _select_color[0] = 1., _select_color[1] = 0., _select_color[2] = 0.;
}

void PalmGrammar::reset() {
    _palm_grid = MatrixXi::Constant(_num_rows, _num_cols, '-');
    _symbol_dir = MatrixXi::Constant(_num_rows, _num_cols, -1);
    _palm_grid(0, _num_cols / 2) = _starting_symbol;
    _symbol_dir(0, _num_cols / 2) = 0;
    _palm_grid_his.clear();
    _palm_grid_his.push_back(_palm_grid);
    _symbol_dir_his.clear();
    _symbol_dir_his.push_back(_symbol_dir);
}

void PalmGrammar::undo() {
    if (_palm_grid_his.size() > 1) {
        _palm_grid_his.pop_back();
        _palm_grid = _palm_grid_his[_palm_grid_his.size() - 1];
        _symbol_dir_his.pop_back();
        _symbol_dir = _symbol_dir_his[_symbol_dir_his.size() - 1];
    }
}

std::vector<int> PalmGrammar::get_available_rules() {
    std::vector<int> available_rules;
    available_rules.clear();
    for (int i = 0;i < _rules.size();i++) {
        _rules[i]._match_location(0) = _rules[i]._match_location(1) = -1;
        bool success = false;
        for (int j = 0;j < _num_rows && !success;j++) 
            for (int k = 0;k < _num_cols;k++) {
                bool flag = true;    
                for (int jj = 0;jj < _rules[i]._h && flag;jj++) 
                    for (int kk = 0;kk < _rules[i]._w;kk++) {
                        if (_rules[i]._LHS(jj, kk) != '*') {
                            if (j + jj < 0 || j + jj >= _num_rows || k + kk < 0 || k + kk >= _num_cols
                                || _rules[i]._LHS(jj, kk) != _palm_grid(j + jj, k + kk)) {
                                flag = false;
                                break;
                            }
                        }
                        if (_rules[i]._RHS(jj, kk) != '*') {
                            if (j + jj < 0 || j + jj >= _num_rows || k + kk < 0 || k + kk >= _num_cols) {
                                flag = false;
                                break;
                            }
                        }
                    }
                if (flag) {
                    success = true;
                    _rules[i]._match_location = Vector2i(j, k);
                    break;
                }
            }
        if (success) {
            available_rules.push_back(i);
        }
    }

    return available_rules;
}

bool PalmGrammar::apply_rule(int rule_idx) {
    if (_rules[rule_idx]._match_location(0) == -1 || _rules[rule_idx]._match_location(1) == -1) {
        std::cerr << "cannot apply rule" << std::endl;
        return false;
    }

    int i = rule_idx;
    int j = _rules[rule_idx]._match_location(0), k = _rules[rule_idx]._match_location(1);
    bool flag = true;    
    for (int jj = 0;jj < _rules[i]._h && flag;jj++) 
        for (int kk = 0;kk < _rules[i]._w;kk++) {
            if (_rules[i]._LHS(jj, kk) != '*') {
                if (j + jj < 0 || j + jj >= _num_rows || k + kk < 0 || k + kk >= _num_cols
                    || _rules[i]._LHS(jj, kk) != _palm_grid(j + jj, k + kk)) {
                    flag = false;
                    break;
                }
            }
            if (_rules[i]._RHS(jj, kk) != '*') {
                if (j + jj < 0 || j + jj >= _num_rows || k + kk < 0 || k + kk >= _num_cols) {
                    flag = false;
                    break;
                }
            }
        }
    if (!flag) {
        std::cerr << "cannot apply rule" << std::endl;
        return false;
    }

    // determine the direction for symbols
    int dir_x[4] = {-1, 0, 1, 0};
    int dir_y[4] = {0, -1, 0, 1};
    for (int jj = 0;jj < _rules[i]._h;jj++)
        for (int kk = 0;kk < _rules[i]._w;kk++) {
            if (_rules[i]._RHS(jj, kk) != '*' && _palm_grid(_rules[i]._match_location(0) + jj, _rules[i]._match_location(1) + kk) == '-') {
                for (int dir = 0;dir < 4;dir++) {
                    int nej = jj + dir_x[dir];
                    int nek = kk + dir_y[dir];
                    if (nej >= 0 && nej < _rules[i]._h && nek >= 0 && nek < _rules[i]._w) {
                        if (_rules[i]._LHS(nej, nek) != '*' && _rules[i]._LHS(nej, nek) != '-') {
                            _symbol_dir(_rules[i]._match_location(0) + jj, _rules[i]._match_location(1) + kk) = dir;
                            break;
                        }
                    }
                }
            }
        }

    for (int jj = 0;jj <_rules[i]._h;jj++)
        for (int kk = 0;kk < _rules[i]._w;kk++) {
            if (_rules[i]._RHS(jj, kk) != '*')
                _palm_grid(_rules[i]._match_location(0) + jj, _rules[i]._match_location(1) + kk)
                    = _rules[i]._RHS(jj, kk);
        }

    _palm_grid_his.push_back(_palm_grid);
    _symbol_dir_his.push_back(_symbol_dir);

    return true;
}

void PalmGrammar::construct_palm_meshes(std::vector<TriMesh>& meshes, std::vector<Vector3>& colors) {
    int center_row = 0, center_col = _num_cols / 2;
    _cell_index = MatrixXi::Ones(_num_rows, _num_cols) * -1;
    _cell_widths.clear();
    _cell_heights.clear();
    _fixed_cell_heights.clear();
    _fixed_cell_widths.clear();
    for (int i = 0;i < _num_rows;i++) {
        _cell_heights.push_back((float)_component_height);
        bool fixed = false;
        for (int j = 0;j < _num_cols;j++)
            if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                fixed = true;
                break;
            }
        _fixed_cell_heights.push_back(fixed);
    }
    for (int j = 0;j < _num_cols;j++) {
        _cell_widths.push_back((float)_component_width);
        bool fixed = false;
        for (int i = 0;i < _num_rows;i++) {
            if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                fixed = true;
                break;
            }
        }
        _fixed_cell_widths.push_back(fixed);
    }

    for (int i = 0;i < _palm_grid.rows();i++) {
        for (int j = 0;j < _palm_grid.cols();j++) {
            if (_palm_grid(i, j) != 'e' && _palm_grid(i, j) != '-') {
                dtype center_x = (i - center_row) * _component_height;
                dtype center_y = (j - center_col) * _component_width;
                if (_symbol_visual_meshes.find(_palm_grid(i, j)) != _symbol_visual_meshes.end()) {
                    TriMesh mesh = _symbol_visual_meshes[_palm_grid(i, j)];
                    mesh.rotate(Vector3::UnitZ(), pi / 2. * _symbol_dir(i, j));
                    mesh.translate(Vector3(center_x, center_y, 0.));
                    meshes.push_back(mesh);
                    if (_terminal_symbols.find(_palm_grid(i, j)) != _terminal_symbols.end())
                        // colors.push_back(Vector3(255.0/255.0,228.0/255.0,58.0/255.0));
                        colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
                    else 
                        // colors.push_back(Vector3(0., 1., 0.));
                        colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));
                    _cell_index(i, j) = meshes.size() - 1;
                } else {
                    throw_error("Symbol visual mesh not found!");
                }
            }
        }
    }
}

void PalmGrammar::construct_palm_meshes(std::vector<TriMesh>& meshes, 
        std::vector<Vector3>& colors, 
        std::vector<TriMesh>& parent_meshes, 
        std::vector<TriMesh>& child_meshes, 
        std::vector<bool>& connector_masks,
        std::vector<WeightHandles>& handles_list,
        std::vector<std::pair<int, int> >& cell_coord) {

    int center_row = 0, center_col = _num_cols / 2;
    connector_masks.clear();
    handles_list.clear();
    cell_coord.clear();
    _cell_widths.clear();
    _cell_heights.clear();
    _fixed_cell_heights.clear();
    _fixed_cell_widths.clear();
    for (int i = 0;i < _num_rows;i++) {
        _cell_heights.push_back((float)_component_height);
        bool fixed = false;
        for (int j = 0;j < _num_cols;j++)
            if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                fixed = true;
                break;
            }
        _fixed_cell_heights.push_back(fixed);
    }
    for (int j = 0;j < _num_cols;j++) {
        _cell_widths.push_back((float)_component_width);
        bool fixed = false;
        for (int i = 0;i < _num_rows;i++) {
            if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                fixed = true;
                break;
            }
        }
        _fixed_cell_widths.push_back(fixed);
    }

    _cell_index = MatrixXi::Ones(_num_rows, _num_cols) * -1;
    for (int i = 0;i < _palm_grid.rows();i++) {
        for (int j = 0;j < _palm_grid.cols();j++) {
            if (_palm_grid(i, j) != 'e' && _palm_grid(i, j) != '-') {
                dtype center_x = (i - center_row) * _component_height;
                dtype center_y = (j - center_col) * _component_width;
                if (_symbol_visual_meshes.find(_palm_grid(i, j)) != _symbol_visual_meshes.end()) {
                    TriMesh mesh = _symbol_visual_meshes[_palm_grid(i, j)];
                    mesh.rotate(Vector3::UnitZ(), pi / 2. * _symbol_dir(i, j));
                    mesh.translate(Vector3(center_x, center_y, 0.));
                    meshes.push_back(mesh);
                    if (_terminal_symbols.find(_palm_grid(i, j)) != _terminal_symbols.end())
                        colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
                    else 
                        colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));

                    if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                        connector_masks.push_back(true);
                    } else {
                        connector_masks.push_back(false);
                    }

                    if (_connector_symbols.find(_palm_grid(i, j)) != _connector_symbols.end()) {
                        TriMesh child_mesh = _child_meshes[_palm_grid(i, j)];
                        child_mesh.rotate(Vector3::UnitZ(), pi / 2. * _symbol_dir(i, j));
                        child_mesh.translate(Vector3(center_x, center_y, 0.));
                        child_meshes.push_back(child_mesh);
                    } else {
                        child_meshes.push_back(TriMesh());
                    }

                    {
                        TriMesh parent_mesh = _parent_meshes[_palm_grid(i, j)];
                        parent_mesh.rotate(Vector3::UnitZ(), pi / 2. * _symbol_dir(i, j));
                        parent_mesh.translate(Vector3(center_x, center_y, 0.));
                        parent_meshes.push_back(parent_mesh);
                    }

                    WeightHandles handles = _handles_template;
                    MatrixX new_handle_positions(handles.handle_positions.rows(), 3);
                    for (int k = 0;k < handles.handle_positions.rows();k++) {
                        Vector3 pos = handles.handle_positions.row(k).transpose();
                        pos = Eigen::AngleAxis<dtype>(pi / 2. * _symbol_dir(i, j), Vector3::UnitZ()).matrix() * pos + Vector3(center_x, center_y, 0.);
                        new_handle_positions.row(k) = pos.transpose();
                    }
                    handles.set_handle_positions(new_handle_positions);
                    handles_list.push_back(handles);

                    cell_coord.push_back(std::make_pair(i, j));

                    _cell_index(i, j) = meshes.size() - 1;
                } else {
                    throw_error("Symbol visual mesh not found!");
                }
            }
        }
    }
}

void PalmGrammar::get_colors_when_hover(std::vector<Vector3>& colors, int hover_rule_idx) {
    if (hover_rule_idx != -1) {
        int j = _rules[hover_rule_idx]._match_location(0), k = _rules[hover_rule_idx]._match_location(1);

        colors.clear();
        for (int jj = 0;jj < _palm_grid.rows();jj++) {
            for (int kk = 0;kk < _palm_grid.cols();kk++) {
                if (_palm_grid(jj, kk) != 'e' && _palm_grid(jj, kk) != '-') {
                    if (jj >= j && jj < j + _rules[hover_rule_idx]._h && kk >= k && kk < k + _rules[hover_rule_idx]._w) {
                        colors.push_back(Vector3(_select_color[0], _select_color[1], _select_color[2]));
                    }
                    else if (_terminal_symbols.find(_palm_grid(jj, kk)) != _terminal_symbols.end()) {
                        colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
                    } else {
                        colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));
                    }
                }
            }
        }
    } else {
        colors.clear();
        for (int jj = 0;jj < _palm_grid.rows();jj++) {
            for (int kk = 0;kk < _palm_grid.cols();kk++) {
                if (_palm_grid(jj, kk) != 'e' && _palm_grid(jj, kk) != '-') {
                    if (_terminal_symbols.find(_palm_grid(jj, kk)) != _terminal_symbols.end()) {
                        colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
                    } else {
                        colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));
                    }
                }
            }
        }
    }
}

void PalmGrammar::get_palm_info(
    TriMesh& palm_mesh, 
    ConnectionFace& connection_parent, 
    std::vector<ConnectionFace>& connection_children) {

    std::vector<TriMesh> meshes;
    std::vector<Vector3> colors;
    construct_palm_meshes(meshes, colors);

    palm_mesh = meshes[0];
    for (int i = 1;i < meshes.size();i++)
        palm_mesh.merge(meshes[i]);

    connection_parent._origin = Vector3::Zero();
    connection_parent._R = Matrix3::Identity();

    int center_row = 0, center_col = _num_cols / 2;
    connection_children.clear();
    for (int jj = 0;jj < _palm_grid.rows();jj++)
        for (int kk = 0;kk < _palm_grid.cols();kk++) {
            // TODO: change here if need to add more connectors
            if (_connector_symbols.find(_palm_grid(jj, kk)) != _connector_symbols.end()) {
                dtype center_x = (jj - center_row) * _component_height;
                dtype center_y = (kk - center_col) * _component_width;
                // Vector3d origin(16., 0., 0.);

                // Eigen::Matrix3d R;
                // R.col(0) = Vector3::UnitX();
                // R.col(1) = Vector3::UnitZ();
                // R.col(2) = Vector3::UnitY();
                // Eigen::Matrix3d extra_rotation = Eigen::AngleAxis<double>(pi / 2. * _symbol_dir(jj, kk), Vector3::UnitZ()).matrix();

                // origin = extra_rotation * origin + Vector3(center_x, center_y, 0.);
                // R = extra_rotation * R;

                Eigen::Matrix3d extra_rotation = Eigen::AngleAxis<double>(pi / 2. * _symbol_dir(jj, kk), Vector3::UnitZ()).matrix();
                Vector3 origin = _connection_origin + Vector3(center_x, center_y, 0.);
                Matrix3 R = extra_rotation * _connection_R;
                
                connection_children.push_back(ConnectionFace(origin, R));
            }
        }
}