#include "linear_blend_skinning_ui.h"

void CameraConfig::get(igl::opengl::glfw::Viewer* viewer) {
    for (int k = 0;k < 3;k++) {
        _camera_base_translation[k] = viewer->core().camera_base_translation[k];
        _camera_translation[k] = viewer->core().camera_translation[k];
    }
    _camera_base_zoom = viewer->core().camera_base_zoom;
    _camera_zoom = viewer->core().camera_zoom;
    for (int k = 0;k < 4;k++)
        _trackball_angle[k] = viewer->core().trackball_angle.coeffs()[k];
}

void CameraConfig::set(igl::opengl::glfw::Viewer* viewer) {
    for (int k = 0;k < 3;k++) {
        viewer->core().camera_base_translation[k] = _camera_base_translation[k];
        viewer->core().camera_translation[k] = _camera_translation[k];
    }
    viewer->core().camera_base_zoom = _camera_base_zoom;
    viewer->core().camera_zoom = _camera_zoom;
    for (int k = 0;k < 4;k++)
        viewer->core().trackball_angle.coeffs()[k] = _trackball_angle[k];
}

void CameraConfig::save_to_file(std::string path) {
    FILE* fp = fopen(path.c_str(), "w");
    for (int k = 0;k < 3;k++)
        fprintf(fp, "%.6f ", _camera_base_translation[k]);
    fprintf(fp, "\n");
    for (int k = 0;k < 3;k++)
        fprintf(fp, "%.6f ", _camera_translation[k]);
    fprintf(fp, "\n");
    fprintf(fp, "%.6f\n", _camera_base_zoom);
    fprintf(fp, "%.6f\n", _camera_zoom);
    for (int k = 0;k < 4;k++)
        fprintf(fp, "%.6f ", _trackball_angle[k]);
    fprintf(fp, "\n");
    fclose(fp);
}

void CameraConfig::load_from_file(std::string path) {
    FILE* fp = fopen(path.c_str(), "r");
    for (int k = 0;k < 3;k++)
        fscanf(fp, "%f ", &_camera_base_translation[k]);
    for (int k = 0;k < 3;k++)
        fscanf(fp, "%f ", &_camera_translation[k]);
    fscanf(fp, "%f", &_camera_base_zoom);
    fscanf(fp, "%f", &_camera_zoom);
    for (int k = 0;k < 4;k++)
        fscanf(fp, "%f ", &_trackball_angle[k]);
    fclose(fp);
}

LinearBlendSkinningUI::LinearBlendSkinningUI() :
    _draw_handle(1),
    _weight_type(0),
    menu()
{
    _objects.clear();
    _selected_object_id = -1;
    _palm_grammar = PalmGrammar(PROJECT_SOURCE_DIR "/data/grammar/palm_grammar/palm_grammar.txt");
    _finger_grammar = FingerGrammar(PROJECT_SOURCE_DIR "/data/grammar/finger_grammar/grammar.txt");
    
    _finger_grammar.print_grammar();
    
    menu.set_menu_name("Design Panel");
    _stage = PALM;
    _hover_rule_id = -1;
    _naive_deformation_mode = 1;
    _background_color[0] = 1., _background_color[1] = 1., _background_color[2] = 1.;
    _terminal_color[0] = 1., _terminal_color[1] = 228.0 / 255.0, _terminal_color[2] = 58.0 / 255.0;
    _nonterminal_color[0] = 0., _nonterminal_color[1] = 1., _nonterminal_color[2] = 0.;
    _select_color[0] = 1., _select_color[1] = 0., _select_color[2] = 0.;

    _design_scale = _old_design_scale = 1.;

    // load camera configs
    std::string path = PROJECT_SOURCE_DIR "/data/camera_config/";
    for (const auto & entry : std::experimental::filesystem::directory_iterator(path)) {
        CameraConfig config;
        config.load_from_file(entry.path().generic_string());
        _saved_camera_configs.push_back(config);
    }
}

void LinearBlendSkinningUI::init(igl::opengl::glfw::Viewer* _viewer) {
    igl::opengl::glfw::ViewerPlugin::init(_viewer);
    viewer->core().set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);

    viewer->plugins.push_back(&menu);

#ifdef WIN32
    menu.callback_draw_viewer_window = [&]() {draw_menu(); };
#else
    menu.callback_draw_viewer_menu = [&]() {draw_menu(); };
#endif

    viewer->data().show_lines = false;
    
    reset_palm();
    update_viewer();    
}

void LinearBlendSkinningUI::reset_palm() {
    _palm_grammar.reset();
    construct_meshes();
}

void LinearBlendSkinningUI::reset_finger() {
    _finger_grammar.reset();

    construct_meshes();
}

void LinearBlendSkinningUI::construct_meshes() {
    _meshes.clear();
    _colors.clear();
    if (_stage == PALM) {
        _palm_grammar.construct_palm_meshes(_meshes, _colors);
    } else if (_stage >= FINGER) {
        construct_finger_objects(_finger_grammar._graph);
    }
}

void LinearBlendSkinningUI::construct_finger_objects(Tree &graph) {
    graph.get_deformation_objects(_objects, _object_types, &_palm_grammar);
    for (int i = 0;i < _objects.size();i++) {
        TriMesh mesh(_objects[i]->_V[0].transpose(), _objects[i]->_F[0].transpose());
        for (int j = 1;j < _objects[i]->_V.size();j++)
            mesh.merge(_objects[i]->_V[j].transpose(), _objects[i]->_F[j].transpose());
        _meshes.push_back(mesh);
        if (_object_types[i] == Node::NodeType::TERMINAL)
            _colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
        else
            _colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));
    }
}

void LinearBlendSkinningUI::update_object_mesh(
    std::vector<MatrixX3d>& V, 
    std::vector<MatrixX4i>& T, 
    std::vector<MatrixX3i>& F) {

    int num_vertices = _V.rows();
    int num_tets = _T.rows();
    int num_faces = _F.rows();
    
    for (int i = 0;i < V.size();i++) {
        _V.conservativeResize(num_vertices + V[i].rows(), 3);
        _V.block(num_vertices, 0, V[i].rows(), 3) = V[i];

        _T.conservativeResize(num_tets + T[i].rows(), 4);
        _T.block(num_tets, 0, T[i].rows(), 4) = T[i] + MatrixXi::Ones(T[i].rows(), 4) * num_vertices;

        _F.conservativeResize(num_faces + F[i].rows(), 3);
        _F.block(num_faces, 0, F[i].rows(), 3) = F[i] + MatrixXi::Ones(F[i].rows(), 3) * num_vertices;

        num_vertices += V[i].rows();
        num_tets += T[i].rows();
        num_faces += F[i].rows();
    }

    _scale = max(_V.col(0).maxCoeff() - _V.col(0).minCoeff(), 
                max(_V.col(1).maxCoeff() - _V.col(1).minCoeff(), _V.col(2).maxCoeff() - _V.col(2).minCoeff()));
}

void LinearBlendSkinningUI::update_viewer_mesh() {
    if (_meshes.size() > 0) {
        TriMesh mesh = _meshes[0];
        for (int i = 1;i < _meshes.size();i++) {
            mesh.merge(_meshes[i]);
        }
        _V = mesh._V.transpose();
        _F = mesh._F.transpose();
        _vert_colors = MatrixX3d(mesh._V.cols(), 3);
        int idx = 0;
        for (int i = 0;i < _meshes.size();i++) {
            for (int j = 0;j < _meshes[i]._V.cols();j++)
                _vert_colors.row(idx++) = _colors[i].transpose();
        }

        _scale = max(_V.col(0).maxCoeff() - _V.col(0).minCoeff(), 
                max(_V.col(1).maxCoeff() - _V.col(1).minCoeff(), _V.col(2).maxCoeff() - _V.col(2).minCoeff()));
    }
}

void LinearBlendSkinningUI::update_colors() {
    if (_stage == PALM) {
        _palm_grammar.get_colors_when_hover(_colors, _hover_rule_id);
    } else if (_stage == FINGER) {
        _colors.clear();
        for (int i = 0;i < _objects.size();i++)
            if (_object_types[i] == Node::NodeType::TERMINAL)
                _colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
            else
                _colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));

        if (_hover_rule_id != -1) {
            Node* application_node = _finger_grammar._graph.find_node(_finger_grammar._rules[_hover_rule_id]);
            if (application_node->_symbol != "P" && application_node->_symbol != "palm") {
                _colors[application_node->_id_in_objects](0) = _select_color[0];
                _colors[application_node->_id_in_objects](1) = _select_color[1];
                _colors[application_node->_id_in_objects](2) = _select_color[2];
            }
        }

    } else if (_stage == DEFORMATION) {
        _colors.clear();
        for (int i = 0;i < _objects.size();i++) {
            if (i == _selected_object_id)
                _colors.push_back(Vector3(1., 0., 0.));
            else if (_object_types[i] == Node::NodeType::TERMINAL)
                _colors.push_back(Vector3(_terminal_color[0], _terminal_color[1], _terminal_color[2]));
            else
                _colors.push_back(Vector3(_nonterminal_color[0], _nonterminal_color[1], _nonterminal_color[2]));
        }
    }
}

void LinearBlendSkinningUI::update_viewer(bool update_mesh) {
    update_viewer_mesh();

    // add mesh
    if (update_mesh) {
        viewer->data().clear();
        viewer->data().set_mesh(_V, _F);
    }
    viewer->data().set_colors(_vert_colors);

    // draw vertices
    draw_mesh();
    // draw handles
    draw_handles();
}

void LinearBlendSkinningUI::start_finger_grammar_stage() {
    _finger_grammar._starting_symbol = "P";

    // infer node information from palm grammar
    Node::NodeType symbol_type = Node::NodeType::NONTERMINAL;
    std::string symbol = "P";
    Node::PartType part_type = Node::PartType::BODY;
    Node::MeshType mesh_type = Node::MeshType::SINGLE;
    TriMesh palm_mesh;
    ConnectionFace connection_parent;
    std::vector<ConnectionFace> connection_children;
    _palm_grammar.get_palm_info(palm_mesh, connection_parent, connection_children);
    std::vector<TriMesh> meshes;

    meshes.push_back(palm_mesh);
    _finger_grammar._template_nodes["P"] = 
        Node(Node::NodeType::NONTERMINAL, "P", 
            Node::PartType::BODY, std::make_pair(Vector3d::Zero(), Vector3d::Zero()),
            Node::MeshType::SINGLE, meshes,
            nullptr,
            connection_parent, connection_children,
            DeformationType::SCALE, std::vector<Eigen::Vector3d>(),
            WeightHandles());
    _finger_grammar._template_nodes["palm"] = 
        Node(Node::NodeType::TERMINAL, "palm", 
            Node::PartType::BODY, std::make_pair(Vector3d::Zero(), Vector3d::Zero()),
            Node::MeshType::SINGLE, meshes,
            nullptr,
            connection_parent, connection_children,
            DeformationType::SCALE, std::vector<Eigen::Vector3d>(),
            WeightHandles());
    
    _finger_grammar.finalize_rule_nodes();

    _finger_grammar.reset();

    _finger_grammar._graph.print_graph();

    _stage = FINGER;
}

void LinearBlendSkinningUI::start_deformation_stage() {
    _stage = DEFORMATION;
    _hover_rule_id = -1;
    for (int i = 0;i < _objects.size();i++)
        _objects[i]->compute_weights();

    if (_naive_deformation_mode == 1) 
        for (int i = 0;i < _objects.size();i++)
            _objects[i]->reset_naive_mode_parameters();
}

void LinearBlendSkinningUI::scale_design() {
    _design_scale = max(_design_scale, 0.1);
    dtype scale_ratio = _design_scale / _old_design_scale;
    for (int i = 0;i < _objects.size();i++)
        _objects[i]->scale_design(scale_ratio);
    for (int i = 0;i < _palm_grammar._cell_widths.size();i++)
        _palm_grammar._cell_widths[i] *= scale_ratio;
    for (int i = 0;i < _palm_grammar._cell_heights.size();i++)
        _palm_grammar._cell_heights[i] *= scale_ratio;
    _old_design_scale = _design_scale;
    update_viewer(true);
}

// deform the mesh of each object
MatrixX3d LinearBlendSkinningUI::transform_mesh() {
    int start_idx = 0;
    MatrixX3d V(_V);
    for (int i = 0;i < _objects.size();i++) {
        std::vector<MatrixX3d> sub_V = _objects[i]->transform_mesh(_naive_deformation_mode);
        for (int j = 0;j < sub_V.size();j++) {
            V.block(start_idx, 0, sub_V[j].rows(), sub_V[j].cols()) = sub_V[j];
            start_idx += sub_V[j].rows();
        }
    }
    return V;
}

void LinearBlendSkinningUI::draw_mesh() {
    MatrixX3d vertices = transform_mesh();
    viewer->data().set_vertices(vertices);
}

void LinearBlendSkinningUI::draw_handles() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd points, point_colors, line_colors;
    Eigen::MatrixXi lines;
    
    viewer->data().point_size = 10;

    if (_draw_handle == 1) {
        if (_selected_object_id != -1 && _objects[_selected_object_id]->_handles.positions().rows() > 0) {
            _objects[_selected_object_id]->visualize_handles(V, F, points, point_colors, lines, line_colors, true, _naive_deformation_mode);
            if (_moving_handle_id >= 0) {
                point_colors.row(_moving_handle_id) = Eigen::RowVector3d(0.0, 1.0, 0.0);
            }
        }
    }
    viewer->data().clear_labels();
    viewer->data().set_points(points, point_colors);
    viewer->data().set_edges(points, lines, line_colors);
}

void LinearBlendSkinningUI::draw_menu() {
    // ImGui::ShowDemoWindow();

    // Helper for making checkboxes
    auto make_checkbox = [&](const char* label, unsigned int& option)
    {
        return ImGui::Checkbox(label,
            [&]() { return viewer->core().is_set(option); },
            [&](bool value) { return viewer->core().set(option, value); }
        );
    };

    // set camera
    if (_stage == PALM) {
        viewer->core().camera_base_translation = Vector3f(0., 0.0, 1.5);
        viewer->core().camera_base_zoom = 0.064516;
        viewer->core().camera_translation = Vector3f(-145., -11.3, 1.25);
        viewer->core().orthographic = true;
        viewer->core().rotation_type = igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION;
        viewer->core().trackball_angle.coeffs() = Vector4f(- sqrt(2.) / 2., - sqrt(2.) / 2., 0., 0.);

        viewer->core().camera_zoom = 0.19;
    } else if (_stage >= FINGER) {
        // viewer->core().orthographic = true;
        viewer->core().orthographic = false;
        viewer->core().rotation_type = igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL;
    }

    if (_stage == PALM) {
        menu.set_menu_name("Design Stage I: Palm Grammar");
        
        char buf[100];
        sprintf(buf, "Palm Grammar");
        ImGui::Text(buf);
        if (ImGui::Button("Load Palm Grammar")) {
            auto filename = igl::file_dialog_open();
            if (filename.size() > 0) {
                _palm_grammar = PalmGrammar(filename);
                reset_palm();
                update_viewer();
            }
        }

        ImGui::Separator();

        {
            ImGui::Text("Applicable Rules");
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            
            int rule_per_row = 2;
            ImVec2 button_sz((window_width - style.ItemSpacing.x * (rule_per_row + 1)) / rule_per_row, 40);
            std::vector<int> available_rules = _palm_grammar.get_available_rules();
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(3. / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(3.0 / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(3.0 / 7.0f, 0.8f, 0.8f));

            int cnt = 0;
            int hover_rule_id = -1;
            for (int rule_id : available_rules) {
                MatrixXi lhs = _palm_grammar._rules[rule_id]._LHS;
                MatrixXi rhs = _palm_grammar._rules[rule_id]._RHS;
                
                std::string rule_description;
                for (int i = lhs.rows() - 1;i >= 0;i--) {
                    for (int j = 0;j < lhs.cols();j++)
                        rule_description = rule_description + char(lhs(i, j));
                    if (i == lhs.rows() / 2) {
                        rule_description = rule_description + " -> ";
                    } else {
                        rule_description = rule_description + "    ";
                    }
                    for (int j = 0;j < rhs.cols();j++)
                        rule_description = rule_description + char(rhs(i, j));
                    if (i > 0) 
                        rule_description = rule_description + '\n';
                }
                
                if (ImGui::Button(rule_description.c_str(), button_sz)) {
                    _palm_grammar.apply_rule(rule_id);
                    
                    construct_meshes();
                    update_viewer();
                }

                if (ImGui::IsItemHovered()) {
                    hover_rule_id = rule_id;
                }
            
                if (cnt % rule_per_row < rule_per_row - 1 && cnt < available_rules.size() - 1)
                    ImGui::SameLine();

                cnt += 1;
            }

            if (_hover_rule_id != hover_rule_id) {
                _hover_rule_id = hover_rule_id;
                update_colors();
                update_viewer(false);
            }

            ImGui::PopStyleColor(3);
        }

        ImGui::Separator();

        {
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            ImVec2 button_sz((window_width - style.ItemSpacing.x * 3) / 2, 20);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
            if (ImGui::Button("reset", button_sz)) {
                reset_palm();
                update_viewer();
            }
            ImGui::PopStyleColor(3);
            ImGui::SameLine();
            if (ImGui::Button("undo", button_sz)) {
                _palm_grammar.undo();
                construct_meshes();
                update_viewer();
            }
            
        }

        ImGui::Separator();

        {
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            ImVec2 button_sz(window_width - 2. * style.ItemSpacing.x, 20);
            if (_palm_grammar.get_available_rules().size() == 0) {
                ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(2. / 7.f, 0.6f, 0.6f));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(2. / 7.f, 0.7f, 0.7f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(2. / 7.f, 0.8f, 0.8f));
                if (ImGui::Button("finish", button_sz)) {
                    start_finger_grammar_stage();
                    // _finger_grammar.construct_hand_meshes(_meshes, _colors);
                    construct_meshes();
                    update_viewer();
                }
            } else {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::Button("finish", button_sz);
            }
            ImGui::PopStyleColor(3);
        }

    } else if (_stage == FINGER) {
        menu.set_menu_name("Design Stage II: Finger Grammar");
        
        char buf[100];
        sprintf(buf, "Finger Grammar");
        ImGui::Text(buf);
        if (ImGui::Button("Load Finger Grammar")) {
            auto filename = igl::file_dialog_open();
            if (filename.size() > 0) {
                _finger_grammar = FingerGrammar(filename);
                reset_finger();
                update_viewer();
            }
        }

        ImGui::Separator();

        {
            ImGui::Text("Applicable Rules");

            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            
            int rule_per_row = 1;
            ImVec2 button_sz((window_width - style.ItemSpacing.x * (rule_per_row + 1)) / rule_per_row, 45);

            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(3. / 7.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(3.0 / 7.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(3.0 / 7.0f, 0.8f, 0.8f));

            std::vector<int> available_rules = _finger_grammar.get_available_rules();
            int hover_rule_id = -1;
            for (int rule_id : available_rules) {
                std::string lhs = _finger_grammar._rules[rule_id]._LHS;
                std::string rhs = "";
                for (int i = 0;i < _finger_grammar._rules[rule_id]._RHS.size();i++)
                    rhs = rhs + ' ' + _finger_grammar._rules[rule_id]._RHS_nodes[i]._symbol;
                char rule_description[200];
                // sprintf(rule_description, "apply rule %d: %s -> %s", rule_id, lhs.c_str(), rhs.c_str());
                sprintf(rule_description, "%s -> %s", lhs.c_str(), rhs.c_str());
                if (ImGui::Button(rule_description)) {
                    _finger_grammar.apply_rule(rule_id);
                    _finger_grammar._graph.print_graph();
                    construct_meshes();
                    update_viewer();
                }

                if (ImGui::IsItemHovered()) {
                    hover_rule_id = rule_id;
                }
            }

            if (_hover_rule_id != hover_rule_id) {
                _hover_rule_id = hover_rule_id;
                update_colors();
                update_viewer(false);
            }

            ImGui::PopStyleColor(3);
        }

        ImGui::Separator();

        {
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            ImVec2 button_sz((window_width - style.ItemSpacing.x * 3) / 2, 20);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
            if (ImGui::Button("reset", button_sz)) {
                reset_finger();
                update_viewer();
            }
            ImGui::PopStyleColor(3);
            ImGui::SameLine();
            if (ImGui::Button("undo", button_sz)) {
                _finger_grammar.undo();
                _finger_grammar._graph.print_graph();
                construct_meshes();
                update_viewer();
            }
        }

        ImGui::Separator();

        {
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            ImVec2 button_sz(window_width - 2. * style.ItemSpacing.x, 20);
            if (_finger_grammar.get_available_rules().size() == 0) {
                ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(2. / 7.f, 0.6f, 0.6f));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(2. / 7.f, 0.7f, 0.7f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(2. / 7.f, 0.8f, 0.8f));
                if (ImGui::Button("finish", button_sz)) {
                    construct_meshes();
                    start_deformation_stage();
                    update_viewer();
                }
            } else {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.7, 0.7, 0.7, 1.0));
                ImGui::Button("finish", button_sz);
            }
            ImGui::PopStyleColor(3);
        }
    } else {
        menu.set_menu_name("Design Stage III: Shape Deformation");
        if (_naive_deformation_mode) {
            if (_selected_object_id != -1) {
                if (_objects[_selected_object_id]->_symbol == "fork") {

                } else if (_objects[_selected_object_id]->_symbol == "palm") {
                    int cell_row = _objects[_selected_object_id]->_cell_row;
                    int cell_col = _objects[_selected_object_id]->_cell_col;
                    if (!_palm_grammar._fixed_cell_widths[cell_col])
                        ImGui::DragFloat("cell width", &(_palm_grammar._cell_widths[cell_col]), 0.1, 5., 50.0);
                    if (!_palm_grammar._fixed_cell_heights[cell_row])
                        ImGui::DragFloat("cell height", &(_palm_grammar._cell_heights[cell_row]), 0.1, 5., 50.0);
                } else if (_objects[_selected_object_id]->_symbol == "phalanx" || _objects[_selected_object_id]->_symbol == "tip"
                            || _objects[_selected_object_id]->_symbol == "solid" || _objects[_selected_object_id]->_symbol == "adapter") {
                    
                    ImGui::DragFloat("length", &_objects[_selected_object_id]->_length, 0.1, 5., 50.0);
                    ImGui::DragFloat("shear x", &_objects[_selected_object_id]->_shear_x, 0.1, -10.0, 10.0);
                    ImGui::DragFloat("shear y", &_objects[_selected_object_id]->_shear_y, 0.1, -10.0, 10.0);
                    
                    if (_objects[_selected_object_id]->_children.size() == 0 || _objects[_selected_object_id]->_children[0]->_symbol != "fork") {
                        ImGui::DragFloat("face width", &_objects[_selected_object_id]->_child_face_width, 0.1, 5., 50.0);
                        ImGui::DragFloat("face height", &_objects[_selected_object_id]->_child_face_height, 0.1, 5., 50.0);
                    }
                } else if (_objects[_selected_object_id]->_symbol == "joint") {
                    if (_objects[_selected_object_id]->_parent != nullptr && _objects[_selected_object_id]->_parent->_symbol != "palm" 
                        && _objects[_selected_object_id]->_parent->_symbol != "fork" 
                        && (_objects[_selected_object_id]->_children.size() == 0 || _objects[_selected_object_id]->_children[0]->_symbol != "folk")) {
                        ImGui::DragFloat("length", &_objects[_selected_object_id]->_length, 0.1, 5., 50.0);
                        ImGui::DragFloat("face width", &_objects[_selected_object_id]->_child_face_width, 0.1, 5., 50.0);
                    }
                } else {
                    throw_error("unrecognized symbol type for deformation: " + _objects[_selected_object_id]->_symbol);
                }
                _objects[_selected_object_id]->update_naive_deformation(_palm_grammar, _objects);
                draw_mesh();
                draw_handles();
            }
            ImGui::Separator();
            if (ImGui::InputDouble("scale", &_design_scale, 0.01f, 1.0f, "%.4f")) {
                scale_design();
            }
        }

        ImGui::Separator();

        {
            ImGui::LabelText("", "design name: ");
            ImGui::SameLine();
            ImGui::InputText("", _design_name);
            ImGuiStyle& style = ImGui::GetStyle();
            float window_width = ImGui::GetWindowWidth();
            ImVec2 button_sz(window_width - 2. * style.ItemSpacing.x, 20);
            ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(2. / 7.f, 0.6f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(2. / 7.f, 0.7f, 0.7f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(2. / 7.f, 0.8f, 0.8f));
            if (ImGui::Button("Export Design", button_sz)) {
                export_design(_design_name);
            }
            ImGui::PopStyleColor(3);
        }
    }

    ImGui::NewLine();

    unsigned int _draw_handle_old = _draw_handle;
    if (ImGui::CollapsingHeader("Options")) {
        make_checkbox("Wireframe", viewer->data().show_lines);
        make_checkbox("Fill", viewer->data().show_faces);
        make_checkbox("Naive Mode", _naive_deformation_mode);
    }

    if (ImGui::CollapsingHeader("Colors")) {
        ImGuiStyle& style = ImGui::GetStyle();
        float window_width = ImGui::GetWindowWidth();
        ImGui::PushItemWidth(window_width - 2. * style.ItemSpacing.x - 90);

        ImGui::ColorEdit3("background", _background_color);

        viewer->core().background_color(0) = _background_color[0];
        viewer->core().background_color(1) = _background_color[1];
        viewer->core().background_color(2) = _background_color[2];

        ImGui::ColorEdit3("nonterminal", _nonterminal_color);

        ImGui::ColorEdit3("terminal", _terminal_color);

        ImGui::ColorEdit3("select item", _select_color);
        
        bool changed = false;
        for (int ii = 0;ii < 3;ii++) {
            // std::cerr << fabs(_palm_grammar._select_color[ii] - _select_color[ii]) << std::endl;
            if (fabs(_palm_grammar._nonterminal_color[ii] - _nonterminal_color[ii]) > 1e-4) {
                changed = true;
            }
            if (fabs(_palm_grammar._terminal_color[ii] - _terminal_color[ii]) > 1e-4) {
                changed = true;
            }
            if (fabs(_palm_grammar._select_color[ii] - _select_color[ii]) > 1e-4) {
                changed = true;
            }
            _palm_grammar._nonterminal_color[ii] = _nonterminal_color[ii];
            _palm_grammar._terminal_color[ii] = _terminal_color[ii];
            _palm_grammar._select_color[ii] = _select_color[ii];
        }

        if (changed) {
            update_colors();
            update_viewer(false);
        }

        ImGui::PopItemWidth();
    }

    if (ImGui::CollapsingHeader("Camera")) {
        ImGuiStyle& style = ImGui::GetStyle();
        float window_width = ImGui::GetWindowWidth();
        ImGui::PushItemWidth(window_width - 2. * style.ItemSpacing.x - 120);

        _camera_config.get(viewer);

        ImGui::InputFloat3("base translation", _camera_config._camera_base_translation);
        ImGui::InputFloat3("translation", _camera_config._camera_translation);
        ImGui::InputFloat("base zoom", &_camera_config._camera_base_zoom);
        ImGui::InputFloat("zoom", &_camera_config._camera_zoom);
        ImGui::InputFloat4("angle", _camera_config._trackball_angle);

        ImGui::PopItemWidth();

        ImVec2 button_sz(window_width - 2. * style.ItemSpacing.x, 20);
        ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(2. / 7.f, 0.6f, 0.6f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(2. / 7.f, 0.7f, 0.7f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(2. / 7.f, 0.8f, 0.8f));
        if (ImGui::Button("save camera", button_sz)) {
            std::string camera_config_folder = PROJECT_SOURCE_DIR "/data/camera_config/";
            if (!std::experimental::filesystem::exists(camera_config_folder.c_str())) {
                std::experimental::filesystem::create_directories(camera_config_folder.c_str());
            }
            std::string camera_config_path = camera_config_folder + "/" + std::to_string(_saved_camera_configs.size()) + ".txt";
            _camera_config.save_to_file(camera_config_path);
            _saved_camera_configs.push_back(_camera_config);
        }
    
        ImGui::PopStyleColor(3);

        {   
            ImVec2 button_sz((window_width - 3. * style.ItemSpacing.x) / 2., 20);
            
            ImGui::Separator();
            for (int k = 0;k < _saved_camera_configs.size();k++) {
                std::string text = "camera " + std::to_string(k);
                if (ImGui::Button(text.c_str(), button_sz)) {
                    _saved_camera_configs[k].set(viewer);
                }    
                if (k % 2 == 0) {
                    ImGui::SameLine();
                }
            }
        }
    }
}

int LinearBlendSkinningUI::get_closest_mesh_vertex() {
    int index = -1;
    int fid;
    Eigen::Vector3f bc;
    double x = viewer->current_mouse_x;
    double y = viewer->core().viewport(3) - viewer->current_mouse_y;
    Eigen::RowVector3f last_mouse(x, y, 0);
    if (igl::unproject_onto_mesh(
        last_mouse.head(2),
        viewer->core().view,
        viewer->core().proj,
        viewer->core().viewport,
        viewer->data().V,
        viewer->data().F,
        fid,
        bc))
    {
        const Eigen::MatrixX3i& F = viewer->data().F;
        int coord;
        bc.maxCoeff(&coord);
        index = F(fid, coord);
    }
    return index;
}

int LinearBlendSkinningUI::get_closest_handle_id(int object_id) {
    double x = viewer->current_mouse_x;
    double y = viewer->core().viewport(3) - viewer->current_mouse_y;
    Eigen::RowVector3f last_mouse(x, y, 0);
    int handle_id = _objects[object_id]->_handles.select_handle(
                            last_mouse.head(2),
                            viewer->core().view,
                            viewer->core().proj,
                            viewer->core().viewport,
                            _scale / 50.,
                            _naive_deformation_mode);
    return handle_id;
}

int LinearBlendSkinningUI::query_object_by_vertex(int vertex_id) {
    for (int i = 0;i < _objects.size();i++) {
        for (int j = 0;j < _objects[i]->_V.size();j++) {
            if (vertex_id < _objects[i]->_V[j].rows()) {
                return i;
            }
            vertex_id -= _objects[i]->_V[j].rows();
        }
    }
    return -1;
}

bool LinearBlendSkinningUI::mouse_down(int button, int modifier) {
    if (_stage == DEFORMATION) {
        if (button == GLUT_LEFT_BUTTON) {
            if ((modifier & GLUT_ACTIVE_SHIFT) && !(modifier & GLUT_ACTIVE_CTRL)) {
                if (_selected_object_id != -1) {
                    _moving_handle_id = get_closest_handle_id(_selected_object_id);
                    draw_handles();
                    if (_moving_handle_id != -1) {
                        Eigen::MatrixXd H;
                        _objects[_selected_object_id]->_handles.transformed_handles(H);
                        _sel_pos = H.row(_moving_handle_id);
                    }
                }
            }
        }
    }
    return false;
}

bool LinearBlendSkinningUI::mouse_up(int button, int modifier)
{
    if (_moving_handle_id != -1) {
        _moving_handle_id = -1;
        draw_handles();
        for (int i = 0;i < _objects.size();i++)
            _objects[i]->finalize_operation();
    }

    if (_stage == DEFORMATION) {
        if (button == GLUT_LEFT_BUTTON) {
            if (modifier & GLUT_ACTIVE_CTRL) {
                if (viewer->data().V.rows() > 0 && viewer->data().F.rows() > 0) {
                    int v = get_closest_mesh_vertex();
                    if (v != -1) {
                        _selected_object_id = query_object_by_vertex(v);
                    } else {
                        _selected_object_id = -1;
                    }
                    update_colors();
                    update_viewer(false);
                    return true;
                }
            }
        }
    }

    return false;
}

bool LinearBlendSkinningUI::mouse_move(int mouse_x, int mouse_y) {
    if (_selected_object_id != -1 && _moving_handle_id != -1) {
        float x = viewer->current_mouse_x;
        float y = viewer->core().viewport(3) - viewer->current_mouse_y;

        Eigen::RowVector3f orig_pos = _sel_pos.cast<float>();

        Eigen::RowVector3f orig_screen_pos;

        igl::project(
            orig_pos,
            viewer->core().view,
            viewer->core().proj,
            viewer->core().viewport,
            orig_screen_pos
        );

        Eigen::RowVector3f new_screen_pos((float)x, (float)y, orig_screen_pos(2));
        Eigen::RowVector3f new_pos;

        igl::unproject(
            new_screen_pos,
            viewer->core().view,
            viewer->core().proj,
            viewer->core().viewport,
            new_pos
        );

        Eigen::RowVector3d pos = new_pos.cast<double>();

        if (_naive_deformation_mode == 0)
            _objects[_selected_object_id]->move_handle(_moving_handle_id, pos);

        draw_mesh();
        draw_handles();
        return true;
    }
    return false;
}

bool LinearBlendSkinningUI::key_down(int key, int modifiers) {
    if (key == 'c' || key == 'C') {
        if (_selected_object_id != -1) {
            _handle_mode = 1;
            return true;
        }
    } else if (key == 'd' || key == 'D') {
        if (_selected_object_id != -1) {
            _handle_mode = 2;
            return true;
        }
    }
    return false;
}

bool LinearBlendSkinningUI::key_up(int key, int modifiers) {
    _handle_mode = 0;

    return false;
}

void LinearBlendSkinningUI::export_design(std::string name) {
    // export 3d printing mesh
    {
        // create foldre
        std::string folder = PROJECT_SOURCE_DIR "/exported _designs/" + name + "/printing_meshes/";
        if (!std::experimental::filesystem::exists(folder.c_str())) {
            std::experimental::filesystem::create_directories(folder.c_str());
        }

        std::vector<TriMesh> meshes;
        meshes.clear();

        // export palm
        {
            TriMesh mesh;
            for (int i = 0;i < _objects.size();i++)
                if (_objects[i]->_symbol == "palm") {
                    MatrixX3d transformed_vertices = _objects[i]->transform_parent_mesh(_naive_deformation_mode);
                    mesh.merge(transformed_vertices.transpose(), _objects[i]->_parent_mesh._F);
                }

            // add mount
            LinearBlendSkinningObject object;
            
            object._handles = _palm_grammar._mount_handles_template;
            
            std::vector<Eigen::Matrix3Xd> V;
            std::vector<Eigen::Matrix3Xi> F;
            V.push_back(_palm_grammar._symbol_visual_meshes['m']._V);
            F.push_back(_palm_grammar._symbol_visual_meshes['m']._F);
            object.set_mesh(V, F);

            object.compute_weights();

            // find W id
            int w_id = _palm_grammar._cell_index(0, _palm_grammar._num_cols / 2);

            for (int i = 0;i < _palm_grammar._mount_connection_pairs.size();i++) {
                Vector3 pos = _objects[w_id]->_handles.get_handle_position(_palm_grammar._mount_connection_pairs[i].first);
                object._handles.move_handle(_palm_grammar._mount_connection_pairs[i].second, pos);
            }

            std::vector<MatrixX3d> transformed_vertices = object.transform_mesh(1);
            TriMesh mount_mesh(transformed_vertices[0].transpose(), object._F[0].transpose());
            mesh.merge(mount_mesh);

            meshes.push_back(mesh);
        }

        // export fingers
        for (int i = 0;i < _finger_grammar._graph._root->_children.size();i++) {
            int object_id = _finger_grammar._graph._root->_children[i]->_id_in_objects;
            _objects[object_id]->export_printing_meshes(meshes);
        }

        std::cerr << "export mesh size = " << meshes.size() << std::endl;
        for (int i = 0;i < meshes.size();i++)
            meshes[i].export_obj(folder + "/" + to_string(i) + ".obj");
    }

    // export tactile quad mesh
    {
        // create foldre
        std::string folder = PROJECT_SOURCE_DIR "/exported _designs/" + name + "/tactile_mesh/";
        if (!std::experimental::filesystem::exists(folder.c_str())) {
            std::experimental::filesystem::create_directories(folder.c_str());
        }
        std::vector<QuadMesh> meshes;
        meshes.clear();
        for (int i = 0;i < _finger_grammar._graph._root->_children.size();i++) {
            int object_id = _finger_grammar._graph._root->_children[i]->_id_in_objects;
            _objects[object_id]->export_tactile_mesh(meshes, std::vector<Vector3>(), std::vector<Vector4i>());
        }
        for (int i = 0;i < meshes.size();i++) {
            meshes[i].export_obj(folder + "/" + to_string(i) + ".obj");
        }
    }

}