#include "linear_blend_skinning_ui.h"
#include <string>

int main(int argc, char** argv) {
    igl::opengl::glfw::Viewer viewer;
    LinearBlendSkinningUI ui;
    viewer.plugins.push_back(&ui);
    viewer.launch();

	return 0;
}