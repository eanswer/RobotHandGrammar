#pragma once

#include <igl/boundary_conditions.h>
#include <Eigen/Core>
#include <fstream>
#include <igl/project.h>
#include <igl/lbs_matrix.h>
#include <igl/unproject_ray.h>
#include <limits>

struct WeightHandles {
public:
	WeightHandles() :
		handle_positions(0, 3),
		point_handle_indices(0),
		bone_edge_indices(0, 2),
		cage_edge_points(0, 2),
        cage_tri_mesh_elements(0, 3),
        _stale(true)
	{
		// Nothing to do
	}

	Eigen::MatrixXd& positions() {
		return handle_positions;
	}

    Eigen::MatrixXi& cage_tri_mesh() {
        return cage_tri_mesh_elements;
    }

	void load_handle_file(std::string filename)
	{
		std::ifstream file(filename);
		read_mat(handle_positions, file);
		read_vec(point_handle_indices, file);
		read_mat(bone_edge_indices, file);
		read_mat(cage_edge_points, file);
        read_cage_faces(file);
		file.close();

		handle_transforms.resize(4 * handle_positions.rows(), 3);
		reset_handle_transformations();
        _handle_translation.resize(handle_positions.rows(), 3);
        reset_handle_translation();
		int n = handle_positions.rows();
		igl::lbs_matrix(handle_positions, Eigen::MatrixXd::Identity(n, n), handle_lbs_matrix);
        _stale = true;
	}

	void save_handle_file(std::string filename)
	{
		std::ofstream file(filename);
		write_mat(handle_positions, file);
		write_vec(point_handle_indices, file);
		write_mat(bone_edge_indices, file);
		write_mat(cage_edge_points, file);
        write_cage_faces(file);
		file.close();
	}

	
	void add_point_handle(const Eigen::Vector3d point)
	{
		int n = handle_positions.rows();
		int m = point_handle_indices.rows();
		handle_positions.conservativeResize(n + 1, Eigen::NoChange);
		handle_positions.row(n) = point;
		point_handle_indices.conservativeResize(m + 1, Eigen::NoChange);
		point_handle_indices(m) = n;

		handle_transforms.conservativeResize(4 * handle_positions.rows(), 3);
		handle_transforms.block(4 * n, 0, 4, 3) <<
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
			0, 0, 0;
        _handle_translation.conservativeResize(handle_positions.rows(), 3);
        _handle_translation.row(n).setZero();

		n = handle_positions.rows();
		igl::lbs_matrix(handle_positions, Eigen::MatrixXd::Identity(n, n), handle_lbs_matrix);
        _stale = true;
	}

    void remove_point_handle(const int handle_id) {
        // int n = handle_positions.rows();
        // int m = point_handle_indices.rows();
        // Eigen::MatrixXd new_handle_positions(n - 1, 3);
    }

	void add_bone_edge(const Eigen::Vector3d& tip, const Eigen::Vector3d& tail)
	{
		int n = handle_positions.rows();
		int m = bone_edge_indices.rows();
		handle_positions.conservativeResize(n + 2, Eigen::NoChange);
		handle_positions.row(n) = tip;
		handle_positions.row(n + 1) = tail;
		bone_edge_indices.conservativeResize(m + 1, Eigen::NoChange);
		bone_edge_indices(m, 0) = n;
		bone_edge_indices(m, 1) = n + 1;

		handle_transforms.conservativeResize(4 * handle_positions.rows(), 3);

		handle_transforms.block(4 * n, 0, 4, 3) <<
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
			0, 0, 0;

		handle_transforms.block(4 * (n+1), 0, 4, 3) <<
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
			0, 0, 0;
        
        // TODO: update handle_translation

		n = handle_positions.rows();
		igl::lbs_matrix(handle_positions, Eigen::MatrixXd::Identity(n, n), handle_lbs_matrix);
        _stale = true;
	}

	void add_cage_edge(const Eigen::Vector3d& tip, const Eigen::Vector3d& tail)
	{
		int n = point_handle_indices.rows();
		int m = cage_edge_points.rows();
		add_point_handle(tip);
		add_point_handle(tail);
		cage_edge_points.conservativeResize(m + 1, Eigen::NoChange);
		cage_edge_points(m, 0) = n;
		cage_edge_points(m, 1) = n + 1;
        _stale = true;
	}

	bool boundary_conditions(
		const Eigen::MatrixXd& V, 
		const Eigen::MatrixXi& F, 
		Eigen::VectorXi& b, 
		Eigen::MatrixXd& bc) 
	{
		return igl::boundary_conditions(
			V,
			F,
			handle_positions,
			point_handle_indices,
			bone_edge_indices,
			cage_edge_points,
			b,
			bc
		);
	}

	// Get the Nearest Handle within
	int select_handle(
		const Eigen::Vector2f& pos,
		const Eigen::Matrix4f& model,
		const Eigen::Matrix4f& proj,
		const Eigen::Vector4f& viewport,
		double threshold = std::numeric_limits<double>::max(),
        const unsigned int naive_deformation_mode = 1) 
	{
		
		Eigen::Vector3f s, dir; 
		Eigen::Vector3d src, dst;
		
		igl::unproject_ray(pos, model, proj, viewport, s, dir);
		
		src = s.cast<double>();
		dst = s.cast<double>() + dir.cast<double>();
        Eigen::MatrixXd H;
        transformed_handles(H);
		int closest = -1;
		float best = threshold * threshold;
        if (naive_deformation_mode == 0) {
            for (int i = 0; i < H.rows(); ++i) {
                double t, sqrd;
                igl::project_to_line(
                    H(i, 0), H(i, 1), H(i, 2),
                    src(0), src(1), src(2),
                    dst(0), dst(1), dst(2),
                    t, sqrd
                );
                if (sqrd < best) {
                    best = sqrd;
                    closest = i;
                }
            }
        } else {
            for (int i = 0; i < 8; ++i) {
                double t, sqrd;
                igl::project_to_line(
                    H(i, 0), H(i, 1), H(i, 2),
                    src(0), src(1), src(2),
                    dst(0), dst(1), dst(2),
                    t, sqrd
                );
                if (sqrd < best) {
                    best = sqrd;
                    closest = i;
                }
            }
        }

		return closest;
	}

    void set_handle_positions(Eigen::MatrixXd handle_positions) {
        this->handle_positions = handle_positions;
        reset_handle_transformations();
        reset_handle_translation();
        int n = this->handle_positions.rows();
		igl::lbs_matrix(this->handle_positions, Eigen::MatrixXd::Identity(n, n), this->handle_lbs_matrix);
        _stale = true;
    }

	void reset_handle_transformations() {
		handle_transforms.resize(4 * handle_positions.rows(), 3);
		for (int i = 0; i < handle_positions.rows(); ++i) {
			handle_transforms.block(4*i, 0, 4, 3) <<
				1, 0, 0,
				0, 1, 0,
				0, 0, 1,
				0, 0, 0;
		}
        _stale = true;
	}

    void reset_handle_translation() {
        _handle_translation.resize(handle_positions.rows(), 3);
        _handle_translation.setZero();
    }

	void move_handle(
		int handle_id,
		const Eigen::RowVector3d pos)
	{
		Eigen::RowVector3d old_pos = handle_positions.row(handle_id);
		Eigen::RowVector3d translation = pos - old_pos;
		handle_transforms.block(4 * handle_id + 3, 0, 1, 3) = translation;
        _stale = true;
	}

    void translate_handle(
        int handle_id,
        const Eigen::RowVector3d translation)
    {
        handle_transforms.block(4 * handle_id + 3, 0, 1, 3) = translation;
        _stale = true;
    }

    Eigen::Vector3d get_handle_position(int handle_id) {
        Eigen::Vector3d pos = (handle_positions.row(handle_id) + handle_transforms.block(4 * handle_id + 3, 0, 1, 3)).transpose();
        return pos;
    }

    void finalize_operation() {
        for (int i = 0;i < handle_positions.rows();i++)
            _handle_translation.block(i, 0, 1, 3) = handle_transforms.block(4 * i + 3, 0, 1, 3);
    }

	Eigen::MatrixXd& transform() {
        _stale = false;
		return handle_transforms;
	}

	void transformed_handles(Eigen::MatrixXd& H) {
		H = handle_lbs_matrix * handle_transforms;
	}

	// void visualize_handles(
    //     Eigen::MatrixXd& V,
    //     Eigen::MatrixXi& F,
	// 	Eigen::MatrixXd& points, 
	// 	Eigen::MatrixXd& point_colors, 
	// 	Eigen::MatrixXi& lines, 
	// 	Eigen::MatrixXd& line_colors,
	// 	bool transform = false)
	// {   
    //     int num_points = points.rows();
    //     int num_lines = lines.rows();
        
	// 	Eigen::Vector3d point_handle_color(1.0, 0.7, 0.3); // Red
	// 	Eigen::Vector3d bone_handle_color(0.0, 0.0, 1.0); // Green
	// 	Eigen::Vector3d cage_handle_color(0.0, 1.0, 1.0); // Blue
        
    //     points.conservativeResize(num_points + handle_positions.rows(), 3);
    //     lines.conservativeResize(num_lines + bone_edge_indices.rows() + cage_edge_points.rows() + cage_tri_mesh_elements.rows() * 3, 2);

	// 	if (transform) {
	// 		points.block(num_points, 0, handle_positions.rows(), 3) = handle_lbs_matrix * handle_transforms;
	// 	}
	// 	else {
	// 		points.block(num_points, 0, handle_positions.rows(), 3) = handle_positions;
	// 	}
        
	// 	point_colors.conservativeResizeLike(points);
	// 	lines.block(num_lines, 0, bone_edge_indices.rows(), 2) = bone_edge_indices;
	// 	line_colors.conservativeResize(lines.rows(), 3);

	// 	for (int i = 0; i < point_handle_indices.size(); ++i) {
	// 		point_colors.row(num_points + point_handle_indices(i)) = point_handle_color;
	// 	}
	// 	for (int i = 0; i < bone_edge_indices.rows(); ++i) {
	// 		point_colors.row(num_points + bone_edge_indices(i, 0)) = bone_handle_color;
	// 		point_colors.row(num_points + bone_edge_indices(i, 1)) = bone_handle_color;
	// 		line_colors.row(num_points + i) = bone_handle_color;
	// 	}
	// 	for (int i = 0; i < cage_edge_points.rows(); ++i) {
    //         int u = point_handle_indices(cage_edge_points(i, 0));
    //         int v = point_handle_indices(cage_edge_points(i, 1));
	// 		point_colors.row(num_points + u) = cage_handle_color;
	// 		point_colors.row(num_points + v) = cage_handle_color;
	// 		int line_row = num_lines + bone_edge_indices.rows() + i;
	// 		lines.row(line_row) = Eigen::Vector2i(u, v);
	// 		line_colors.row(line_row) = cage_handle_color;
	// 	}
    //     for (int i = 0;i < cage_tri_mesh_elements.rows();i++) {
    //         for (int j = 0;j < 3;j++) {
    //             int u = cage_tri_mesh_elements(i, j);
    //             int v = cage_tri_mesh_elements(i, (j + 1) % 3);
    //             point_colors.row(num_points + u) = cage_handle_color;
    //             point_colors.row(num_points + v) = cage_handle_color;
    //             int line_row = num_lines + bone_edge_indices.rows() + cage_edge_points.rows() + i * 3 + j;
    //             lines.row(line_row) = Eigen::Vector2i(u, v);
    //             line_colors.row(line_row) = cage_handle_color;
    //         }
    //     }
	// }

    Eigen::MatrixXd handle_transforms;
	Eigen::MatrixXd handle_lbs_matrix;

	Eigen::MatrixXd handle_positions;
	Eigen::VectorXi point_handle_indices; // indices into handle_positions
	Eigen::MatrixXi bone_edge_indices;    // indices into handle_positions
	Eigen::MatrixXi cage_edge_points;     // indices into point_handle_indices
    Eigen::MatrixXi cage_tri_mesh_elements;   // triangle cage mesh for mean value weights

    std::vector<std::vector<int> > cage_faces;

    Eigen::MatrixXd _handle_translation; // the real handle positions before current movement.

    bool _stale;

private:

	template <typename T>
	void read_mat(Eigen::Matrix<T, -1, -1>& mat, std::ifstream& file) {
		int n, m;
		T v;
		file >> n;
		file >> m;
		mat.resize(n, m);
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				file >> v;
				mat(i, j) = v;
			}
		}
	}

	template <typename T>
	void write_mat(const Eigen::Matrix<T, -1, -1> & mat, std::ofstream& file) {
		file << mat.rows() << " ";
		file << mat.cols() << std::endl;
		for (int i = 0; i < mat.rows(); ++i) {
			for (int j = 0; j < mat.cols(); ++j) {
				file << mat(i, j) << " ";
			}
            file << std::endl;
		}
	}

	void read_vec(Eigen::VectorXi& vec, std::ifstream& file) {
		int n;
		int v;
		file >> n;
		vec.resize(n);
		for (int i = 0; i < vec.size(); ++i) {
			file >> v;
			vec(i) = v;
		}
	}

	void write_vec(const Eigen::VectorXi& vec, std::ofstream& file) {
		file << vec.size() << std::endl;
		for (int i = 0; i < vec.size(); ++i) {
			file << vec(i) << " ";
		}
        file << std::endl;
	}

    void read_cage_faces(std::ifstream& file) {
        if (file.eof()) {
            return;
        }

        cage_faces.clear();
        int n;
        file >> n;
        for (int i = 0;i < n;i++) {
            cage_faces.push_back(std::vector<int>());
            int m;
            file >> m;
            for (int j = 0;j < m;j++) {
                int vid;
                file >> vid;
                cage_faces[i].push_back(vid);
            }
        }

        build_cage_tri_mesh();
    }

    void write_cage_faces(std::ofstream& file) {
        file << cage_faces.size() << std::endl;
        for (int i = 0;i < cage_faces.size();i++) {
            file << cage_faces[i].size();
            for (int j = 0;j < cage_faces[i].size();j++) 
                file << " " << cage_faces[i][j];
            file << std::endl;
        }
    }

    void build_cage_tri_mesh() {
        // assume each cage face is convex, thus the centroid is inside the face
        int n = handle_positions.rows();
        handle_positions.conservativeResize(n + cage_faces.size(), 3);

        std::vector<Eigen::Vector3i> ele;
        ele.clear();
        for (int i = 0;i < cage_faces.size();i++) {
            Eigen::Vector3d centroid;
            centroid.setZero();
            for (int j = 0;j < cage_faces[i].size();j++) {
                centroid += handle_positions.row(cage_faces[i][j]);
                ele.push_back(Eigen::Vector3i(cage_faces[i][j], cage_faces[i][(j + 1) % cage_faces[i].size()], n + i));
            }
            centroid /= cage_faces[i].size();
            handle_positions.row(n + i) = centroid; 
        }
        cage_tri_mesh_elements = Eigen::MatrixXi(ele.size(), 3);
        for (int i = 0;i < ele.size();i++)
            cage_tri_mesh_elements.row(i) = ele[i];
    }
};