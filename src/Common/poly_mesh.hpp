// Tao Du
// taodu@csail.mit.edu
// Nov 11, 2016
#pragma once
// Define a base mesh class, which consists of:
// - vertex_: a m x n matrix where each column is a m-dimensional point.
// - element_: a p x q matrix where each column represents an element in the
// mesh.
// - edge_in_elemnet_: a 2 x r matrix which defines the edges in an element.
//
// Example:
// - A 2-D square represented by a triangle mesh:
// vertex_ = [0 0; 0 1; 1 0; 1 1]';
// element_ = [0 1 2; 1 2 3]';
// edge_in_element_ = [0 1; 0 2; 1 2]';
//
// - A 2-D square represented by a quad mesh:
// vertex_ = [0 0; 0 1; 1 0; 1 1]';
// element_ = [0 1 3 2]';
// edge_in_element_ = [0 1; 1 2; 2 3; 3 0]';
#include <assert.h>
#include <fstream>
#include <sys/stat.h>
#include "Eigen/Dense"
#include <math.h>
#include <float.h>
#include "Common/typedefs.hpp"
#ifdef WIN32
#include <direct.h>
#endif

namespace materials {

    template<int vertex_dim, typename T>
    class PolyMesh {
    protected:
        // Name convention: replace the number with 'Dim' in Matrix2Xd, Vector2d,
        // Vector2i.
        typedef Eigen::Matrix<T, vertex_dim, Eigen::Dynamic> MatrixDimXT;
        typedef Eigen::Matrix<T, vertex_dim, 1> VectorDimT;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;



    public:
        PolyMesh(const MatrixDimXT& vertex, const Eigen::MatrixXi& element,
                 const Eigen::Matrix2Xi& edge_in_element)
                : vertex_(vertex), element_(element), edge_in_element_(edge_in_element) {}
        explicit PolyMesh(const std::string& file_name) {
            ReadFromFile(file_name);
        }
        PolyMesh(const PolyMesh<vertex_dim, T>& poly_mesh)
                : vertex_(poly_mesh.vertex_), element_(poly_mesh.element_),
                  edge_in_element_(poly_mesh.edge_in_element_) {}
        PolyMesh<vertex_dim, T>& operator=(const PolyMesh& poly_mesh) {
            vertex_ = poly_mesh.vertex_;
            element_ = poly_mesh.element_;
            edge_in_element_ = poly_mesh.edge_in_element_;
            return *this;
        }
        virtual ~PolyMesh() {}

        const MatrixDimXT& vertex() const { return vertex_; }
        const VectorDimT vertex(const int index) const {
            assert(index >= 0 && index < static_cast<int>(vertex_.cols()));
            return vertex_.col(index);
        }
        const Eigen::MatrixXi& element() const { return element_; }
        const Eigen::VectorXi element(const int index) const {
            assert(index >= 0 && index < static_cast<int>(element_.cols()));
            return element_.col(index);
        }
        const Eigen::Matrix2Xi& edge_in_element() const {
            return edge_in_element_;
        }
        const MatrixDimXT vertex_in_element(const int index) const {
            assert(index >= 0 && index < static_cast<int>(element_.cols()));
            const int e_dim = static_cast<int>(element_.rows());
            MatrixDimXT v = MatrixXT::Zero(vertex_dim, e_dim);
            for (int i = 0; i < e_dim; ++i) {
                v.col(i) = vertex_.col(element_(i, index));
            }
            return v;
        }
        const int NumOfVertex() const { return static_cast<int>(vertex_.cols()); }
        const int NumOfElement() const { return static_cast<int>(element_.cols()); }


        void GetBoundingBox(std::pair<VectorDimT, VectorDimT>& bounding_box) const {


            //loop over all vertices
            //find min and max in each direction
            for (int i = 0; i < vertex_dim; ++i) {
                T min_dim = DBL_MAX;
                T max_dim = -DBL_MAX;

                for (int j = 0; j < NumOfVertex(); j++) {
                    T val = vertex_(i, j);
                    if (val < min_dim) {
                        min_dim = val;
                    }
                    if (val > max_dim) {
                        max_dim = val;
                    }
                    bounding_box.first(i) = min_dim;
                    bounding_box.second(i) = max_dim;
                }
            }
        }

        void GetScale(VectorDimT& scale) const {
            std::pair<VectorDimT, VectorDimT> bounding_box;
            GetBoundingBox(bounding_box);

            scale = bounding_box.second - bounding_box.first;
        }

        virtual void WriteToFile(const std::string& file_name) const {
            std::size_t found = file_name.rfind("/");
            if (found != std::string::npos) {
                const std::string folder_name = file_name.substr(0, found + 1);
                size_t pos = 0;
                do {
                    pos = folder_name.find_first_of('/', pos + 1);
#ifdef WIN32
					_mkdir(folder_name.substr(0, pos).c_str());
#else
					mkdir(folder_name.substr(0, pos).c_str(), S_IRWXU);
#endif
                } while (pos != std::string::npos);
            }
            std::ofstream fout;
            fout.open(file_name);
            const int vertex_num = static_cast<int>(vertex_.cols());
            const int element_num = static_cast<int>(element_.cols());
            const int element_dim = static_cast<int>(element_.rows());
            const int edge_in_element_dim = 2;
            const int edge_in_element_num = static_cast<int>(edge_in_element_.cols());
            fout << vertex_dim << " " << vertex_num << std::endl
                 << element_dim << " " << element_num << std::endl
                 << edge_in_element_dim << " " << edge_in_element_num << std::endl;
            for (int i = 0; i < vertex_num; ++i) {
                for (int j = 0; j < vertex_dim; ++j) {
                    fout << vertex_(j, i) << " ";
                }
                fout << std::endl;
            }
            for (int i = 0; i < element_num; ++i) {
                for (int j = 0; j < element_dim; ++j) {
                    fout << element_(j, i) << " ";
                }
                fout << std::endl;
            }
            for (int i = 0; i < edge_in_element_num; ++i) {
                for (int j = 0; j < edge_in_element_dim; ++j) {
                    fout << edge_in_element_(j, i) << " ";
                }
                fout << std::endl;
            }
            fout.close();
        }

    protected:
        virtual void ReadFromFile(const std::string& file_name) {
            // file_name should end with '.plm', which stands for 'poly mesh'.
            // The default file format for any mesh:
            // vertex_dim vertex_num
            // element_dim element_num
            // edge_in_element_dim edge_in_element_num
            // v0(0) v0(1) ... v0(vertex_dim - 1)
            // v1(0) v1(1) ... v1(vertex_dim - 1)
            // ...
            // e0(0) e0(1) ... e0(edge_dim - 1)
            // e1(0) e1(1) ... e1(edge_dim - 1)
            // ...
            // ee0(0) ee0(1)
            // ee1(0) ee1(1)
            // ...
            std::ifstream fin;
            fin.open(file_name);
            int v_dim = 0, vertex_num = 0;
            int e_dim = 0, element_num = 0;
            int edge_in_element_dim = 0, edge_in_element_num = 0;
            fin >> v_dim >> vertex_num >> e_dim >> element_num
                >> edge_in_element_dim >> edge_in_element_num;
            assert(v_dim == vertex_dim);
            assert(edge_in_element_dim == 2);
            vertex_ = MatrixXT::Zero(vertex_dim, vertex_num);
            for (int i = 0; i < vertex_num; ++i) {
                for (int j = 0; j < vertex_dim; ++j) {
                    fin >> vertex_(j, i);
                }
            }
            element_ = Eigen::MatrixXi::Zero(e_dim, element_num);
            for (int i = 0; i < element_num; ++i) {
                for (int j = 0; j < e_dim; ++j) {
                    fin >> element_(j, i);
                }
            }
            edge_in_element_ =
                    Eigen::MatrixXi::Zero(edge_in_element_dim, edge_in_element_num);
            for (int i = 0; i < edge_in_element_num; ++i) {
                for (int j = 0; j < edge_in_element_dim; ++j) {
                    fin >> edge_in_element_(j, i);
                }
            }
            fin.close();
        }

        MatrixDimXT vertex_;
        Eigen::MatrixXi element_;
        Eigen::Matrix2Xi edge_in_element_;

    };




}
