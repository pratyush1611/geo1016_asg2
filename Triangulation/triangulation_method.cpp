/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <iostream>

using namespace easy3d;


/// convert a 3 by 3 matrix of type 'Matrix<double>' to mat3
mat3 to_mat3(Matrix<double> &M) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}


/// convert M of type 'matN' (N can be any positive integer) to type 'Matrix<double>'
template<typename mat>
Matrix<double> to_Matrix(const mat &M) {
    const int num_rows = M.num_rows();
    const int num_cols = M.num_columns();
    Matrix<double> result(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}

Matrix<double> norm_transform(std::vector<vec3> pts)
{
    Matrix<double> T (3,3, 0.0);
    //1st calculate the centroid

    float avgx = 0.0;
    float avgy = 0.0;
    for(auto & i : pts)
    {
        avgx += i.x;
        avgy += i.y;
    }

    avgx/=pts.size() ;
    avgy/=pts.size() ;

    vec3 centrd {avgx, avgy, 1.0};

    //translate the coordinates
    //get mean distances
//    std::vector<vec3> translated;
    double mean_dist=0.0;
    for(auto & j : pts)
    {
        vec3 temp = j - centrd;
        mean_dist += temp.length();
//        translated.push_back(temp);
    }
    mean_dist/=pts.size();
    mean_dist = sqrt(2.0)/mean_dist;

    T.set_row({mean_dist, 0.0, 0.0},0);
    T.set_row({0.0, mean_dist, 0.0},1);
    T.set_row({-1.0*centrd.x, -1.0*centrd.y, 1.0},2);

    return T.transpose();
}
/*
std::vector<vec3> centroidinator(std::vector<vec3> points_vec, std::vector<vec3> &translated)
{
    //1st calculate the centroid
    vec3 centrd;
    double avgx,avgy = 0.0, 0.0 ;
    for(auto & i : points_vec)
    {
        avgx += i.x;
        avgy += i.y;
    }

    avgx/=points_vec.size() ;
    avgy/=points_vec.size() ;

    centrd= {avgx, avgy, 1};

    //run a translation
    for(auto & j : points_vec)
    {
        vec3 temp = j - centrd;
        translated.push_back(temp);
    }

    return translated;
}
*/
/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'.
 */
bool Triangulation::triangulation(
        float fx, float fy,     /// input: the focal lengths (same for both cameras)
        float cx, float cy,     /// input: the principal point (same for both cameras)
        const std::vector<vec3> &points_0,    /// input: image points (in homogenous coordinates) in the 1st image.
        const std::vector<vec3> &points_1,    /// input: image points (in homogenous coordinates) in the 2nd image.
        std::vector<vec3> &points_3d,         /// output: reconstructed 3D points
        mat3 &R,   /// output: recovered rotation of 2nd camera (used for updating the viewer and visual inspection)
        vec3 &t    /// output: recovered translation of 2nd camera (used for updating the viewer and visual inspection)
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for each sub-task. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or feel free to put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any data structure and function offered by Easy3D, in particular the following two\n"
                 "\tfiles for vectors and matrices:\n"
                 "\t    - easy3d/core/mat.h  Fixed-size matrices and related functions.\n"
                 "\t    - easy3d/core/vec.h  Fixed-size vectors and related functions.\n"
                 "\tFor matrices with unknown sizes (e.g., when handling an unknown number of corresponding points\n"
                 "\tstored in a file, where their sizes can only be known at run time), a dynamic-sized matrix data\n"
                 "\tstructure is necessary. In this case, you can use the templated 'Matrix' class defined in\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\tPlease refer to the corresponding header files for more details of these data structures.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please refer to\n"
                 "\t'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations. \n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (original code framework + your implementation).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without any modification.\n\n" ;//<< std::flush;

    /// Easy3D provides fixed-size matrix types, e.g., mat2 (2x2), mat3 (3x3), mat4 (4x4), mat34 (3x4).
    /// To use these matrices, their sizes should be known to you at the compile-time (i.e., when compiling your code).
    /// Once defined, their sizes can NOT be changed.
    /// In 'Triangulation/matrix.h', another templated 'Matrix' type is also provided. This type can have arbitrary
    /// dimensions and their sizes can be specified at run-time (i.e., when executing your program).
    /// Below are a few examples showing some of these data structures and related APIs.


    // TODO: delete all above demo code in the final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    std::cout<<"checking print statement\n";
//    if(cx != cy) return false;
//    if(fx != fy) return false;
    if(points_0.size() <8 || points_1.size()<8 || points_0.size() != points_1.size() ) return false;
    std::cout<<"if doesnt return \n";

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //comput W matrix
    // w has a size of Nx9; N = number of points
    int no_pt = points_0.size();
    Matrix<double> W(no_pt, 9, 0.0); // all entries initialized to 0.0.
    // iterate through the vector of points and craete rows of W on the go

    // to normalize get T and T'
    Matrix<double> T = norm_transform(points_0);
    Matrix<double> T_prime = norm_transform(points_1);

    std::cout<<"\nT is \n"<<T<<'\n';
    std::cout<<"\nT' is \n"<<T_prime<<'\n';

/*
    std::cout<<T+ Matrix<double> (3,3,2.2);
    std::cout<<T* Matrix<double> (3,3,2.2);
    Matrix<double> Z(3,1,0.0);
    Z.set_column( {10,10,10},0);
    std::cout<<Z<<'\n'<< points_1[0].data();
    std::cout<<T* Z;
  */
    //tranform coordinates
    std::vector<vec3> q, q_prime;
    for(auto &i: points_0)
    {
        Matrix<double> i_mat(3,1, 0.0);

        i_mat.set_column( {i[0],i[1],i[2]} , 0);
//        std::cout<<T * i_mat;
        auto temp_ = T * i_mat
        q.push_back({temp_[0][0] , temp_[1][0], temp_[2][0]}  );
//        std::cout<<i_mat;
    }
//    for(auto &i: points_1) q_prime.push_back(T_prime * i);
    std::cout<<"q is \n"<<q;

    //we find Fq
    for(int i =0 ; i<no_pt; i++)
    {
        std::vector<double> row_to_set = {(points_0[i].x * points_1[i].x), (points_0[i].y * points_1[i].x), (points_1[i].x ),
                                          (points_0[i].x * points_1[i].y), (points_0[i].y * points_1[i].y), (points_1[i].y ),
                                           points_0[i].x,                    points_0[i].y,                  1.0
                                         } ;

        W.set_row(row_to_set, i);
    }

    //find F from T T' and Fq

    // estimate F
    // F is estimated as last column of Vt in SVD of W
    //SVD decomposition of the above W-matrix
    Matrix<double> U(no_pt, no_pt, 0.0);
    Matrix<double> S(no_pt, 9, 0.0);
    Matrix<double> Vt(9, 9, 0.0);
    svd_decompose(W, U, S, Vt);

    // We get m by taking the last column of V
//    std::vector<double> f = Vt.get_column(8);//get last col
    std::vector<double> f = Vt.get_column(Vt.cols() - 1);

    //      - compute the essential matrix E;

    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you to check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       However, there are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
//    return points_3d.size() > 0;
    return true;
}

/// ----------- fixed-size matrices
/*
/// define a 3 by 4 matrix M (you can also define 3 by 4 matrix similarly)
mat34 M(1.0f);  /// entries on the diagonal are initialized to be 1 and others to be 0.

/// set the first row of M
M.set_row(0, vec4(1,1,1,1));    /// vec4 is a 4D vector.

/// set the second column of M
M.set_col(1, vec4(2,2,2,2));

/// get the 3 rows of M
vec4 M1 = M.row(0);
vec4 M2 = M.row(1);
vec4 M3 = M.row(2);

/// ----------- fixed-size vectors

/// how to quickly initialize a std::vector
std::vector<double> rows = {0, 1, 2, 3,
                            4, 5, 6, 7,
                            8, 9, 10, 11};
/// get the '2'-th row of M
const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

/// get the '1'-th column of M
const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

/// modify the element value at row 2 and column 1 (Note the 0-based indices)
M(2, 1) = b.x;

/// apply transformation M on a 3D point p (p is a 3D vector)
vec3 p(222, 444, 333);
vec3 proj = M * vec4(p, 1.0f);  // use the homogenous coordinates. result is a 3D vector

/// the length of a vector
float len = p.length();
/// the squared length of a vector
float sqr_len = p.length2();

/// the dot product of two vectors
float dot_prod = dot(p, proj);

/// the cross product of two vectors
vec3 cross_prod = cross(p, proj);

/// normalize this vector
cross_prod.normalize();

/// a 3 by 3 matrix (all entries are intentionally NOT initialized for efficiency reasons)
mat3 F;
/// ... here you compute or initialize F.
/// compute the inverse of K
mat3 invF = inverse(F);

/// ----------- dynamic-size matrices

/// define a non-fixed size matrix
Matrix<double> W(2, 3, 0.0); // all entries initialized to 0.0.

/// set its first row by a 3D vector (1.1, 2.2, 3.3)
W.set_row({ 1.1, 2.2, 3.3 }, 0);   // here "{ 1.1, 2.2, 3.3 }" is of type 'std::vector<double>'
*/
// get the last column of a matrix
//    std::vector<double> last_column = W.get_column(W.cols() - 1);
