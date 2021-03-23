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

void norm_transform(mat3 &t , std::vector<vec3> &pts)
{
//    Matrix<double> T (3,3, 0.0);
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
//    make translation matrix
    mat3 T( 1.0);
    T.set_col(2, vec3(-avgx, -avgy, 1.0));

    //translate the coordinates
    //get mean distances

    double mean_dist=0.0;
    for(auto & j : pts)
    {
//        vec3 temp = j - centrd;
//        mean_dist += temp.length();
        mean_dist += j.length() ;
    }
    mean_dist/=pts.size();

    float s = sqrt(2.0)/mean_dist;
    mat3 S = mat3::scale (s, s,1.0);
/*
    T.set_row({mean_dist, 0.0, 0.0},0);
    T.set_row({0.0, mean_dist, 0.0},1);
    T.set_row({-1.0*centrd.x, -1.0*centrd.y, 1.0},2);

    t= to_mat3( T.transpose() ); //returns matrix for scaling and translating

    T.set_row({mean_dist, mean_dist, 1.0},0);
    T.set_row({mean_dist, mean_dist, 1.0},1);
    T.set_row({-1.0*mean_dist*centrd.x, -1.0*mean_dist*centrd.y, 1.0 },2);


    t= to_mat3( T.transpose() ); //returns matrix for scaling and translating
*/
    t = S*T;

    for(auto &i: pts)
    {
        i = t * i;
    }

}


void Essentialator(float fx, float fy,     float cx, float cy,
                   const std::vector<vec3> &points_0, const std::vector<vec3> &points_1,
                   Matrix<double> &E, Matrix<double> &K)
{
    //make copies of points
    std::vector<vec3> q=points_0;
    std::vector<vec3> q_prime = points_1;
    // to normalize get T and T'
    mat3 T, T_prime;
    norm_transform(T, q);
    norm_transform(T_prime, q_prime);

//    std::cout<<"\nT  is \n"<<T<<'\n';
//    std::cout<<"\nT' is \n"<<T_prime<<'\n';

    /*
    //tranform coordinates
    std::vector<vec3> q, q_prime;
    for(auto &i: points_0)
    {
        vec3 q_vect = T * i ;
        q.push_back( q_vect );
    }
    for(auto &i: points_1)
    {
        vec3 q_vec_prime =T_prime * i ;
        q_prime.push_back(q_vec_prime);
    }
    */
    //comput W matrix
    // w has a size of Nx9; N = number of points
    int no_pt = points_0.size();
    Matrix<double> W(no_pt, 9, 0.0); // all entries initialized to 0.0.
    // iterate through the vector of points and craete rows of W on the go
    for(int i =0 ; i<no_pt; i++)
    {
        std::vector<double> row_to_set = {(q[i].x * q_prime[i].x), (q[i].y * q_prime[i].x), (q_prime[i].x ),
                                          (q[i].x * q_prime[i].y), (q[i].y * q_prime[i].y), (q_prime[i].y ),
                                          q[i].x,                    q[i].y,                  1.0
                                         } ;
        W.set_row(row_to_set, i);
    }

// find F_q wrt new coordinates
    // F is estimated as last column of Vt in SVD of W
    //SVD decomposition of the above W-matrix
    Matrix<double> U(no_pt, no_pt, 0.0);
    Matrix<double> S(no_pt, 9, 0.0);
    Matrix<double> Vt(9, 9, 0.0);
    svd_decompose(W, U, S, Vt);
    //estimate F as last col of Vt
//    std::vector<double> f = Vt.get_column(Vt.cols() - 1);
//    std::cout<<"USV" << W<<'\n'<<S<<'\n'<<Vt<<'\n';
    // before rank 2 setting in F, with normalised coordinates
    mat3 F_q_unrank;
    int ind=0;
    for(int i=0; i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            F_q_unrank(i,j) = Vt(ind,8) ;
            ind++;
        }

    }
    /*
    vec3 r1 (f[0],f[1],f[2]);
    vec3 r2 (f[3],f[4],f[5]);
    vec3 r3 (f[6],f[7],f[8]); //last value set to 1.0 because so is asked in assignment

    F_q_unrank.set_row(0, r1);
    F_q_unrank.set_row(1, r2);
    F_q_unrank.set_row(2, r3);
    */

    Matrix<double> Ua (3,3,0.0);
    Matrix<double> Sa (3,3,0.0);
    Matrix<double> Va_T (3,3,0.0);
    Matrix<double> new_Sa (3,3,0.0);
    Matrix<double> F_q_unrank_mx(3,3,0.0) ;

    F_q_unrank_mx = to_Matrix(F_q_unrank);
    svd_decompose(F_q_unrank_mx, Ua, Sa, Va_T);

    std::cout<<"F_q unrank check \n"<< F_q_unrank_mx <<'\n';
    new_Sa=Sa;
    new_Sa(2,2)=0;

    mat3 F_q = to_mat3(Ua * (new_Sa) * Va_T.transpose() );
    std::cout<< "Fq rank check is \n"<<F_q <<'\n';

    //find F from T T' and Fq : denormalize
    mat3 F = to_mat3(to_Matrix(T_prime).transpose()) * F_q * T;
    //rescale wrt F(2,2) aka set last element as 1
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            F(i,j) = F(i,j) / F(2,2);
        }
    }
    std::cout<<"F after denormalisation and normalising last value to 1 is \n"<< F <<'\n';

    //      - compute the essential matrix E;
    // assuming that we can make the intrinsic matrix from the cx cy fx fy
    //make vectors k and k'
    K = Matrix<double> (3,3,
                        std::vector<double> {fx,0.0,cx,
                                                  0.0,fy,cy,
                                                  0.0,0.0,1.0});
    E = K.transpose() * to_Matrix(F) * K;
    std::cout<<"E\n"<<E<<std::endl;
//    std::cout<<"K\n"<<K<<std::endl;

}

Matrix<double> get_KRT( Matrix<double> R, std::vector<double> T, Matrix<double> K )
{
    // create M matrix for given R T and K in the equation AP=0
    // M = K [R | T]
    Matrix<double> RT (3,4,  0.0);
    RT.set_column( R.get_column(0) , 0);
    RT.set_column( R.get_column(1) , 1);
    RT.set_column( R.get_column(2) , 2);

    RT.set_column( T , 3);

    Matrix<double> M = K* RT;
    return M;
}

int rt_decider( Matrix<double> R, std::vector<double> T, Matrix<double> K , std::vector<vec3> pts_0, std::vector<vec3> pts_1)
{
    //vec3 p, vec3 p_prime ,
    /* create M matrix for given R T and K in the equation AP=0
     * then runs SVD to compute for P for a given point
     */
    // pt_0 for m2 [i|o]
    // pt_1 for m1;
    Matrix<double> M1 = get_KRT(R,T,K);
    Matrix<double> M2 (3,4,0.0);
    std::vector<double> I_row1 {0.0, 0.0, 1.0, 0.0};
    std::vector<double> I_row2 {0.0, 0.0, 1.0, 0.0};
    std::vector<double> I_row3 {0.0, 0.0, 1.0, 0.0};
    M2.set_row( I_row1, 0);
    M2.set_row( I_row2, 1);
    M2.set_row( I_row3, 2);
    std::cout<<"m2 set\n";
    //randomly pick points to choose the side of the orientation
    //for now deal with 30 at interval of 5

    Matrix<double> U(3,3, 0.0);
    Matrix<double> S(3,3, 0.0);
    Matrix<double> Vt(3,3, 0.0);
    int s=0;
    for(int i=0; i <= pts_0.size(); i+=5)
    {
        //compute A matrix for a selected point
        Matrix<double> A (4,4,0.0);
        A.set_row( pts_0[i].x * M2.get_row(2) - M2.get_row(0), 0);
        A.set_row( pts_0[i].y * M2.get_row(2) - M2.get_row(1), 1);
        A.set_row( pts_1[i].x * M1.get_row(2) - M1.get_row(0), 2);
        A.set_row( pts_1[i].y * M1.get_row(2) - M1.get_row(1), 3);

        //sun SVD of A, to get P matrix contained in last col of V vector
        svd_decompose(A, U, S, Vt);
        std::vector<double> p = U.get_column(U.cols() - 1);
//        Matrix<double> P (4,1, p);
        //check if p.z is mositive
        std::cout<<"wir sind hier!\n";
        if( p[2] > 0 ) s++;
    }

    return s;

//    std::cout<<"\nR is " <<R<<'\n';
//    std::cout<<"\nT is " <<T<<'\n';
//    std::cout<<"\nRT is " <<RT<<'\n';
//    std::cout<<"K is \n"<<K<<'\n';
}

void points_3d_gen( Matrix<double> R, std::vector<double> T, Matrix<double> K , std::vector<vec3> pts_0, std::vector<vec3> pts_1, std::vector<vec3> & pt3d)
{
    //vec3 p, vec3 p_prime ,
    /* create M matrix for given R T and K in the equation AP=0
     * then runs SVD to compute for P for a given point
     */
    // pt_0 for m2 [i|o]
    // pt_1 for m1;
    Matrix<double> M1 = get_KRT(R,T,K);
    Matrix<double> M2 (3,4,0.0);

    M2.set_column( std::vector<double> {1.0, 0.0, 0.0} , 0);
    M2.set_column( std::vector<double> {0.0, 1.0, 0.0} , 1);
    M2.set_column( std::vector<double> {0.0, 0.0, 1.0} , 2);

    //randomly pick points to choose the side of the orientation
    //for now deal with 30 at interval of 5

    Matrix<double> U(4,4, 0.0);
    Matrix<double> S(4,4, 0.0);
    Matrix<double> Vt(4,4, 0.0);

    int s=0;

    for(int i=0; i<pts_0.size(); i++)
    {
        Matrix<double> A (4,4,0.0);
        A.set_row( pts_0[i].y * M2.get_row(2) - M2.get_row(1), 1);
        A.set_row( pts_1[i].x * M1.get_row(2) - M1.get_row(0), 2);

        A.set_row( pts_0[i].x * M2.get_row(2) - M2.get_row(0), 0);
        A.set_row( pts_1[i].y * M1.get_row(2) - M1.get_row(1), 3);

        // sun SVD of A, to get P matrix contained in last col of V vector
        svd_decompose(A, U, S, Vt);
        std::vector<double> p = Vt.get_column(Vt.cols() - 1);
        vec3 temporary (p[0], p[1], p[2]) ;
        pt3d.push_back( temporary );
    }

}



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
    // TODO: delete all above demo code in the final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    std::cout<<"checking print statement\n";
    if(points_0.size() <8 || points_1.size()<8 || points_0.size() != points_1.size() ) return false;
    std::cout<<"if doesnt return \n";

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //NORMalize, SVD get f, SVD rank 2 check, DEnormalize
    Matrix<double> E, K;
    Essentialator(fx, fy, cx, cy, points_0, points_1, E, K);

    //      - recover rotation R and t.
    //set values for W and Z
    Matrix<double> We ( 3, 3, std::vector<double>{0.0, -1.0, 0.0,
                                                1.0, 0.0, 0.0,
                                                0.0, 0.0, 1.0});
    Matrix<double> Z ( 3, 3, std::vector<double>{0.0, 1.0, 0.0,
                                                 -1.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0});
    //run SVD on E
    Matrix<double> Ue(3,3, 0.0);
    Matrix<double> Se(3,3, 0.0);
    Matrix<double> Vt_e(3,3, 0.0);

    svd_decompose(E, Ue, Se, Vt_e);

    std::vector<double> t_1 = Ue.get_column(Ue.cols() - 1);
    std::vector<double> t_2 = -1.0*Ue.get_column(Ue.cols() - 1);
    std::vector< std::vector<double> > t_types {t_1, t_2};

    Matrix<double> R1 = determinant(Ue * We * Vt_e.transpose()) * Ue * We * Vt_e.transpose();
    Matrix<double> R2 = determinant(Ue * We.transpose() * Vt_e.transpose()) * Ue * We.transpose() * Vt_e.transpose();
    std::vector< Matrix<double> > R_types {R1, R2 };

//    std::cout<<"translation vectors\n"<< t_1<<'\n'<< t_2<<'\n';
    std::cout<<"R1\n"<<R1<<"\nR2\n"<<R2<<"\ndetR1\n"<<determinant(R1)<<"\ndetR2\n"<<determinant(R2)<<'\n';


    //print
//    std::cout<<"Ue, Se, Vte are \n"<<Ue<<'\n'<<Se<<'\n'<<Vt_e<<'\n';
    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    /* do triangulation //AP = 0
     * run SVD on A to find values of P
     * run on at least 30 random points  choose majority with +ve side as correct orientation
     * check rank of positive returned for the 4 combinations
     * once returned, check for
     */
//    for(int i=0; i< 30; i++)
//    {
////        points_0[i];
//        tiangulator( R1, t_1, K );
//    }
    Matrix<double> r_store = R1;
    std::vector<double> t_store = t_1;
    int no_of_positive_points =0;
    for(auto &i : t_types )
    {
        for( auto &j: R_types)
        {
            std::cout<<"in loop again \n";
            int ret = rt_decider(j, i, K, points_0, points_1);
            if( no_of_positive_points < ret )
            {
                no_of_positive_points = ret;
                //also store the r and t
                t_store = i;
                r_store = j;
            }
        }
    }

    //final rstore and tstore should be here
    std::cout<<"the best r and t are \n"<<r_store<<'\n'<<t_store <<'\n';
    points_3d_gen( r_store, t_store, K , points_0, points_1 , points_3d);


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
    return points_3d.size() > 0;
//    return true;
}

void liangliang() {
//for all the stuff liang liang gave us

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
}