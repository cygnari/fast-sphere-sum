#ifndef struct_H
#define struct_H

using namespace std;

struct interaction_pair {
    int lev_target; //icosahedron level of target/source
    int lev_source;
    int curr_target; // index of target/source
    int curr_source;
    int count_target;
    int count_source;
    string type; // pp, pc, cp, or cc
};

// struct icos_struct {
//     int levels;
//     double radius;
//     vector<vector<double>> verts; // list of icosahedron vertices
//     vector<vector<vector<double>>> tri_info; // information about each triangle forming the faces
//     vector<vector<vector<int>>> tri_verts; // length 3 int vector listing the vertex indices
//     vector<vector<vector<int>>> tri_points; // points inside each triangle
//     vector<interaction_pair> interactions; // interactions
//     vector<vector<int>> point_locs; // triangle each point is in
// };

// struct interp_struct {
//     int degree;
//     int point_count;
//     int info;
//     vector<double> matrix;
//     vector<vector<double>> interp_points;
//     vector<int> ipiv;
// };

struct mpi_struct {
    int mpi_ranks; // total processors
    int rank; // current processor ID
    int particle_lb; // even division of particles
    int particle_ub;
    int interact_lb; // even division of fast interactions
    int interact_ub;
};

// struct amr_struct {
//     int max_count;
//     vector<vector<int>> parent_verts;
// };

struct run_information  {
    // run configuration options
    bool use_mpi = false;
    bool use_amr = false;
    bool use_remesh = false;
    bool use_fast = false;
    bool testing = false;
    double radius = 1.0;
    double end_time;

    // dynamics information
    int point_count; // number of points, 4 ^ start_level * 10 + 2
    int tri_count; //number of triangles
    int start_level = 5; // starting number of icosahedron refinements
    string init_cond; // initial vorticity
    vector<double> dynamics_state; // 5 * point_count long, 0 is x coord, 1 is y coord, 2 is z coord, 3 is relative vorticity, 4 is tracer
    vector<vector<int>> dynamics_tris; // vector of length 3 vectors, listing the vertices of each triangle
    vector<double> dynamics_areas; // vector containing the area assigned to each point
    vector<int> dynamic_tri_levels; // refinement that each triangle exists at

    // amr information
    int max_count;
    int max_level;
    vector<vector<int>> amr_parent_verts; // vector of length 2 vectors, listing the two vertices that each point comes from, -1 for primitive point
    vector<vector<int>> amr_vert_tris; // vert_tris[i] is the int vector containing the indices of the triangles adjacent to point i

    // fast sum info
    int many_count; // threshold for cluster
    int tree_levels = 1; // how many levels to use in the tree refinement
    double theta = 0.7; // distance ratio threshold
    vector<vector<double>> tree_verts; // list of icosahedron vertices
    vector<vector<vector<double>>> tree_tri_info; // information about each triangle forming the faces
    vector<vector<vector<int>>> tree_tri_verts; // length 3 int vector listing the vertex indices
    vector<vector<vector<int>>> tree_tri_points; // points inside each triangle
    vector<interaction_pair> tree_interactions; // interactions
    vector<vector<int>> tree_point_locs; // triangle each point is in

    // interpolation information
    int degree; // degree of interpolation
    int info; // 0 is LU decomp worked
    int interp_count; // number of interpolation points
    vector<double> interp_matrix; // vandermone matrix for the interpolation problem
    vector<vector<double>> interp_points; // vector of length 2 vectors of barycentric coordinates of interpolation points
    vector<int> ipiv; // pivots for LU decomp

    // mpi info
    int mpi_ranks = 1; // total processors
    int rank = 0; // current processor ID
    int particle_lb; // even division of particles
    int particle_ub;
    int interact_lb; // even division of fast interactions
    int interact_ub;
};

#endif
