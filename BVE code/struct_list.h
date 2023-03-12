#ifndef struct_H
#define struct_H

struct interaction_pair {
    int lev_target; //icosahedron level of target/source
    int lev_source;
    int curr_target; // index of target/source
    int curr_source;
    int count_target;
    int count_source;
    string type; // pp, pc, cp, or cc
};

struct icos_struct {
    int levels;
    double radius;
    vector<vector<double>> verts; // list of icosahedron vertices
    vector<vector<vector<double>>> tri_info; // information about each triangle forming the faces
    vector<vector<vector<int>>> tri_verts; // length 3 int vector listing the vertex indices
    vector<vector<vector<int>>> tri_points; // points inside each triangle
    vector<interaction_pair> interactions; // interactions
    vector<vector<int>> point_locs; // triangle each point is in
};

struct interp_struct {
    int degree;
    int point_count;
    int info;
    vector<double> matrix;
    vector<vector<double>> interp_points;
    vector<int> ipiv;
};

struct mpi_struct {
    int mpi_ranks; // total processors
    int rank; // current processor ID
    int particle_lb; // even division of particles
    int particle_ub;
    int interact_lb; // even division of fast interactions
    int interact_ub;
};

struct run_config {
    bool use_mpi = false;
    bool use_amr = false;
    bool use_remesh = false;
    bool use_fast = false;
    bool testing = false;
    int mpi_ranks = 1;
    int point_count;
    int many_count;
    int levels = 1;
    double radius;
    double end_time;
    double theta;
    string init_cond;
    icos_struct icos;
    interp_struct interp;
    //
};

#endif
