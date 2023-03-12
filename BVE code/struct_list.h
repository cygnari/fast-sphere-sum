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
};

struct interp_struct {
    int degree;
    int point_count;
    int info;
    vector<double> matrix;
    vector<vector<double>> interp_points;
    vector<int> ipiv;
};

struct run_config {
    bool use_mpi = false;
    bool use_amr = false;
    bool use_remesh = false;
    bool use_fast = false;
    int mpi_ranks = 1;
    int point_count;
    double radius;
    string init_cond;
    icos_struct icos;
    interp_struct interp;
    //
};

#endif
