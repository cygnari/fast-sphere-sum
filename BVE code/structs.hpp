#ifndef struct_H
#define struct_H

using namespace std;

struct run_config {
    // primitive configuration options
    // bool use_mpi = false;
    bool use_amr = false;
    bool use_remesh = false;
    bool use_fast = false;
    bool use_fixer = false;
    bool write_output = false;
    bool write_tris = false;
    double radius = 1.0;
    double end_time;
    double delta_t; // time step size
    int dynamics_levels_min; // initial icosahedron refinement levels, 0 = base icosahedron
    int dynamics_levels_max; // max icosahedron refinement if using amr
    string initial_vor_condition; // initial vorticity distribution
    int fast_sum_cluster_thresh; // threshold for a triangle being a cluster
    int fast_sum_tree_levels; // icosahedron levels for fast summation
    double fast_sum_theta; // well separated threshold
    int interp_degree; // interpolation degree
    int interp_point_count; // number of interpolation points
    int info_per_point; // how many doubles each point is, for example, storing x y z vor tracer = 5
    double amr_circ_thresh = 0.005; // threshold for circulation in amr
    double amr_vor_thresh = 0.4; // threshold for vorticity difference in amr

    // derived run config info
    int time_steps; // number of time steps
    int dynamics_initial_points; // initial number of dynamics points
    int dynamics_max_points; // maximum number of dynamics points
    int dynamics_initial_triangles; // initial number of dynamics triangles
    int dynamics_max_triangles; // maximum number of dynamics triangles
    int dynamics_curr_point_count; // current number of points
    int dynamics_curr_tri_count; // current number of triangles
    int tracer_count; // number of tracers

    // mpi info
    int particle_lb;
    int particle_ub;
    int interaction_lb;
    int interaction_ub;
};

struct interaction_pair {
    int lev_target; //icosahedron level of target/source
    int lev_source;
    int curr_target; // index of target/source
    int curr_source;
    int count_target; // number of particles in target/source
    int count_source;
    string type; // pp, pc, cp, or cc
};

#endif
