#ifndef mpiu_H
#define mpiu_H

#include <mpi.h>

void sync_state_direct(vector<double>& state, int ID, int P, vector<int>& lb, vector<int>& counts, MPI_Win *win) { // syncs state for direct sum
    MPI_Win_fence(0, *win);
    for (int i = 0; i < P; i++) {
        if (i != ID) {
            MPI_Get(&state[5 * lb[i]], 5 * counts[i], MPI_DOUBLE, i, 5 * lb[i], 5 * counts[i], MPI_DOUBLE, *win);
        }
    }
    MPI_Win_fence(0, *win);
}

void sync_state_fast(vector<double>& state, int point_count, int ID, vector<int>& particle_thread, MPI_Win *win) { // syncs state for fast sum
    MPI_Win_fence(0, *win);
    int thread;
    for (int i = 0; i < point_count; i++) {
        thread = particle_thread[i];
        if (thread != ID) {
            MPI_Get(&state[5 * i], 5, MPI_DOUBLE, thread, 5 * i, 5, MPI_DOUBLE, *win);
        }
    }
    MPI_Win_fence(0, *win);
}

#endif
