#include <iostream>
#include <cmath>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <chrono>
#include <mpi.h>

using namespace std;

bool test_is_same(int x) {
    int p[2];
    p[0] = -x;
    p[1] = x;
    MPI_Allreduce(MPI_IN_PLACE, p, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return (p[0] == -p[1]);
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P, ID;
    MPI_Status status;
    MPI_Win win1;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    // int point_count = 10243;
    // int lbid, ubid;
    //
    // vector<double> particles (P, int(point_count / P));
    // vector<double> lb (P, 0);
    // vector<double> ub (P, 0);
    // int total = P * int(point_count / P);
    // int gap = point_count - total;
    // for (int i = 1; i < gap + 1; i++) {
    //     particles[i] += 1;
    // }
    // total = 0;
    // for (int i = 0; i < P; i++) {
    //     total += particles[i];
    // }
    // // assertm(total == point_count, "Particle count not correct");
    //
    // // for (int )
    // ub[0] = particles[0];
    // for (int i = 1; i < P; i++) {
    //     lb[i] = ub[i-1];
    //     ub[i] = lb[i] + particles[i];
    // }
    // lbid = lb[ID];
    // ubid = ub[ID];
    // cout << "Process: " << ID << " lb: " << lbid << " ub: " << ubid << " particles: " << ubid - lbid << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Processors: " << P << " rank: " << ID << endl;
    int length = 100000000;
    vector<double> testvec(length, 0);
    vector<double> testvec2(length, 1);
    chrono::steady_clock::time_point begin_time, end_time;
    if (ID == 0) begin_time = chrono::steady_clock::now();

    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Win_create(&testvec[0], length * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win1);
    MPI_Win_fence(0, win1);
    // if (ID != 0) {
    // for (int i = 0; i < P; i++) {
    //     MPI_Accumulate(&testvec2[0], length, MPI_DOUBLE, i, 0, length, MPI_DOUBLE, MPI_SUM, win1);
    // }
    MPI_Accumulate(&testvec2[0], length, MPI_DOUBLE, 0, 0, length, MPI_DOUBLE, MPI_SUM, win1);
    // first 0 is target rank, second 0 is displacement
    // }
    // testvec = testvec2;
    MPI_Win_fence(0, win1);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Bcast(&testvec[0], length, MPI_DOUBLE, 0, MPI_COMM_WORLD); //  1515801 microseconds for bcast
    // MPI_Barrier(MPI_COMM_WORLD);
    if (ID != 0) {
        MPI_Get(&testvec[0], length, MPI_DOUBLE, 0, 0, length, MPI_DOUBLE, win1); // 1139048 microseconds for gets
    }
    // if (ID == 0) {
    //     for (int i = 1; i < P; i++) {
    //         MPI_Put(&testvec[0], length, MPI_DOUBLE, i, 0, length, MPI_DOUBLE, win1); // 1102637.2 for puts
    //     }
    // }
    MPI_Win_fence(0, win1);
    if (ID == 0) {
        chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
        cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end_time - begin_time).count() << " microseconds" << endl;
    }
    double sum = 0;
    // sum = 0;
    for (int i = 0; i < length; i++) {
        sum += testvec[i];
    }
    // cout
    // bool issame;
    // issame = test_is_same(sum);
    // cout << issame << endl;
    // MPI_Barrier(MPI_COMM_WORLD);
    cout << "Processors: " << P << " rank: " << ID << " sum of vec: " << sum << endl;

    // free(&win1);
    MPI_Win_free(&win1);
    MPI_Finalize();
    return 0;
}
