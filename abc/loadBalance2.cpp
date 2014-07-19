#include <mpi.h>
#include <vector>
#define DIETAG     10000000
#define MASTER_RANK 0

using namespace std;

void master();
void slave();

int main(int argc, char* argv[]) {
    int myrank;
    MPI_Init(&argc, &argv);   /* initialize MPI */
    MPI_Comm_rank(
      MPI_COMM_WORLD,   /* always use this */
      &myrank);      /* process rank, 0 thru N-1 */
    if (myrank == 0) {
        master();
    } else {
        slave();
    }
    MPI_Finalize();       /* cleanup MPI */
    return 0;
}

void master() {
    int ntasks = 10;
    int scheduledTasks=0;
    int mpi_size, rank, work;
    vector<double> result(ntasks);
    double buffer;
    MPI_Status status;
    MPI_Comm_size(
      MPI_COMM_WORLD,   /* always use this */
      &mpi_size);          /* #processes in application */
    /*
     * Seed the slaves.
     */

    int workINT = 0;
    for (rank = 1; rank < mpi_size; ++rank) {
        work = workINT++ ;      /* get_next_work_request */;
        MPI_Send(&work,         /* message buffer */
                1,              /* one data item */
                MPI_INT,        /* data item is an integer */
                rank,           /* destination process rank */
                scheduledTasks, /* user chosen message tag */
                MPI_COMM_WORLD);/* always use this */
        scheduledTasks++;
    }

    /*
     * Receive a result from any slave and dispatch a new work
     * request work requests have been exhausted.
     */
    //work = /* get_next_work_request */;
    work = workINT++;
    while ( workINT <= ntasks ) { 
        MPI_Recv(&buffer,       /* message buffer */
                1,              /* one data item */
                MPI_DOUBLE,     /* of type double real */
                MPI_ANY_SOURCE, /* receive from any sender */
                MPI_ANY_TAG,    /* any type of message */
                MPI_COMM_WORLD, /* always use this */
                &status);       /* received message info */

        result[status.MPI_TAG] = buffer;

        cerr << "work: " << work << endl;
        cerr << "task: " << status.MPI_TAG << " " << buffer << endl;
        MPI_Send(&work, 1, MPI_INT, status.MPI_SOURCE,
                scheduledTasks, MPI_COMM_WORLD);
        scheduledTasks++;
        work = workINT++; /* get_next_work_request */;
    }
    /*
     * Receive results for outstanding work requests.
     */
    for (rank = 1; rank < mpi_size; ++rank) {
        MPI_Recv(&buffer, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
                MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        result[status.MPI_TAG] = buffer;
        cerr << "task: " << status.MPI_TAG << " " << buffer << endl;
    }
    /*
     * Tell all the slaves to exit.
     */
    cerr << "almost done\n";
    for (rank = 1; rank < mpi_size; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }
//    for (auto x: result) {cout << x;}
//    cout << endl;
}

void slave() {
    double              result;
    int                 work;
    MPI_Status          status;
    while (1) {
        MPI_Recv(&work, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);
        /*
         * Check the tag of the received message.
         */
        if (status.MPI_TAG == DIETAG) {
            return;
        }
        result = work*work; /* do the work */;
//        cerr << result << endl;
        MPI_Send(&result, 1, MPI_DOUBLE, MASTER_RANK, status.MPI_TAG, MPI_COMM_WORLD);
    }
}
