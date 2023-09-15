//  ==
//  ||
//  ||   C L A S S:    m p i I n f o
//  ||
//  ==

class mpiInfo
{
public:
  int myPE;
  int numPE;
  int nRealx, nRealy;
  int nx, ny;
  int nPEx, nPEy;
  int iPE, jPE;
  int iMin, iMax, jMin, jMax; // The global i-j numbers on this processor
  int nei_n, nei_s, nei_e, nei_w;
  int nei_ne, nei_nw, nei_se, nei_sw;

  int countx, county;

  MPI_Status status;
  int err;
  int tag;
  MPI_Request request;

  //  -
  //  |
  //  |   GridDecomposition: Set up PE numbering system in figure below and
  //  |                      establish communication arrays.
  //  |
  //  |                      nPEx -- number of PEs in the x-direction
  //  |                      nPEy -- number of PEs in the y-direction
  //  |                      numPE = total number of PEs
  //  |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |       |       |         | numPE |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       | nPEx  | nPEx+1|         |       |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |   0   |   1   |         | nPEx-1|
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |
  //  |
  //  -

  void GridDecomposition(int _nPEx, int _nPEy, int nCellx, int nCelly)
  {

    nRealx = nCellx + 1;
    nRealy = nCelly + 1;

    // Store and check incoming processor counts

    nPEx = _nPEx;
    nPEy = _nPEy;

    if (nPEx * nPEy != numPE)
    {
      if (myPE == 0)
        cout << "Fatal Error:  Number of PEs in x-y directions do not add up to numPE" << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
    }

    // Get the i-j location of this processor, given its number.  See figure above:

    jPE = int(myPE / nPEx);
    iPE = myPE - jPE * nPEx;

    // Set neighbor values

    nei_n = nei_s = nei_e = nei_w = -1;

    if (iPE > 0)
    {
      nei_w = myPE - 1;
    }
    if (jPE > 0)
    {
      nei_s = myPE - nPEx;
    }
    if (iPE < nPEx - 1)
    {
      nei_e = myPE + 1;
    }
    if (jPE < nPEy - 1)
    {
      nei_n = myPE + nPEx;
    }

    nei_nw = nei_sw = nei_ne = nei_se = -1;

    if (iPE > 0 && jPE > 0)
      // nei_sw = /* TO-DO */;
      nei_sw = myPE - nPEx - 1;
    if (iPE < nPEx - 1 && jPE > 0)
      // nei_se = /* TO-DO */;
      nei_se = myPE - nPEx + 1;
    if (iPE > 0 && jPE < nPEy - 1)
      // nei_nw = /* TO-DO */;
      nei_nw = myPE + nPEx - 1;
    if (iPE < nPEx - 1 && jPE < nPEy - 1)
      // nei_ne = /* TO-DO */;
      nei_ne = myPE + nPEx + 1;

    // Acquire memory for the communication between adjacent processors:
    countx = nRealx + 2;
    county = nRealy + 2;

    tag = 0;
  }

  //  ==
  //  ||
  //  ||  ParticlesExchange
  //  ||
  //  ||  Exchange particles between processors
  //  ||
  //  ||
  //  ==

  void ParticleExchange(VI &ptcl_send_list, VI &ptcl_send_PE, particles &PTCL)
  {
    MPI_Status status;
    MPI_Request request;

    // (1) Get the max number particles to be sent by any particular processor, and make sure all processors  know that number.

    int numToSend = ptcl_send_list.size();
    int maxToSend;

    /* TO-DO Place MPI call here to complete (1) */
    // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
    // Combines values from all processes and distributes the result back to all processes.
    MPI_Allreduce(&numToSend, &maxToSend, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // (2) Allocate contributions to the upcoming Gather operation.  Here, "C" for "Contribution" to be Gathered

    int *Cptcl_PE;
    Cptcl_PE = new int[maxToSend /* TO-DO */]; // Particles' destination PEs
    double *Cptcl_x;                           // x - postion
    Cptcl_x = new double[maxToSend /* TO-DO */];
    double *Cptcl_y; // y - postion
    Cptcl_y = new double[maxToSend /* TO-DO */];
    double *Cptcl_vx; // x - velocity
    Cptcl_vx = new double[maxToSend /* TO-DO */];
    double *Cptcl_vy; // y - volocity
    Cptcl_vy = new double[maxToSend /* TO-DO */];

    // (3) Populate contributions on all processors for the upcoming Gather operation

    for (int i = 0; i < maxToSend; ++i)
    {
      Cptcl_PE[i] = -1;
      Cptcl_x[i] = 0.;
      Cptcl_y[i] = 0.;
      Cptcl_vx[i] = 0.;
      Cptcl_vy[i] = 0.;
    }

    // (4) Populate with all the particles on this PE.  Note that some/most processors will have left-over space in the C* arrays.

    for (int i = 0; i < ptcl_send_list.size(); ++i)
    {
      int id = ptcl_send_list[i /* TO-DO */];
      Cptcl_PE[i] = ptcl_send_PE[i /* TO-DO */];
      Cptcl_x[i] = PTCL.x[id /* TO-DO */];
      Cptcl_y[i] = PTCL.y[id /* TO-DO */];
      Cptcl_vx[i] = PTCL.vx[id /* TO-DO */];
      Cptcl_vy[i] = PTCL.vy[id /* TO-DO */];
    }

    // (5) Allocate and initialize the arrays for upcoming Gather operation to PE0.  The sizeOfGather takes
    //     into account the number of processors, like this figure:
    //
    //     |<----------------------------- sizeOfGather ------------------------------>|
    //     |                                                                           |
    //     |                                                                           |
    //     |<- maxToSend    ->|<- maxToSend    ->|<- maxToSend    ->|<- maxToSend    ->|
    //     +------------------+------------------+------------------+------------------+
    //             PE0               PE1                PE2               PE3

    int sizeOfGather = maxToSend * numPE /* TO-DO */;

    int *Gptcl_PE;
    Gptcl_PE = new int[sizeOfGather /* TO-DO */];
    double *Gptcl_x;
    Gptcl_x = new double[sizeOfGather /* TO-DO */];
    double *Gptcl_y;
    Gptcl_y = new double[sizeOfGather /* TO-DO */];
    double *Gptcl_vx;
    Gptcl_vx = new double[sizeOfGather /* TO-DO */];
    double *Gptcl_vy;
    Gptcl_vy = new double[sizeOfGather /* TO-DO */];

    for (int i = 0; i < sizeOfGather; ++i)
    {
      Gptcl_PE[i] = -1;
      Gptcl_x[i] = 0.;
      Gptcl_y[i] = 0.;
      Gptcl_vx[i] = 0.;
      Gptcl_vy[i] = 0.;
    }

    // (6)  Gather "Contributions" ("C" arrays) from all PEs onto all PEs into these bigger arrays so all PE will know what particles
    //      need to go where.

    // https://rookiehpc.org/mpi/docs/mpi_barrier/index.html
    MPI_Barrier(MPI_COMM_WORLD); // Synchronization between MPI processes

    MPI_Iallgather(/* TO-DO */); // MPI_Allgather, MPI_Iallgather - Gathers data from all processes and distributes it to all processes

    // https://rookiehpc.org/mpi/docs/mpi_iallgather/index.html
    // MPI_Allgather, MPI_Iallgather - Gathers data from all processes and distributes it to all processes
    // int MPI_Iallgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)

    // Cptcl_PE
    // Cptcl_x
    // Cptcl_y
    // Cptcl_vx
    // Cptcl_vy

    MPI_Iallgather(Cptcl_PE, maxToSend, MPI_INT, Gptcl_PE, maxToSend, MPI_INT, MPI_COMM_WORLD, &request /* TO-DO */);
    MPI_Wait(&request, &status);

    MPI_Iallgather(Cptcl_x, maxToSend, MPI_DOUBLE, Gptcl_x, maxToSend, MPI_DOUBLE, MPI_COMM_WORLD, &request /* TO-DO */);
    MPI_Wait(&request, &status);

    MPI_Iallgather(Cptcl_y, maxToSend, MPI_DOUBLE, Gptcl_y, maxToSend, MPI_DOUBLE, MPI_COMM_WORLD, &request /* TO-DO */);
    MPI_Wait(&request, &status);

    MPI_Iallgather(Cptcl_vx, maxToSend, MPI_DOUBLE, Gptcl_vx, maxToSend, MPI_DOUBLE, MPI_COMM_WORLD, &request /* TO-DO */);
    MPI_Wait(&request, &status);

    MPI_Iallgather(Cptcl_vy, maxToSend, MPI_DOUBLE, Gptcl_vy, maxToSend, MPI_DOUBLE, MPI_COMM_WORLD, &request /* TO-DO */);
    MPI_Wait(&request, &status);

    MPI_Barrier(MPI_COMM_WORLD);

    // (7) Put in vector form so they can be added to PTCL.  These arrays are 1-based.

    int Np = 0;
    for (int i = 0; i < sizeOfGather /* TO-DO */; ++i)
      if (Gptcl_PE[i] == myPE)
        ++Np;

    VD std_add_x;
    std_add_x.resize(Np + 1); //make size to length of Np plus one.
    VD std_add_y;
    std_add_y.resize(Np + 1);
    VD std_add_vx;
    std_add_vx.resize(Np + 1);
    VD std_add_vy;
    std_add_vy.resize(Np + 1);


    // Because of Np + 1. Count starts at 1. Why? We'll never know 
    // Dennis: :) 
    // Paul: ;)
    // Both :'(

    int count = 1;
    for (int i = 0; i < sizeOfGather; ++i)
      if (Gptcl_PE[i] == myPE/* TO-DO */)
      {
        std_add_x[count /* TO-DO */] = Gptcl_x[i];
        std_add_y[count /* TO-DO */] = Gptcl_y[i];
        std_add_vx[count /* TO-DO */] = Gptcl_vx[i];
        std_add_vy[count /* TO-DO */] = Gptcl_vy[i];
        ++count /* TO-DO */;
      }

    PTCL.add(std_add_x, std_add_y, std_add_vx, std_add_vy/* TO-DO */);

    // (8) Free up memory

    if (maxToSend > 0)
    {
      delete[] Cptcl_PE;
      delete[] Cptcl_x;
      delete[] Cptcl_y;
      delete[] Cptcl_vx;
      delete[] Cptcl_vy;
    }
    if (sizeOfGather > 0)
    {
      delete[] Gptcl_PE;
      delete[] Gptcl_x;
      delete[] Gptcl_y;
      delete[] Gptcl_vx;
      delete[] Gptcl_vy;
    }
  }

  int pid(int i, int j) { return (i + 1) + (j) * (nRealx + 2); }
};
