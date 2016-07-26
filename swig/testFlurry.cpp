#include <valgrind/callgrind.h>

#include "flurry_interface.hpp"
#include "tiogaInterface.h"

#include "mpi.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
CALLGRIND_STOP_INSTRUMENTATION;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int gridID = 0;
  int nGrids = 2;
  if (size == 8)
  {
    gridID = (rank >= 3);
  }
  else if (size == 3)
  {
    gridID = (rank > 0);
  }
  else if (size == 1)
  {
    gridID = 1;
    nGrids = 1;
  }
  else
  {
    gridID = rank%nGrids; //(rank > (size/2));
  }

  gridID = rank>0;
//  // 2-sphere test case
//  nGrids = 3;
//  if (rank == 0) gridID = 0;
//  if (rank == 1) gridID = 1;
//  if (rank > 1) gridID = 2;

  bool sphereTest = true;
  if (sphereTest)
  {
    if (nGrids == 2)
      gridID = (rank>0);
    else if (nGrids == 3)
      gridID = (rank>0);
  }

  MPI_Comm gridComm;

  bool oneGrid = false;
  if (oneGrid)
  {
    nGrids = 1; gridID = 0;
    gridComm = MPI_COMM_WORLD;
  }
  else
  {
    MPI_Comm_split(MPI_COMM_WORLD, gridID, rank, &gridComm);
  }
  cout << "Rank " << rank << ", GridID = " << gridID << ", nproc = " << size << endl;

  // Setup the Flurry solver object
  char inputFile[] = "input_sphere_overset";
  flurry::initialize(gridComm, inputFile, nGrids, gridID);

  Flurry *fr = flurry::get_flurry_object();
  fr->setup_solver();

  BasicGeo geo = flurry::get_basic_geo_data();
  ExtraGeo geoAB = flurry::get_extra_geo_data();
  CallbackFuncs cbs = flurry::get_callback_funcs();
  input &inp = fr->get_input();

  double *U_spts = flurry::get_q_spts();
  double *U_fpts = flurry::get_q_fpts();

  // Setup the TIOGA connectivity object
  tioga_init_(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
//  if (rank==0)
//  {
//    std::cout << "Before TIOGA setup - Press enter to continue... " << std::flush;
//    std::cin.clear();
//    std::cin.get();
//  }
  MPI_Barrier(MPI_COMM_WORLD);

  Timer tg_time;
  tg_time.startTimer();

  tioga_registergrid_data_(geo.btag, geo.nnodes, geo.xyz, geo.iblank,
      geo.nwall, geo.nover, geo.wallNodes, geo.overNodes,
      geo.nCellTypes, geo.nvert_cell, geo.nCells_type, geo.c2v);

  tioga_setcelliblank_(geoAB.iblank_cell);

  tioga_register_face_data_(geoAB.f2c,geoAB.c2f,geoAB.iblank_face,
      geoAB.nOverFaces,geoAB.nMpiFaces,geoAB.overFaces,geoAB.mpiFaces,
      geoAB.procR,geoAB.mpiFidR,geoAB.nFaceTypes,geoAB.nvert_face,
      geoAB.nFaces_type,geoAB.f2v);

  tioga_set_highorder_callback_(cbs.get_nodes_per_cell,
      cbs.get_receptor_nodes, cbs.donor_inclusion_test,
      cbs.donor_frac, cbs.convert_to_modal);

  tioga_set_ab_callback_(cbs.get_nodes_per_face, cbs.get_face_nodes,
      cbs.get_q_index_face, cbs.get_q_spt, cbs.get_q_fpt);

  if (nGrids > 1)
  {
    tioga_preprocess_grids_();
    tioga_performconnectivity_();
  }
  tg_time.stopTimer();

  // Output initial solution and grid
  fr->write_solution();

  MPI_Barrier(MPI_COMM_WORLD);
  CALLGRIND_START_INSTRUMENTATION;

  // Run the solver loop now
  Timer runTime("Compute Time: ");
  Timer tgTime("Interp Time: ");
  Timer mpiTime("MPI Wait Time: ");
  Timer totTime("Total run time: ");
//  inp.waitTimer = mpiTime;

  totTime.startTimer();

  while (inp.iter < inp.iterMax and inp.time < inp.maxTime)
  {
    tgTime.startTimer();
    if (nGrids > 1)
      tioga_dataupdate_ab(5,U_spts,U_fpts);
    tgTime.stopTimer();

    runTime.startTimer();
    fr->do_step();
    runTime.stopTimer();

    if (inp.iter == inp.initIter + 1) {fr->write_residual(); fr->write_error();}
    if (inp.iter%inp.monitorResFreq == 0 or inp.iter == inp.iterMax) fr->write_residual();
    if (inp.iter%inp.plotFreq == 0 or inp.iter == inp.iterMax) fr->write_solution();
    if (inp.monitorErrFreq > 0 and (inp.iter%inp.monitorErrFreq == 0 or inp.iter == inp.iterMax)) fr->write_error();
  }
  CALLGRIND_STOP_INSTRUMENTATION;

  totTime.stopTimer();

  fr->write_solution();

//  if (rank == 0)
//  {
    std::cout << "Preprocessing/Connectivity Time: ";
    tg_time.showTime(2);

    tgTime.showTime(2);
    runTime.showTime(2);
//  }

//    inp.waitTimer.showTime();

  MPI_Barrier(MPI_COMM_WORLD);
  totTime.showTime();

  flurry::finalize();
  tioga_delete_();

  MPI_Finalize();

  return 0;
}
