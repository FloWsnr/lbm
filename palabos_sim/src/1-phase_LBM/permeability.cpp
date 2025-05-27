#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

typedef double T;
#define DESCRIPTOR plb::descriptors::D3Q19Descriptor

// This function object returns a zero velocity, and a pressure which decreases
// linearly in x-direction. It is used to initialize the particle populations.
class PressureGradient {
public:
  PressureGradient(T deltaP_, plb::plint nx_) : deltaP(deltaP_), nx(nx_) {}

  void operator()(plb::plint iX, plb::plint iY, plb::plint iZ, T &density,
                  plb::Array<T, 3> &velocity) const {
    velocity.resetToZero();
    density = (T)1 - deltaP * DESCRIPTOR<T>::invCs2 / (T)(nx - 1) * (T)iX;
  }

private:
  T deltaP;
  plb::plint nx;
};

// This function grabs the appropriate geometry for the single-phase simulation
void readGeometry(const std::string &fNameIn, const std::string &fNameOut,
                  plb::MultiScalarField3D<int> &geometry, plb::plint run,
                  plb::plint runnum, bool vtk_out,
                  const std::string &GeometryName) {
  const plb::plint nx = geometry.getNx();
  const plb::plint ny = geometry.getNy();
  const plb::plint nz = geometry.getNz();
  plb::plint run_diff = ((runnum - 1) / 2) + 1;
  std::string fNameIn_temp;

  plb::pcout << "\nRun: " << run << std::endl;

  plb::Box3D sliceBox(0, 0, 0, ny - 1, 0, nz - 1);

  // Selects between the original geometry or the fluid 1,2 final config from
  // multiphase sim
  if (run == 1) { // original geometry - absolute permeability
    fNameIn_temp = fNameIn + GeometryName + ".dat";
    plb::pcout << "Running absolute permeability" << std::endl;
  } else if (run > run_diff) { // Fluid 1 : Krnw
    const plb::plint runner = run - run_diff;
    fNameIn_temp = fNameIn + "f1_for_kr_" + std::to_string(runner) + ".dat";
    plb::pcout << "Running kr_f1" << std::endl;
  } else { // Fluid 2 : Krw
    fNameIn_temp = fNameIn + "f2_for_kr_" + std::to_string(run - 1) + ".dat";
    plb::pcout << "Running kr_f2" << std::endl;
  }

  plb::pcout << "The geometry name is " << fNameIn_temp << std::endl;

  std::unique_ptr<plb::MultiScalarField3D<int>> slice =
      plb::generateMultiScalarField<int>(geometry, sliceBox);
  plb::plb_ifstream geometryFile(fNameIn_temp.c_str());

  for (plb::plint iX = 0; iX < nx - 1; ++iX) {
    if (!geometryFile.is_open()) {
      plb::pcout << "Error: could not open the geometry file " << fNameIn_temp
                 << std::endl;
      std::exit(EXIT_FAILURE);
    }

    geometryFile >> *slice;
    plb::copy(*slice, slice->getBoundingBox(), geometry,
              plb::Box3D(iX, iX, 0, ny - 1, 0, nz - 1));
  }

  if (vtk_out) {
    plb::VtkImageOutput3D<T> vtkOut(plb::createFileName("PorousMedium", run, 6),
                                    1.0);
    vtkOut.writeData<float>(
        *plb::copyConvert<int, T>(geometry, geometry.getBoundingBox()), "tag",
        1.0);
  }
}

void porousMediaSetup(
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
    plb::OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition,
    plb::MultiScalarField3D<int> &geometry, T deltaP) {
  const plb::plint nx = lattice.getNx();
  const plb::plint ny = lattice.getNy();
  const plb::plint nz = lattice.getNz();

  plb::pcout << "Definition of inlet/outlet." << std::endl;

  plb::Box3D inlet(0, 0, 1, ny - 2, 1, nz - 2);
  boundaryCondition->addPressureBoundary0N(inlet, lattice);
  plb::setBoundaryDensity(lattice, inlet, (T)1.);

  plb::Box3D outlet(nx - 1, nx - 1, 1, ny - 2, 1, nz - 2);
  boundaryCondition->addPressureBoundary0P(outlet, lattice);
  plb::setBoundaryDensity(lattice, outlet,
                          (T)1. - deltaP * DESCRIPTOR<T>::invCs2);

  // Where "geometry" evaluates to 1, use bounce-back.
  plb::defineDynamics(lattice, geometry, new plb::BounceBack<T, DESCRIPTOR>(),
                      1);
  // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
  plb::defineDynamics(lattice, geometry, new plb::NoDynamics<T, DESCRIPTOR>(),
                      2);

  plb::initializeAtEquilibrium(lattice, lattice.getBoundingBox(),
                               PressureGradient(deltaP, nx));
  lattice.initialize();
  delete boundaryCondition;
}

void writeGifs(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
               plb::plint run) {
  const plb::plint nx = lattice.getNx();
  const plb::plint ny = lattice.getNy();
  const plb::plint nz = lattice.getNz();

  const plb::plint imSize = 600;
  plb::ImageWriter<T> imageWriter("leeloo");

  // Write velocity-norm at x=1.
  imageWriter.writeScaledGif(
      plb::createFileName("ux_inlet", run, 6),
      *plb::computeVelocityNorm(lattice,
                                plb::Box3D(2, 2, 0, ny - 1, 0, nz - 1)),
      imSize, imSize);

  // Write velocity-norm at x=nx/2.
  imageWriter.writeScaledGif(
      plb::createFileName("ux_half", run, 6),
      *plb::computeVelocityNorm(
          lattice, plb::Box3D(nx / 2, nx / 2, 0, ny - 1, 0, nz - 1)),
      imSize, imSize);
}

void writeVTK(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
              plb::plint run) {
  plb::VtkImageOutput3D<T> vtkOut(plb::createFileName("vtk_vel", run, 6), 1.);
  vtkOut.writeData<float>(*plb::computeVelocityNorm(lattice), "velocityNorm",
                          1.);
  vtkOut.writeData<3, float>(*plb::computeVelocity(lattice), "velocity", 1.);
}

void computePermeability(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice, T nu,
                         T deltaP, plb::Box3D domain, T &perm, T &meanU) {
  // Compute only the x-direction of the velocity (direction of the flow).
  plb::plint xComponent = 0;
  plb::plint nx = lattice.getNx();
  plb::plint ny = lattice.getNy();
  plb::plint nz = lattice.getNz();
  plb::Box3D domain1(0, nx - 1, 0, ny - 1, 0, nz - 1);

  meanU = plb::computeAverage(
      *plb::computeVelocityComponent(lattice, domain1, xComponent));
  perm = nu * meanU / (deltaP / (T)(nx - 1));

  plb::pcout << "Average velocity = " << meanU << std::endl;
  plb::pcout << "Permeability = " << perm << std::endl;
}

int main(int argc, char **argv) {
  plb::plbInit(&argc, &argv);

  std::string fNameIn, fNameOut;
  plb::plint nx, ny, nz;
  T deltaP;
  T Run;
  bool nx_p, ny_p, nz_p;
  bool vtk_out;
  std::string GeometryName;
  plb::plint maxT;
  T conv;

  std::string xmlFname;
  try {
    plb::global::argv(1).read(xmlFname);
  } catch (plb::PlbIOException &exception) {
    plb::pcout << "Wrong parameters; the syntax is: "
               << (std::string)plb::global::argv(0) << " input-file.xml"
               << std::endl;
    return -1;
  }

  // Read input parameters from the XML file.
  plb::pcout << "Reading inputs from xml file \n";
  try {
    plb::XMLreader document(xmlFname);
    document["geometry"]["file_geom"].read(GeometryName);
    document["geometry"]["size"]["x"].read(nx);
    document["geometry"]["size"]["y"].read(ny);
    document["geometry"]["size"]["z"].read(nz);
    document["geometry"]["per"]["x"].read(nx_p);
    document["geometry"]["per"]["y"].read(ny_p);
    document["geometry"]["per"]["z"].read(nz_p);
    document["folder"]["out_f"].read(fNameOut);
    document["folder"]["in_f"].read(fNameIn);
    document["simulations"]["press"].read(deltaP);
    document["simulations"]["num"].read(Run);
    document["simulations"]["iter"].read(maxT);
    document["simulations"]["conv"].read(conv);
    document["simulations"]["vtk_out"].read(vtk_out);
  } catch (plb::PlbIOException &exception) {
    plb::pcout << exception.what() << std::endl;
    return -1;
  }

  std::string inputF = fNameIn;
  plb::global::directories().setOutputDir(fNameOut + "/");
  plb::global::directories().setInputDir(inputF + "/");

  const T omega = 1.0;
  const T nu = ((T)1 / omega - (T)0.5) / DESCRIPTOR<T>::invCs2;
  const plb::plint runnum = Run;
  plb::plint run_diff = ((runnum - 1) / 2) + 1;

  std::vector<T> perm(runnum + 1);
  std::vector<T> meanU(runnum + 1);
  std::vector<T> rel_perm(runnum + 1);
  T Perm, Vel;

  plb::pcout << "Total simulations: " << runnum << std::endl;
  plb::pcout << "The convergence threshold is: " << conv << " %" << std::endl;

  for (plb::plint run = 1; run <= runnum; ++run) {
    plb::MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        nx, ny, nz, new plb::BGKdynamics<T, DESCRIPTOR>(omega));

    lattice.periodicity().toggle(0, nx_p);
    lattice.periodicity().toggle(1, ny_p);
    lattice.periodicity().toggle(2, nz_p);

    plb::MultiScalarField3D<int> geometry(nx, ny, nz);
    readGeometry(fNameIn, fNameOut, geometry, run, runnum, vtk_out,
                 GeometryName);

    porousMediaSetup(lattice,
                     plb::createLocalBoundaryCondition3D<T, DESCRIPTOR>(),
                     geometry, deltaP);

    plb::pcout << "Simulation begins" << std::endl;
    plb::plint iT = 0;
    T new_avg_f = 0, old_avg_f = 0, relE_f1;
    lattice.toggleInternalStatistics(false);

    for (; iT < maxT; ++iT) {
      if (iT % 250 == 0 && iT > 0) {
        lattice.toggleInternalStatistics(true);
        lattice.collideAndStream();
        new_avg_f = plb::getStoredAverageEnergy(lattice);
        lattice.toggleInternalStatistics(false);
        relE_f1 = std::fabs(old_avg_f - new_avg_f) * 100 / old_avg_f / 250;

        if (iT % 10000 == 0) {
          plb::pcout << "Iteration " << iT << std::endl;
          plb::pcout << "-----------------" << std::endl;
          plb::pcout << "Relative difference of Energy: "
                     << std::setprecision(3) << relE_f1 << " %" << std::endl;
          plb::pcout << "The preliminary permeability is: " << std::endl;
          computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(),
                              Perm, Vel);
          plb::pcout << "**********************************************"
                     << std::endl;
        }

        if (relE_f1 < conv) {
          break;
        }
        old_avg_f = new_avg_f;
      }
    }

    plb::pcout << "End of simulation at iteration " << iT << " for Run " << run
               << std::endl;
    computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(), Perm,
                        Vel);

    writeGifs(lattice, run);
    std::string outDir = fNameOut + "/";
    std::string vel_name = outDir + GeometryName + "_vel.dat";
    plb::plb_ofstream ofile3(vel_name.c_str());
    ofile3 << std::setprecision(1) << *plb::computeVelocity(lattice)
           << std::endl;

    perm[run] = Perm;
    meanU[run] = Vel;
    rel_perm[run] = perm[run] / perm[1];

    if (run == 1) {
      plb::pcout << "Absolute Permeability = " << perm[run] << std::endl;
    }
    plb::pcout << "Relative Permeability = " << rel_perm[run] << std::endl;

    if (vtk_out) {
      plb::pcout << "Writing VTK file ..." << std::endl;
      writeVTK(lattice, run);
    }
  }

  plb::pcout << "Printing outputs" << std::endl;
  std::string outDir = fNameOut + "/";
  std::string output = outDir + "relPerm&vels.txt";
  plb::plb_ofstream ofile(output.c_str());
  ofile << "Outputs\n\n";
  ofile << "Krw from run: \nKrnw from run: " << (run_diff + 1) << std::endl;

  for (plb::plint runs = 1; runs <= runnum; ++runs) {
    ofile << "Run = " << runs << std::endl;
    if (runs == 1) {
      ofile << "Absolute Permeability = " << perm[runs] << std::endl;
    }
    ofile << "Relative Permeability = " << rel_perm[runs] << std::endl;
    ofile << "Mean Velocity = " << meanU[runs] << std::endl;
  }

  return 0;
}