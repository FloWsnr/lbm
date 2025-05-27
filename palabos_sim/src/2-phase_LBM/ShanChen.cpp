#include "palabos3D.h"
#include "palabos3D.hh"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>

// Use double-precision arithmetics for all floating point calculations
using T = double;

// Use a grid which additionally to the f's stores two variables for the
// external force term.
#define DESCRIPTOR plb::descriptors::ForcedShanChenD3Q19Descriptor

// creates the gifs
void writeGif_f1(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1,
                 plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2,
                 std::string runs, plb::plint iT) {
  const plb::plint imSize = 600;
  const plb::plint nx = lattice_fluid2.getNx();
  const plb::plint ny = lattice_fluid2.getNy();
  const plb::plint nz = lattice_fluid2.getNz();
  plb::Box3D slice(0, nx, 0, ny, nz / 2, nz / 2);

  std::string im_name;

  im_name = "rho_f1_";
  im_name.append(runs);
  im_name.append("_");

  plb::ImageWriter<T> imageWriter("leeloo.map");
  imageWriter.writeScaledGif(plb::createFileName(im_name, iT, 8),
                             *plb::computeDensity(lattice_fluid1, slice),
                             imSize, imSize);

  return;
}

/**
 * Write velocity data to VTK file
 *
 * Parameters
 * ----------
 * lattice_fluid : plb::MultiBlockLattice3D<T, DESCRIPTOR>
 *     The lattice containing velocity data
 * im_name : std::string
 *     Base name for the output file
 * iter : plb::plint
 *     Iteration number for the output file
 * write_full_velocity : bool, optional
 *     If true, writes both velocity norm and full vector. If false, writes only
 * x-component append_suffix : bool, optional If true, appends "1_" to the
 * filename
 */
void writeVTK_vel(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid,
                  std::string im_name, plb::plint iter,
                  bool write_full_velocity = false, bool append_suffix = true) {
  const plb::plint nx = lattice_fluid.getNx();
  const plb::plint ny = lattice_fluid.getNy();
  const plb::plint nz = lattice_fluid.getNz();
  plb::Box3D domain(0, nx - 1, 0, ny - 1, 0, nz - 1);

  if (append_suffix) {
    im_name.append("1_");
  }

  plb::VtkImageOutput3D<T> vtkOut(plb::createFileName(im_name, iter, 8), 1.);

  if (write_full_velocity) {
    vtkOut.writeData<float>(*plb::computeVelocityNorm(lattice_fluid),
                            "velocityNorm", 1.);
    vtkOut.writeData<3, float>(*plb::computeVelocity(lattice_fluid), "velocity",
                               1.);
  } else {
    plb::plint xComponent = 0;
    vtkOut.writeData<T>(
        (*plb::computeVelocityComponent(lattice_fluid, domain, xComponent)),
        "Velocity", 1.);
  }
}

void writeVTK_rho(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid,
                  std::string im_name, plb::plint iter) {
  const plb::plint nx = lattice_fluid.getNx();
  const plb::plint ny = lattice_fluid.getNy();
  const plb::plint nz = lattice_fluid.getNz();

  im_name.append("1_");

  plb::VtkImageOutput3D<T> vtkOut(plb::createFileName(im_name, iter, 8), 1.);
  vtkOut.writeData<T>((*plb::computeDensity(lattice_fluid)), "Density", 1.);

  return;
}

T computeVelocity_f1(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1,
                     T nu_f1) {
  plb::plint xComponent = 0;
  const plb::plint nx = lattice_fluid1.getNx();
  const plb::plint ny = lattice_fluid1.getNy();
  const plb::plint nz = lattice_fluid1.getNz();

  plb::Box3D domain(0, nx - 1, 0, ny - 1, 0, nz - 1);
  T meanU1 = plb::computeAverage(
      *plb::computeVelocityComponent(lattice_fluid1, domain, xComponent));
  plb::pcout << "Average velocity for fluid1 in x direction    = " << meanU1
             << std::endl;

  return meanU1;
}

T computeVelocity_f2(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2,
                     T nu_f2) {
  plb::plint xComponent = 0;
  const plb::plint nx = lattice_fluid2.getNx();
  const plb::plint ny = lattice_fluid2.getNy();
  const plb::plint nz = lattice_fluid2.getNz();

  plb::Box3D domain(0, nx - 1, 0, ny - 1, 0, nz - 1);
  T meanU2 = plb::computeAverage(
      *plb::computeVelocityComponent(lattice_fluid2, domain, xComponent));
  plb::pcout << "Average velocity for fluid2 in x direction    = " << meanU2
             << std::endl;

  return meanU2;
}

T computeCapillaryNumber_f1(
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1, T nu_f1) {

  T meanU1 = computeVelocity_f1(lattice_fluid1, nu_f1);

  // Ca = viscosity * velocity / surface tension... surface tension = 0.15 in
  // the model (See Young-Laplace example)
  T Ca_fluid1 = nu_f1 * meanU1 / 0.15;
  // plb::pcout << "Ca fluid 1    = " << Ca_fluid1 << std::endl;
  //  plb::pcout << "viscosity fluid 1    = " << nu_f1 << std::endl;

  return Ca_fluid1;
}

T computeCapillaryNumber_f2(
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2, T nu_f2) {

  T meanU2 = computeVelocity_f2(lattice_fluid2, nu_f2);

  // Ca = viscosity * velocity / surface tension... surface tension = 0.15 in
  // the model (See Young-Laplace example)
  T Ca_fluid2 = nu_f2 * meanU2 / 0.15;
  // plb::pcout << "Ca fluid 2    = " << Ca_fluid2 << std::endl;
  //  plb::pcout << "viscosity fluid 2    = " << nu_f2 << std::endl;

  return Ca_fluid2;
}

void readGeometry(std::string fNameIn, std::string fNameOut,
                  plb::MultiScalarField3D<int> &geometry) {
  const plb::plint nx = geometry.getNx();
  const plb::plint ny = geometry.getNy();
  const plb::plint nz = geometry.getNz();

  plb::Box3D sliceBox(0, 0, 0, ny - 1, 0, nz - 1);
  std::unique_ptr<plb::MultiScalarField3D<int>> slice =
      plb::generateMultiScalarField<int>(geometry, sliceBox);
  plb::plb_ifstream geometryFile(fNameIn.c_str());
  for (plb::plint iX = 0; iX < nx - 1; ++iX) {
    if (!geometryFile.is_open()) {
      plb::pcout << "Error: could not open geometry file " << fNameIn
                 << std::endl;
      std::exit(EXIT_FAILURE);
    }
    geometryFile >> *slice;
    plb::copy(*slice, slice->getBoundingBox(), geometry,
              plb::Box3D(iX, iX, 0, ny - 1, 0, nz - 1));
  }
  geometryFile.close();

  plb::VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
  vtkOut.writeData<T>(
      *plb::copyConvert<int, T>(geometry, geometry.getBoundingBox()), "tag",
      1.0);

  std::unique_ptr<plb::MultiScalarField3D<T>> floatTags =
      plb::copyConvert<int, T>(geometry, geometry.getBoundingBox());
  std::vector<T> isoLevels;
  isoLevels.push_back(0.5);
  typedef plb::TriangleSet<T>::Triangle Triangle;
  std::vector<Triangle> triangles;
  plb::Box3D domain = floatTags->getBoundingBox().enlarge(-1);
  domain.x0++;
  domain.x1--;
  plb::isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
  plb::TriangleSet<T> set(triangles);
  std::string outDir = fNameOut + "/";
  set.writeBinarySTL(outDir + "porousMedium.stl");

  return;
}

void setboundaryvalue(plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1,
                      plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2,
                      plb::Box3D inlet, plb::Box3D outlet, T rho_f1_inlet,
                      T rho_f2_outlet, T rhoNoFluid) {

  plb::setBoundaryDensity(lattice_fluid1, inlet, rho_f1_inlet);
  plb::setBoundaryDensity(lattice_fluid2, inlet, rhoNoFluid);
  plb::setBoundaryDensity(lattice_fluid1, outlet, rhoNoFluid);
  plb::setBoundaryDensity(lattice_fluid2, outlet, rho_f2_outlet);

  return;
}

void InitializeFluidsFromImage(
    plb::MultiScalarField3D<int> &geometry,
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1,
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2, T &rho_f1,
    T &rho_f2, T &rhoNoFluid, plb::Array<T, 3> &zeroVelocity) {

  const plb::plint nx = geometry.getNx();
  const plb::plint ny = geometry.getNy();
  const plb::plint nz = geometry.getNz();

  for (plb::plint iX = 0; iX < nx; iX++) {
    for (plb::plint iY = 0; iY < ny; iY++) {
      for (plb::plint iZ = 0; iZ < nz; iZ++) {

        plb::plint geom_value = geometry.get(iX, iY, iZ);

        if (geom_value == 0) {
          // If 0, set as fluid 2
          plb::initializeAtEquilibrium(lattice_fluid2,
                                       plb::Box3D(iX, iX, iY, iY, iZ, iZ),
                                       rho_f2, zeroVelocity);

          plb::initializeAtEquilibrium(lattice_fluid1,
                                       plb::Box3D(iX, iX, iY, iY, iZ, iZ),
                                       rhoNoFluid, zeroVelocity);

        } else if (geom_value == 3) {
          // If 3, set as fluid 1
          plb::initializeAtEquilibrium(lattice_fluid1,
                                       plb::Box3D(iX, iX, iY, iY, iZ, iZ),
                                       rho_f1, zeroVelocity);

          plb::initializeAtEquilibrium(lattice_fluid2,
                                       plb::Box3D(iX, iX, iY, iY, iZ, iZ),
                                       rhoNoFluid, zeroVelocity);
        }
      }
    }
  }
  return;
}

void PorousMediaSetup(
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid1,
    plb::MultiBlockLattice3D<T, DESCRIPTOR> &lattice_fluid2,
    plb::MultiScalarField3D<int> &geometry,
    plb::OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition,
    plb::Box3D inlet, plb::Box3D outlet, T rhoNoFluid, T rho_f1, T rho_f2,
    T rho_f1_inlet, T rho_f2_outlet, T Gads_f1_s1, T Gads_f1_s2, T Gads_f1_s3,
    T Gads_f1_s4, T force_f1, T force_f2, T nx1_f1, T nx2_f1, T ny1_f1,
    T ny2_f1, T nz1_f1, T nz2_f1, T nx1_f2, T nx2_f2, T ny1_f2, T ny2_f2,
    T nz1_f2, T nz2_f2, T runs, bool load_state, bool print_geom,
    bool pressure_bc, bool load_fluids_from_geom) {

  plb::pcout << "Definition of the geometry." << std::endl;

  plb::Array<T, 3> zeroVelocity(0., 0., 0.);

  if (pressure_bc == true) {

    // Inlet BC
    boundaryCondition->addPressureBoundary0N(inlet, lattice_fluid1);
    boundaryCondition->addPressureBoundary0N(inlet, lattice_fluid2);
    plb::setBoundaryDensity(lattice_fluid1, inlet,
                            rho_f1_inlet);                      // rho_f1_inlet
    plb::setBoundaryDensity(lattice_fluid2, inlet, rhoNoFluid); // rhoNoFluid

    // Outlet BC
    boundaryCondition->addPressureBoundary0P(outlet, lattice_fluid1);
    boundaryCondition->addPressureBoundary0P(outlet, lattice_fluid2);
    plb::setBoundaryDensity(lattice_fluid1, outlet, rhoNoFluid); // rhoNoFluid
    plb::setBoundaryDensity(lattice_fluid2, outlet,
                            rho_f2_outlet); // rho_f2_outlet
    // delete boundaryCondition;
  }

  // Assign masks to label geometry

  // NoDynamics (computational efficiency, label grains with 2)
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::NoDynamics<T, DESCRIPTOR>(), 2);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::NoDynamics<T, DESCRIPTOR>(), 2);

  // First contact angle (labeled with 1)
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(Gads_f1_s1), 1);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(-Gads_f1_s1), 1);

  // Second contact angle (labeled with 4)
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(Gads_f1_s2), 4);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(-Gads_f1_s2), 4);

  // Mesh contact angle (labeled with 5). Neutral wet
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(0), 5);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(0), 5);

  // Third contact angle (labeled with 6)
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(Gads_f1_s3), 6);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(-Gads_f1_s3), 6);

  // Fourth contact angle (labeled with 7)
  plb::defineDynamics(lattice_fluid1, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(Gads_f1_s4), 7);
  plb::defineDynamics(lattice_fluid2, geometry,
                      new plb::BounceBack<T, DESCRIPTOR>(-Gads_f1_s4), 7);

  // Array<T, 3> zeroVelocity(0., 0., 0.);

  // bool load_fluids_from_image = true; // For testing
  if (load_state == false) {
    plb::pcout << "Initializing Fluids" << std::endl;

    if (load_fluids_from_geom == false) {

      plb::initializeAtEquilibrium(lattice_fluid2,
                                   plb::Box3D(nx1_f2 - 1, nx2_f2 - 1,
                                              ny1_f2 - 1, ny2_f2 - 1,
                                              nz1_f2 - 1, nz2_f2 - 1),
                                   rho_f2, zeroVelocity);

      plb::initializeAtEquilibrium(lattice_fluid1,
                                   plb::Box3D(nx1_f2 - 1, nx2_f2 - 1,
                                              ny1_f2 - 1, ny2_f2 - 1,
                                              nz1_f2 - 1, nz2_f2 - 1),
                                   rhoNoFluid, zeroVelocity);

      plb::initializeAtEquilibrium(
          lattice_fluid1,
          plb::Box3D(nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1), rho_f1,
          zeroVelocity);

      plb::initializeAtEquilibrium(
          lattice_fluid2,
          plb::Box3D(nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1),
          rhoNoFluid, zeroVelocity);

    } else {
      plb::pcout << "Initialize fluid nodes from geom...";
      InitializeFluidsFromImage(geometry, lattice_fluid1, lattice_fluid2,
                                rho_f1, rho_f2, rhoNoFluid, zeroVelocity);
      plb::pcout << "Done!" << std::endl;
    }

    plb::setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(),
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                           plb::Array<T, 3>(force_f1, 0., 0.));
    plb::setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(),
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                           plb::Array<T, 3>(force_f2, 0., 0.));

    lattice_fluid1.initialize();
    lattice_fluid2.initialize();
  }

  // Output geometry dynamics
  if (print_geom == true) {
    plb::VtkImageOutput3D<int> vtkOut(plb::createFileName("vtkgeometry", 1, 1),
                                      1.);
    vtkOut.writeData<int>(geometry, "Dynamics", 1.);
    plb::pcout << "Creating geometry vtk file" << std::endl;
  }

  return;
}

int main(int argc, char *argv[]) {
  // 1. Declaring the variables
  clock_t t;
  t = clock();
  plb::plbInit(&argc, &argv);

  bool load_state;
  std::string fNameOut;
  std::string fNameIn;
  plb::plint nx, ny, nz;

  bool use_plb_bc;                               //
  bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; // periodicity
  bool pressure_bc;
  bool load_fluids_from_geom;
  plb::plint nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1,
      nz2_f1; // fluid1 configuration
  plb::plint nx1_f2, nx2_f2, ny1_f2, ny2_f2, nz1_f2,
      nz2_f2; // fluid2 configuration

  T G;
  T omega_f1;
  T omega_f2;
  T force_f1;
  T force_f2;
  T Gads_f1_s1;
  T Gads_f1_s2;
  T Gads_f1_s3;
  T Gads_f1_s4;

  T rho_f1;
  T rho_f2;

  T rho_f1_inlet;
  T rho_f2_outlet_initial;
  // T rho_f2_outlet_final;
  T rhoNoFluid;
  // T rho_f2_step;
  // T drho_f2;
  T num_pc_steps;
  T min_radius;

  plb::plint max_iter; // max total iterations
  plb::plint it_max_per_p;
  plb::plint it_conv;
  // plint it_info ;
  plb::plint it_vtk;
  plb::plint it_gif;
  plb::plint save_it;

  bool save_sim, rho_vtk, print_geom, print_stl;

  T convergence;

  std::string xmlFname;
  try {
    plb::global::argv(1).read(xmlFname);
  } catch (plb::PlbIOException &exception) {
    plb::pcout << "Wrong parameters; the syntax is: "
               << (std::string)plb::global::argv(0) << " input-file.xml"
               << std::endl;
    return -1;
  }

  // 2. Read input parameters from the XML file.
  plb::pcout << "Reading inputs from xml file \n";
  try {
    plb::XMLreader document(xmlFname);

    document["load_savedstated"].read(load_state);

    document["geometry"]["file_geom"].read(fNameIn);
    document["geometry"]["size"]["x"].read(nx);
    document["geometry"]["size"]["y"].read(ny);
    document["geometry"]["size"]["z"].read(nz);
    document["geometry"]["per"]["fluid1"]["x"].read(px_f1);
    document["geometry"]["per"]["fluid1"]["y"].read(py_f1);
    document["geometry"]["per"]["fluid1"]["z"].read(pz_f1);
    document["geometry"]["per"]["fluid2"]["x"].read(px_f2);
    document["geometry"]["per"]["fluid2"]["y"].read(py_f2);
    document["geometry"]["per"]["fluid2"]["z"].read(pz_f2);

    document["init"]["fluid_from_geom"].read(load_fluids_from_geom);
    document["init"]["fluid1"]["x1"].read(nx1_f1);
    document["init"]["fluid1"]["x2"].read(nx2_f1);
    document["init"]["fluid1"]["y1"].read(ny1_f1);
    document["init"]["fluid1"]["y2"].read(ny2_f1);
    document["init"]["fluid1"]["z1"].read(nz1_f1);
    document["init"]["fluid1"]["z2"].read(nz2_f1);

    document["init"]["fluid2"]["x1"].read(nx1_f2);
    document["init"]["fluid2"]["x2"].read(nx2_f2);
    document["init"]["fluid2"]["y1"].read(ny1_f2);
    document["init"]["fluid2"]["y2"].read(ny2_f2);
    document["init"]["fluid2"]["z1"].read(nz1_f2);
    document["init"]["fluid2"]["z2"].read(nz2_f2);

    document["fluids"]["Gc"].read(G);
    document["fluids"]["omega_f1"].read(omega_f1);
    document["fluids"]["omega_f2"].read(omega_f2);

    document["fluids"]["force_f1"].read(force_f1);
    document["fluids"]["force_f2"].read(force_f2);

    document["fluids"]["G_ads_f1_s1"].read(Gads_f1_s1);
    document["fluids"]["G_ads_f1_s2"].read(Gads_f1_s2);
    document["fluids"]["G_ads_f1_s3"].read(Gads_f1_s3);
    document["fluids"]["G_ads_f1_s4"].read(Gads_f1_s4);

    document["fluids"]["rho_f1"].read(rho_f1);
    document["fluids"]["rho_f2"].read(rho_f2);

    document["fluids"]["pressure_bc"].read(pressure_bc);

    document["fluids"]["rho_f1_i"].read(rho_f1_inlet);
    document["fluids"]["rho_f2_i"].read(rho_f2_outlet_initial);
    // document["fluids"]["rho_f2_f"].read(rho_f2_outlet_final);
    document["fluids"]["rho_d"].read(rhoNoFluid);
    // document["fluids"]["drho_f2"].read(drho_f2);
    document["fluids"]["num_pc_steps"].read(num_pc_steps);
    document["fluids"]["min_radius"].read(min_radius);

    document["output"]["out_folder"].read(fNameOut);
    document["output"]["save_sim"].read(save_sim);
    document["output"]["save_it"].read(save_it);
    document["output"]["convergence"].read(convergence);

    document["output"]["max_iter"].read(max_iter);
    document["output"]["it_max_per_p"].read(it_max_per_p);
    document["output"]["it_conv"].read(it_conv);

    document["output"]["it_gif"].read(it_gif);
    document["output"]["it_vtk"].read(it_vtk);
    document["output"]["rho_vtk"].read(rho_vtk);

    document["output"]["print_geom"].read(print_geom);
    document["output"]["print_stl"].read(print_stl);

  } catch (plb::PlbIOException &exception) {
    plb::pcout << exception.what() << std::endl;
    plb::pcout << exception.what() << std::endl;
    return -1;
  }

  plb::plint runnum = 0;
  if (pressure_bc == true) {
    // If running using pressure BCs, use the number of pressure steps specified
    runnum = num_pc_steps + 1;

    // Old method
    // runnum = ((rho_f2_outlet_initial - rho_f2_outlet_final) / drho_f2) + 1;
  } else {

    runnum = 1; // If not using pressure bc, set to 1 to save a one set of fluid
                // densities and a pressure value
  }

  plb::global::directories().setOutputDir(fNameOut);

  T rho_fluid1[runnum];
  T rho_fluid2[runnum];
  T deltaP[runnum];
  T new_avg_f1;
  T new_avg_f2;
  T old_avg_f1 = 1.0;
  T old_avg_f2 = 1.0;
  T relE_f1;
  T relE_f2;
  // T k1_high;
  // T k2_high;
  // T meanRho1;
  // T meanRho2;
  // T mu1;
  // T mu2;
  // T rho_F1;
  // T rho_F2;
  // T mean_U1[runnum];
  // T mean_U2[runnum];
  // T mean_rho1[runnum];
  // T mean_rho2[runnum];

  std::string outDir = fNameOut;
  std::string Lattice1 = fNameOut + "lattice1.dat";
  std::string Lattice2 = fNameOut + "lattice2.dat";

  if (pressure_bc == true) {

    // Calculating capillary pressure steps
    T cos_theta =
        abs(4 * Gads_f1_s1 /
            (G * (rho_f1_inlet -
                  rhoNoFluid))); // Taking absolute value so that the difference
                                 // in density is always positive
    T sigma = 0.15;              // tuning parameter from docs
    T delta_rho = 6 * sigma * cos_theta / min_radius;
    T step_size =
        (rho_f2_outlet_initial - (rho_f2_outlet_initial - delta_rho)) /
        num_pc_steps; // To calculate densities in the for loop

    for (plb::plint readnum = 0; readnum <= runnum; ++readnum) {
      rho_fluid2[readnum] = rho_f2_outlet_initial - readnum * step_size;
      // rho_fluid2[readnum] = rho_f2_outlet_initial - (readnum - 1) * drho_f2;
      //  plb::pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
      rho_fluid1[readnum] = rho_f1_inlet;
    }

  } else {
    plb::plint index = 0; // If not using pressure bc, set to 0 to fill in
                          // correct density values
    rho_fluid1[index] = rho_f1_inlet;
    rho_fluid2[index] = rho_f2_outlet_initial;
  }

  const T nu_f1 = ((T)1 / omega_f1 - 0.5) / DESCRIPTOR<T>::invCs2;
  const T nu_f2 = ((T)1 / omega_f2 - 0.5) / DESCRIPTOR<T>::invCs2;

  // Use regularized BGK dynamics to improve numerical stability
  // (but note that BGK dynamics works well too).
  plb::MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2(
      nx, ny, nz,
      new plb::ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega_f2));
  plb::MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1(
      nx, ny, nz,
      new plb::ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega_f1));

  lattice_fluid2.periodicity().toggle(0, px_f2);
  lattice_fluid1.periodicity().toggle(0, px_f1);
  lattice_fluid2.periodicity().toggle(1, py_f2);
  lattice_fluid1.periodicity().toggle(1, py_f1);
  lattice_fluid2.periodicity().toggle(2, pz_f2);
  lattice_fluid1.periodicity().toggle(2, pz_f1);

  std::vector<plb::MultiBlockLattice3D<T, DESCRIPTOR> *> blockLattices;
  blockLattices.push_back(&lattice_fluid2);
  blockLattices.push_back(&lattice_fluid1);

  std::vector<T> constOmegaValues;
  constOmegaValues.push_back(omega_f2);
  constOmegaValues.push_back(omega_f1);
  plb::plint processorLevel = 1;

  plb::integrateProcessingFunctional(
      new plb::ShanChenMultiComponentProcessor3D<T, DESCRIPTOR>(
          G, constOmegaValues),
      plb::Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), blockLattices,
      processorLevel);

  plb::pcout << "The convergence set by the user is = " << convergence
             << std::endl;

  if (pressure_bc == true) {
    plb::pcout << "The boundary conditions per run are:" << std::endl;
    for (plb::plint readnum = 0; readnum < runnum; ++readnum) {
      deltaP[readnum] = (rho_fluid1[readnum] - rho_fluid2[readnum]) / 3;
      plb::pcout << "Run number = " << readnum << std::endl;
      plb::pcout << "Rho_no_1 = " << rho_fluid1[readnum] << std::endl;
      plb::pcout << "Rho_no_2 = " << rho_fluid2[readnum] << std::endl;
    }
  }

  plb::pcout << "Reading the geometry file." << std::endl;
  plb::MultiScalarField3D<int> geometry(nx, ny, nz);
  readGeometry(fNameIn, fNameOut, geometry);

  plb::Box3D inlet(1, 2, 1, ny - 2, 1, nz - 2);
  plb::Box3D outlet(nx - 2, nx - 1, 1, ny - 2, 1, nz - 2);

  // Setup or load fluid lattices
  plb::plint current_run_num = 0;
  if (load_state == true) {
    // First check if run_num.dat is there and if so load it
    plb::pcout << "Check run_num.dat for restart info" << std::endl;
    // Check run_num.dat
    std::string runnum_file = outDir + "/run_num.dat";
    plb::plb_ifstream ifile(runnum_file.c_str());

    if (ifile.is_open()) {
      ifile >> current_run_num;
      plb::global::mpi().bCast(
          &current_run_num,
          1); // Broadcast so all the processors don't get confused!
      plb::pcout << "Current Run Number: " << current_run_num << std::endl;

    } else {
      plb::pcout
          << "No run_num.dat file found. Starting simulation from beginning."
          << std::endl;
      load_state = false;
    }
    ifile.close();
  }

  // If load_state still true, check restart files and load if they're there
  if (load_state == true) {

    plb::pcout << "Loading restart files..." << std::endl;
    try {
      plb::loadBinaryBlock(lattice_fluid1, Lattice1); // "lattice_fluid1.dat"
      plb::loadBinaryBlock(lattice_fluid2, Lattice2); // "lattice_fluid2.dat"
    } catch (plb::PlbIOException &exception) {
      throw std::runtime_error("Restart files not found.");
    }

    T rho_f1_inlet_new = rho_fluid1[current_run_num];
    T rho_f2_outlet_new = rho_fluid2[current_run_num];

    PorousMediaSetup(
        lattice_fluid1, lattice_fluid2, geometry,
        plb::createLocalBoundaryCondition3D<T, DESCRIPTOR>(), inlet, outlet,
        rhoNoFluid, rho_f1, rho_f2, rho_f1_inlet_new, rho_f2_outlet_new,
        Gads_f1_s1, Gads_f1_s2, Gads_f1_s3, Gads_f1_s4, force_f1, force_f2,
        nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1, nx1_f2, nx2_f2, ny1_f2,
        ny2_f2, nz1_f2, nz2_f2, current_run_num, load_state, print_geom,
        pressure_bc, load_fluids_from_geom);

    plb::pcout << "Starting the sim!" << std::endl;
  }

  // Otherwise, set-up a new simulation domain
  if (load_state == false) {

    plb::pcout << "Setting up a new simulation domain..." << std::endl;
    T rho_f1_inlet_new = rho_fluid1[current_run_num];
    T rho_f2_outlet_new = rho_fluid2[current_run_num];

    PorousMediaSetup(
        lattice_fluid1, lattice_fluid2, geometry,
        plb::createLocalBoundaryCondition3D<T, DESCRIPTOR>(), inlet, outlet,
        rhoNoFluid, rho_f1, rho_f2, rho_f1_inlet_new, rho_f2_outlet_new,
        Gads_f1_s1, Gads_f1_s2, Gads_f1_s3, Gads_f1_s4, force_f1, force_f2,
        nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1, nx1_f2, nx2_f2, ny1_f2,
        ny2_f2, nz1_f2, nz2_f2, current_run_num, load_state, print_geom,
        pressure_bc, load_fluids_from_geom);

    plb::pcout << "Starting the sim!" << std::endl;
  }

  use_plb_bc = true; // Use Palabos built-in BC
  // Initialize a global iteration counter
  plb::plint total_iter = 0;

  // Loop simulations with varying saturation
  for (plb::plint runs = current_run_num; runs < runnum; ++runs) {

    // turn off stats for efficiency
    lattice_fluid1.toggleInternalStatistics(false);
    lattice_fluid2.toggleInternalStatistics(false);

    // save a str for figure naming
    std::stringstream save_str;
    save_str << std::setw(3) << std::setfill('0') << runs;
    std::string runs_str;
    save_str >> runs_str;

    plb::pcout << "Run number = " << runs << std::endl;

    // re-use the final state of the previous run
    if (runs > current_run_num) {
      plb::pcout << "Using previous simulation state  " << std::endl;

      if (use_plb_bc == true && pressure_bc == true) {
        plb::pcout << "Updating constant bc pressure" << std::endl;
        setboundaryvalue(lattice_fluid1, lattice_fluid2, inlet, outlet,
                         rho_fluid1[runs], rho_fluid2[runs], rhoNoFluid);
      }
    }

    plb::pcout << std::endl
               << "Starting simulation with rho 1:  " << rho_fluid1[runs]
               << std::endl;
    plb::pcout << std::endl
               << "Starting simulation with rho 2:  " << rho_fluid2[runs]
               << std::endl;

    plb::plint checkconv = 0;

    plb::plint iT;
    if (runs == current_run_num && load_state == true) {
      std::string iter_file = outDir + "/iter_num.dat";
      plb::plb_ifstream ifile(iter_file.c_str());

      if (ifile.is_open()) {
        ifile >> iT;
        plb::global::mpi().bCast(
            &iT, 1); // Broadcast so all the processors don't get confused!
        plb::pcout << "Current iteration number: " << iT << std::endl;
        ifile.close();
      }
    } else {
      iT = 0;
    }

    while (checkconv == 0) { // Main loop over time iterations.
      iT = iT + 1;
      total_iter = total_iter + 1;
      // turn on stats to check convergence
      if (iT % it_conv == 0) {
        lattice_fluid1.toggleInternalStatistics(true);
        lattice_fluid2.toggleInternalStatistics(true);
      }

      lattice_fluid1.collideAndStream();
      lattice_fluid2.collideAndStream();

      // save gifs
      if (iT % it_gif == 0) {
        writeGif_f1(lattice_fluid1, lattice_fluid2, runs_str, iT);
      }

      // save vtks
      if (iT % it_vtk == 0) {
        writeVTK_rho(lattice_fluid1, "rho_f1_", iT);
        if (rho_vtk == true) {
          writeVTK_rho(lattice_fluid2, "rho_f2_", iT);
        }
        // Write full velocity data
        writeVTK_vel(lattice_fluid1, "vtk_vel_rho1_", iT, true, false);
      }

      // check total iters
      if (total_iter >= max_iter) {
        plb::pcout << "Reached maximum simulation iterations of " << max_iter
                   << ". Aborting simulation." << std::endl;

        // Save restart files before exiting
        if (save_sim == true) {
          plb::pcout << "Saving final restart files before aborting."
                     << std::endl;
          plb::saveBinaryBlock(lattice_fluid1, Lattice1);
          plb::saveBinaryBlock(lattice_fluid2, Lattice2);

          std::string run_name = outDir + "/run_num.dat";
          plb::plb_ofstream ofile_final_run(run_name.c_str());
          ofile_final_run << runs << std::endl;
          ofile_final_run.close();

          std::string iter_name = outDir + "/iter_num.dat";
          plb::plb_ofstream ofile_final_iter(iter_name.c_str());
          ofile_final_iter << iT << std::endl;
          ofile_final_iter.close();
        }
        exit(0);
      }

      // Check for maximum iterations per pressure step
      if (iT >= it_max_per_p) {
        plb::pcout << "Reached maximum iterations for this pressure step: "
                   << it_max_per_p << ". Moving to the next pressure step."
                   << std::endl;
        // Break the while loop to proceed to the next run
        checkconv = 1;
        // Saving the current state is done below
      }
      // ** Abort Conditions Section End **

      if (iT % it_conv == 0) {
        // calculate average change in mass if bcs == pressure
        new_avg_f1 =
            plb::getStoredAverageDensity(lattice_fluid1) * (nx * ny * nz);
        new_avg_f2 =
            plb::getStoredAverageDensity(lattice_fluid2) * (nx * ny * nz);

        if (pressure_bc == false) {
          // calculate average change in momentum if bcs == force
          new_avg_f1 = plb::getStoredAverageEnergy(lattice_fluid1);
          new_avg_f2 = plb::getStoredAverageEnergy(lattice_fluid2);
        }

        // mean_rho1[runs] = getStoredAverageDensity<T>(lattice_fluid1);
        // mean_rho2[runs] = getStoredAverageDensity<T>(lattice_fluid2);

        lattice_fluid1.toggleInternalStatistics(false);
        lattice_fluid2.toggleInternalStatistics(false);

        // calculate relative difference
        relE_f1 =
            std::fabs(old_avg_f1 - new_avg_f1) * 100 / old_avg_f1 / it_conv;
        relE_f2 =
            std::fabs(old_avg_f2 - new_avg_f2) * 100 / old_avg_f2 / it_conv;

        plb::pcout << "Run num " << runs;
        plb::pcout << ", Iteration " << iT << std::endl;
        plb::pcout << "-----------------" << std::endl;
        plb::pcout << "Relative difference average per iter fluid1: "
                   << std::setprecision(3) << relE_f1 << " %" << std::endl;
        plb::pcout << "Relative difference average per iter fluid2: "
                   << std::setprecision(3) << relE_f2 << " %" << std::endl;
        plb::pcout << "Has fluid 1 converged?: "
                   << ((relE_f1 < convergence) ? "TRUE" : "FALSE") << std::endl;
        plb::pcout << "Has fluid 2 converged?: "
                   << ((relE_f2 < convergence) ? "TRUE" : "FALSE") << std::endl;
        //        plb::pcout << "-----------------" << std::endl;

        // calculate capillary number
        T Ca_1, Ca_2;
        Ca_1 = computeCapillaryNumber_f1(lattice_fluid1, nu_f1);
        Ca_2 = computeCapillaryNumber_f2(lattice_fluid2, nu_f2);
        plb::pcout << "Ca fluid 1 = " << Ca_1 << std::endl;
        plb::pcout << "Ca fluid 2 = " << Ca_2 << std::endl;
        plb::pcout << "-----------------" << std::endl;

        // store new properties
        old_avg_f1 = new_avg_f1;
        old_avg_f2 = new_avg_f2;

        if (relE_f1 < convergence && relE_f2 < convergence) {
          checkconv = 1;
          plb::pcout << "Pressure increment has converged" << std::endl;
        }
      }

      // saves a binary file (heavy) with the sim state
      if (save_sim == true && iT > 0 && iT % save_it == 0) {
        // save only if checkconv == 0
        // this avoids double saving at the end of the sim
        if (checkconv == 0) {
          plb::pcout << "Saving restart files" << std::endl;
          plb::saveBinaryBlock(lattice_fluid1, Lattice1);
          plb::saveBinaryBlock(lattice_fluid2, Lattice2);

          std::string run_name = outDir + "/run_num.dat";
          plb::plb_ofstream ofile1(run_name.c_str());
          ofile1 << runs << std::endl;
          ofile1.close();

          std::string iter_name = outDir + "/iter_num.dat";
          plb::plb_ofstream ofile_iter(iter_name.c_str());
          ofile_iter << iT << std::endl;
          ofile_iter.close();
        }
      }

      // Save final state after convergence or max iterations per pressure step
      if (checkconv == 1) {
        writeGif_f1(lattice_fluid1, lattice_fluid2, runs_str, iT);

        // saves converged state vtks
        if (it_vtk < 100000) {
          writeVTK_rho(lattice_fluid1, "rho_f1_", iT);
          // Write x-component velocity data
          writeVTK_vel(lattice_fluid1, "vel_f1_", iT, false, true);

          if (rho_vtk == true) {
            writeVTK_rho(lattice_fluid2, "rho_f2_", iT);
          }
        }

        // saves a .dat file with the run number (for restarting sim)
        std::string run_name;
        run_name = outDir + "/run_num.dat";
        plb::plb_ofstream ofile1(run_name.c_str());
        ofile1 << runs + 1 << std::endl;
        ofile1.close();

        // saves a .dat file (lightweight) with the density
        std::string rho_name;
        rho_name = outDir + "/rho_f1_" + runs_str + ".dat";
        plb::plb_ofstream ofile2(rho_name.c_str());
        ofile2 << std::setprecision(2) << *plb::computeDensity(lattice_fluid1)
               << std::endl;
        ofile2.close();

        // saves a .dat file (lightweight) with the velocity
        std::string vel_name;
        vel_name = outDir + "/vel_f1_" + runs_str + ".dat";
        plb::plb_ofstream ofile3(vel_name.c_str());
        ofile3 << std::setprecision(1) << *plb::computeVelocity(lattice_fluid1)
               << std::endl;
        ofile3.close();

        // saves a binary file (heavy) with the sim state
        if (save_sim == true && iT > 0) {
          plb::pcout << "Saving restart files" << std::endl;
          plb::saveBinaryBlock(lattice_fluid1, Lattice1);
          plb::saveBinaryBlock(lattice_fluid2, Lattice2);

          std::string run_name = outDir + "/iter_num.dat";
          plb::plb_ofstream ofile1(run_name.c_str());
          ofile1 << 0 << std::endl;
          ofile1.close();

          // Need to save iteration number, save_it, in a file called
          // current_iteration.dat This way we can load stuff in the middle of a
          // pressure step or during steady state. Add this to the conditional
          // statement: && iT % save_it == 0 Uncomment save_it from inputs above
        }

        // Calculate and print velocity here for both fluids in x-direction
        computeVelocity_f1(lattice_fluid1, nu_f1);
        computeVelocity_f2(lattice_fluid2, nu_f2);
      }
    }
  }

  std::string output = outDir + "/output.dat";
  t = clock() - t;
  plb::pcout << "Simulation took seconds:" << ((float)t) / CLOCKS_PER_SEC
             << std::endl;
  plb::plb_ofstream ofile(output.c_str());
  ofile << "Output of the Simulation Run" << "\n\n";
  ofile << "Simulation took seconds =" << ((float)t) / CLOCKS_PER_SEC << "\n"
        << std::endl;

  ofile << "Kinematic viscosity f1 = " << nu_f1 << "\n" << std::endl;
  ofile << "Kinematic viscosity f2 = " << nu_f2 << "\n" << std::endl;
  ofile << "Gads_f1_s1 = " << Gads_f1_s1 << "\n" << std::endl;
  ofile << "Gads_f1_s2 = " << Gads_f1_s2 << "\n" << std::endl;
  ofile << "Gc = " << G << "\n" << std::endl;
  ofile << "Dissolved density = " << rhoNoFluid << "\n" << std::endl;
  ofile << "Inlet density = " << rho_f1_inlet << "\n" << std::endl;
  ofile << "Geometry flow length = " << nx << "\n" << std::endl;

  for (plb::plint runs = 0; runs < runnum; ++runs) {

    plb::pcout << "Run    = " << runs << std::endl;
    plb::pcout << "Pressure difference =  " << deltaP[runs] << std::endl;

    ofile << "Run = " << runs << "\n" << std::endl;
    ofile << "Pressure difference = " << deltaP[runs] << "\n" << std::endl;
  }
  ofile.close();

  return 0;
}
