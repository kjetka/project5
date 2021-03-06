#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System(int numberOfUnitCellsEachDimension_){
 numberOfUnitCellsEachDimension = numberOfUnitCellsEachDimension_;
 int AtomsPerUnitCell = 4;
 m_nrAtoms = pow(numberOfUnitCellsEachDimension,3) * AtomsPerUnitCell;
}


System::~System(){
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(auto& atom : m_atoms){ //python: for i in list...

        for(int i=0;i<3; i++){
            if (atom->position[i] < 0){
                atom->position[i] += m_systemSize[i];
                atom->boundaryJumps[i] += m_systemSize[i];
            }
            if (atom->position[i] > m_systemSize[i]){
                atom->position[i] -= m_systemSize[i];
                atom->boundaryJumps[i] -= m_systemSize[i];
            }
        }
    }
}

void System::removeTotalMomentum(){
    vec3 totalMomentum;
    for(auto& atom : m_atoms){
        totalMomentum += atom->velocity*atom->mass();
    }

    vec3 everyMomentumChange = totalMomentum/m_nrAtoms;
    for(auto& atom : m_atoms){
        atom->velocity -= everyMomentumChange/atom->mass();
    }
    test_removeTotalMomentum();
}

void System::test_removeTotalMomentum(){
    vec3 totalMomentumTest;
    for(auto& atom : m_atoms){
        totalMomentumTest += atom->velocity*atom->mass();
    }
    double almost0 = 5e-12;
    if (totalMomentumTest.length() > almost0){
        std::cout<<   "ERROR: " <<std::endl;
        std::cout<< "   length of totalMomentum greater than "<< almost0<< " after System::removeTotalMomentum was called"    <<std::endl;
        std::cout<< "   length of totalMomentum is "<< totalMomentumTest.length()<<std::endl;
        exit(EXIT_FAILURE);
    }


}



void System::createFCCLattice(double latticeConstant, double temperature) {

    vec3 iHat(1,0,0); iHat *= latticeConstant;
    vec3 jHat(0,1,0); jHat *= latticeConstant;
    vec3 kHat(0,0,1); kHat *= latticeConstant;

    vec3 r1(0,0,0);
    vec3 r2 = iHat/2.0 + jHat/2.0 + kHat*0;
    vec3 r3 = iHat*0 + jHat/2.0 + kHat/2.0;
    vec3 r4 = iHat/2.0  + jHat*0 + kHat/2.0;

    std::vector<vec3> basisvectors = {r1,r2,r3,r4};
    int atomsInEachUnitscell = basisvectors.size();
    Random::randomSeed();
    Random::nextGaussian(1,0.5);

    for(int Nx =0;Nx < numberOfUnitCellsEachDimension; Nx++){
        for(int Ny =0;Ny < numberOfUnitCellsEachDimension; Ny++){
            for(int Nz =0;Nz < numberOfUnitCellsEachDimension; Nz++){

                for(int iUnitcell =0; iUnitcell < atomsInEachUnitscell; iUnitcell++){
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    vec3 whichAtomInUnitCell = basisvectors[iUnitcell];
                    double x = whichAtomInUnitCell[0] + Nx*latticeConstant;
                    double y = whichAtomInUnitCell[1] + Ny*latticeConstant;
                    double z = whichAtomInUnitCell[2] + Nz*latticeConstant;
                    atom->position.set(x,y,z);
                    atom->position0.set(x,y,z);
                    atom->resetVelocityMaxwellian(temperature);
                    atom->boundaryJumps.set(0,0,0);
                    m_atoms.push_back(atom);
                }
            }
        }
    }
    double L = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(L, L, L));
}



void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential_class.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt); //linking to velocity verlet
    m_steps++;
    m_time += dt;
}
