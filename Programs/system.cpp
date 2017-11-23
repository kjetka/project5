#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System(int nrAtoms_){
 nrAtoms = nrAtoms_;
}

// spr ???
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
            if (atom->position[i] < 0) atom->position[i] += m_systemSize[i];
            if (atom->position[i] < m_systemSize[i]) atom->position[i] -= m_systemSize[i];
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec3 totalMomentum;
    for(auto& atom : m_atoms){
        totalMomentum += atom->velocity*atom->mass();
    }
    vec3 everyMomentumChange = totalMomentum/nrAtoms;
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
    double almost0 = 1e-13;
    if (totalMomentumTest.length() > almost0){
        std::cout<<   "ERROR: " <<std::endl;
        std::cout<< "   length of totalMomentum greater than "<< almost0<< " after System::removeTotalMomentum was called"    <<std::endl;
        exit(EXIT_FAILURE);
    }

}



void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for(int i=0; i<nrAtoms; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
}




void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt); //linking to velocity verlet
    m_steps++;
    m_time += dt;
}
