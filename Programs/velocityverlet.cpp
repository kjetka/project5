#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt){
    double dthalf = dt*0.5;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }
    for(Atom *atom : system.atoms()) {
        atom->velocity += dthalf*atom->force/(atom->mass());
        atom->position += atom->velocity*dt;// /atom->mass(); // This comment should not be here!!!!
    }
    //system.applyPeriodicBoundaryConditions();
    system.calculateForces();

    for(Atom *atom : system.atoms()) {
        atom->velocity += dthalf*atom->force/(atom->mass());
    }
}

