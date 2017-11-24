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
        atom->velocity += atom->force*dthalf/atom->mass();
        atom->position += atom->velocity*dt/atom->mass();
    }

    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*dthalf/atom->mass();
    }
}

