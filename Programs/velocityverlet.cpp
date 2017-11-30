#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt){
    double dthalf = dt*0.5;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }
    //std::cout << "initial " <<system.atoms()[1]->velocity;
    //std::cout <<"new step"<<std::endl;
    for(Atom *atom : system.atoms()) {
        //std::cout <<atom->position<<std::endl;

        atom->velocity += dthalf*atom->force/(atom->mass());
        atom->position += atom->velocity*dt;// /atom->mass();
    }
    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += dthalf*atom->force/(atom->mass());
    }
    //std::cout << "  final: " << system.atoms()[1]->velocity<< std::endl;
}
