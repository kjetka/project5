#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"

class Atom
{
private:
    double m_mass;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    vec3 position0;
    vec3 boundaryJumps;
    double MSD=0;

    Atom(double mass);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);
    void resetBoundaryJumps();

    double mass() { return m_mass; }
    void setMass(double mass) { m_mass = mass; }
};
#endif
