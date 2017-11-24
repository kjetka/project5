#include "lennardjones.h"
#include "system.h"
#include <cmath>

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma){
    m_sigma = sigma;
}

double LennardJones::epsilon() const{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon){
    m_epsilon = epsilon;
}


void LennardJones::findMaxForceRadius(){
    double r = m_sigma;
    double Force = 24*m_epsilon*(2*pow(m_sigma/r, 12)-pow(m_sigma/r, 6))/r;
    double limit = 0.04;
    double maxRadii;
    while ((r >= m_sigma+0.2) && (abs(Force)>=limit)){
        r += 0.01;
        double sigmaDivRadii = m_sigma/r;
        Force = 24*m_epsilon*(2*pow(sigmaDivRadii, 12)-pow(sigmaDivRadii, 6))/r;
    }
    maxRadii = r;
    setMaxForceRadius(maxRadii);
    std::cout<< maxRadii<<std::endl;
}

void LennardJones::setMaxForceRadius(double maxRadii){
    m_maxRadii = maxRadii;
}

void LennardJones::calculateForces(System &system){
    m_potentialEnergy = 0; // Remember to compute this in the loop
}
