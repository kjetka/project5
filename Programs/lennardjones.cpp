#include "lennardjones.h"
#include "system.h"

double LennardJones::potentialEnergy() const{
    return m_potentialEnergy;
}

double LennardJones::sigma() const{
    return m_sigma;
}

void LennardJones::setSigma(double sigma){
    m_sigma = sigma;
    m_maxRadii = 5*sigma;
}

double LennardJones::epsilon() const{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon){
    m_epsilon = epsilon;
}


void LennardJones::calculateForces(System &system){

    m_potentialEnergy = 0;
    double epsilon24 = m_epsilon*24;
    double epsilon4 = m_epsilon*4;

    double L = system.systemSize()[0];
    //std::cout <<L<<std::endl;
    int i_length = system.atoms().size();

    for(int i=0; i<i_length; i++){
        Atom *atom_i = system.atoms()[i];
        for(int j=i+1;j<i_length; j++){
            Atom *atom_j = system.atoms()[j];
            vec3 r_vec = (atom_i->position - atom_j->position);
            for(int u =0; u<3;u++){
                if(r_vec[u] > L*0.5) r_vec[u] -=  L;
                if(r_vec[u] <-L*0.5) r_vec[u] +=  L;
            }
            double r = r_vec.length();
            //if (j==1) std::cout <<r << atom_i->position<<std::endl;
            double temp = m_sigma/r;
            double sigmaDivR6 = 1.0;
            for(int gange=0;gange<6;gange++) { sigmaDivR6*=temp;}

            double sigmaDivR12 = sigmaDivR6*sigmaDivR6;
            vec3 force = epsilon24*( 2*sigmaDivR12 - sigmaDivR6  ) * r_vec/(r*r);
            atom_i->force +=  force;
            atom_j->force -=  force;
            m_potentialEnergy += epsilon4*(  sigmaDivR12- sigmaDivR6   );
        }
    }
}
