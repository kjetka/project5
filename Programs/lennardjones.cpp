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
    double epsilon24 = 24*m_epsilon;
    double epsilon4 = 4*m_epsilon;
    double sigmaDivR6 = 1.0;
    //int i =0; int j = 0;

    int i_length = system.atoms().size();

    for(int i=0; i<i_length; i++){
        Atom *atom_i = system.atoms()[i];
        for(int j=i+1;j<i_length; j++){
            Atom *atom_j = system.atoms()[j];
            vec3 r_vec = atom_i->position - atom_j->position;
            for(int u =0; u<3;u++){
                double L = system.systemSize()[u];
                if(r_vec[u] > L*0.5) r_vec[u] -=  L;
                if(r_vec[u] <=-L*0.5) r_vec[u] +=  L;
            }
            double r = r_vec.length();

            //if(r<= m_maxRadii){
            for(int gange=0;gange<6;gange++) { sigmaDivR6*=m_sigma/r;}
            atom_i->force +=  epsilon24*( 2*sigmaDivR6*sigmaDivR6 - sigmaDivR6  ) * r_vec/(r*r);
            m_potentialEnergy += epsilon4*(  sigmaDivR6*sigmaDivR6- sigmaDivR6   );
            //}
        }

    }

   /*




    for(auto& atom_i : system.atoms()){

        for(auto& atom_j : system.atoms()){
            vec3 r_vec = atom_i->position - atom_j->position;
            double r = r_vec.length();
            if(r!=0 || r<= m_maxRadii){

                for(int gange=0;gange<6;gange++) { sigmaDivRij6*=m_sigma/r;}
                atom_i->force +=  epsilon24*( 2*sigmaDivRij6*sigmaDivRij6 - sigmaDivRij6  ) * r_vec/(r*r);
            }
        }
    }


    double rij2 = rij*rij;
    double Fx = epsilon24*( 2*sigmaDivRij6*sigmaDivRij6 - sigmaDivRij6  ) * xij/(rij2);
    double Fy = epsilon24*( 2*sigmaDivRij6*sigmaDivRij6 - sigmaDivRij6  ) * yij/(rij2);
    double Fz = epsilon24*( 2*sigmaDivRij6*sigmaDivRij6 - sigmaDivRij6  ) * zij/(rij2);



    m_potentialEnergy = 0; // Remember to compute this in the loop
*/
}
