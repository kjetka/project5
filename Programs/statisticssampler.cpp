#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler(const char *filename){
    openFileStats(filename);     // spr: How is this  ofstream outfile; outfile.open(filename);

}

void StatisticsSampler::openFileStats(const char *filename) {
    if(m_file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    m_file.open(filename);
    headerToFile();
}


void StatisticsSampler::headerToFile(){
    m_file << "\t\t" << "timesteps" <<
              "\t\t" << "t" <<
              "\t\t" << "Temperature" <<
              "\t\t" << "Kin" <<
              "\t\t" << "Pot" <<
              "\t\t" << "TotalE" <<
              "\t\t" << "MSD" << endl;
}

void StatisticsSampler::closeFile(){
    if(m_file.is_open()) {
        m_file.close();
    }

}

void StatisticsSampler::saveToFile(System &system){
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        //m_file.open("../results/statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }

    m_file <<  "\t\t" << system.steps() <<
               "\t\t" << system.time() <<
               "\t\t" << temperature() <<
               "\t\t" << kineticEnergy() <<
               "\t\t" << potentialEnergy() <<
               "\t\t" << totalEnergy() <<
               "\t\t" << MSD() << "\n";;//\n faster than endl;
}

void StatisticsSampler::sample(System &system){
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusionConst(system);
    testEnergyConservation();
    //sampleDiffusionConst(system);
}
void StatisticsSampler::testEnergyConservation(){
    double conservecriteria = 2e-4;
    if (m_totEnergyPreviousStep!=0){
        if((totalEnergy()/ (double) m_totEnergyPreviousStep-1) > conservecriteria){
            std::cout<<   "ERROR: Energy is not conserved! check StatisticsSampler class" <<std::endl;
            std::cout<< "  current/previous energy - 1 ="<< totalEnergy()/ m_totEnergyPreviousStep -1 <<endl;
            std::cout <<  "  current/previous energy -1 must be less than " <<conservecriteria <<endl;
            std::cout<< "  current energy = "<< totalEnergy()<<std::endl;
            std::cout<< " m_totEnergyPreviousStep = "<<m_totEnergyPreviousStep <<endl;

           exit(EXIT_FAILURE);

        }
    }
    m_totEnergyPreviousStep = totalEnergy();
}


void StatisticsSampler::sampleKineticEnergy(System &system){
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system){
    m_potentialEnergy = system.potential().potentialEnergy();
}


void StatisticsSampler::sampleTemperature(System &system){
    //Instantanous temperature
    m_temperature = (2.0/3.0)*m_kineticEnergy/(double)(system.nrAtoms());
}

void StatisticsSampler::sampleDensity(System &system){
    double V = system.volume();
    int atoms = system.atoms().size();
    m_density = atoms/V;
}

void StatisticsSampler::sampleMSD(System &system){
    double sum_SquaredDisplacement;
    for(Atom *atom : system.atoms()) {
        atom->Displacement=((atom->position-atom->boundaryJumps) - atom->position0).length();
        sum_SquaredDisplacement += atom->Displacement*atom->Displacement;
    }
    m_MSD = sum_SquaredDisplacement/system.nrAtoms();
}

void StatisticsSampler::sampleDiffusionConst(System &system){
    sampleMSD(system);

}



