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
              "\t\t" << "TotalE" << endl;
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
               "\t\t" << totalEnergy() << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    saveToFile(system);
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
    m_temperature = (2.0/3.0)*m_kineticEnergy/(double)(system.nrAtoms());
}

void StatisticsSampler::sampleDensity(System &system){
    double V = system.volume();
    int atoms = system.atoms().size();
    m_density = atoms/V;
}
