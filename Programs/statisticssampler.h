#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>

class System; // Promise the compiler that this is a class even though we haven't included system.h here

class StatisticsSampler{
private:
    std::ofstream m_file;
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_totEnergyPreviousStep = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_MSD = 0;
public:

    StatisticsSampler(const char *filename);
    void openFileStats(const char *filename);
    void headerToFile();
    void closeFile();
    void saveToFile(System &system);

    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleMSD(System &system);
    void sampleDiffusionConst(System &system);

    void testEnergyConservation();

    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double prevTotalEnergy(){return m_totEnergyPreviousStep;}
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double MSD() { return m_MSD; }
};
#endif
