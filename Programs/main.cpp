#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(){

    // Initial values setting up system
    int nrUnitCellsEachDirection =2;
    int timeLimit = 1e4;
    //vector<double> Temperatures_si = {50.0,85.0,300.0,500.0};
    //vector<double> Temperatures_si = {50.0,85.0};
    //vector<double> Temperatures_si = {300.0,500.0};
    vector<double> Temperatures_si = {300};

    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26);
    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.


/*
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

*/
    cout << "discussion: better to have kinetic energy in Lennard Jones class? Atom class?"<<endl;

    for(int temperature_current:Temperatures_si){

        double initialTemperature = UnitConverter::temperatureFromSI(temperature_current); // in Kelvin

        cout << "------------------------------------------------"<<endl;
        cout << "Md temp: "<<initialTemperature<< " Si temp: "<< temperature_current<<endl;
        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature" <<
                setw(20) << "KineticEnergy" <<
                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy"  << endl;



        // setting up system
        System system(nrUnitCellsEachDirection);

        system.createFCCLattice(latticeConstant, initialTemperature);

        system.potential().setEpsilon(UnitConverter::energyFromSI(119.8*UnitConverter::kb));
        system.potential().setSigma(UnitConverter::lengthToAngstroms(3.405));

        system.removeTotalMomentum();

        string movietitle = "../results/movies/movie_T_"+to_string(temperature_current)+".xyz";
        IO movie(movietitle.c_str());
        string txtfilename = "../results/txt/5d_T_"+to_string(temperature_current)+".txt";
        StatisticsSampler statisticsSampler(txtfilename.c_str());

       /* system.calculateForces();
        statisticsSampler.sample(system);
        cout << setw(20) << system.steps()<<
                setw(20) << system.time() <<
                setw(20) << statisticsSampler.temperature() <<
                setw(20) << statisticsSampler.kineticEnergy() <<
                setw(20) << statisticsSampler.potentialEnergy() <<
                setw(20) << statisticsSampler.totalEnergy() << endl;
        system.step(dt);
        statisticsSampler.sample(system);
        cout << setw(20) << system.steps()<<
                setw(20) << system.time() <<
                setw(20) << statisticsSampler.temperature() <<
                setw(20) << statisticsSampler.kineticEnergy() <<
                setw(20) << statisticsSampler.potentialEnergy() <<
                setw(20) << statisticsSampler.totalEnergy() << endl;
*/
        system.calculateForces();
        for(int timestep=0; timestep<timeLimit; timestep++) {

            statisticsSampler.sample(system); // system - same as *this within a object.

            if( timestep % 100 == 0||timestep ==0 ) { //approx 24*4=96 frames second.
              cout << setw(20) << system.steps()<<
                      setw(20) << system.time() <<
                      setw(20) << statisticsSampler.temperature() <<
                      setw(20) << statisticsSampler.kineticEnergy() <<
                      setw(20) << statisticsSampler.potentialEnergy() <<
                      setw(20) << statisticsSampler.totalEnergy() << endl;


                statisticsSampler.saveToFile(system);
                movie.saveState(system);
            }

           system.step(dt);

            // Defining the interval the variables are written to .txt and .xyz file
        } // End integration loop

        //cout << "check if applyPeriodicBoundaryConditions works for diffusion"<< endl;
        movie.close();
        statisticsSampler.closeFile();



    } // end Temperature loop :)












    return 0;
}
