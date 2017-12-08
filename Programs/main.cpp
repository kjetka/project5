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
    int nrUnitCellsEachDirection =5;
    int timeLimit = 2e5;
    vector<double> Temperatures_si = {590, 595, 600, 605, 610, 615};
    //vector<double> Temperatures_si = {100.0};

    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26);
    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.
    int printrate = timeLimit/(double) 1e3;


/*
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

*/
    cout << "discussion: better to have kinetic energy in Lennard Jones class? Atom class?"<<endl;
    cout << "check if applyPeriodicBoundaryConditions works for diffusion"<< endl;

    cout << "writing to file " << timeLimit/printrate << " times "<<endl;

    cout << UnitConverter::temperatureToSI(2.0)<<endl;

    cout << UnitConverter::energyFromSI(119.8*UnitConverter::kb)<<endl;

    for(int temperature_current:Temperatures_si){

        double initialTemperature = UnitConverter::temperatureFromSI(temperature_current); //Kelvin
        cout << "MD temp: "<<initialTemperature<< " SI temp: "<< temperature_current<<endl;
///*
//        cout << "------------------------------------------------"<<endl;
//        cout << setw(20) << "Timestep" <<
//                setw(20) << "Time" <<
//                setw(20) << "Temperature" <<
//                setw(20) << "KineticEnergy" <<
//                setw(20) << "PotentialEnergy" <<
//                setw(20) << "TotalEnergy"  << endl;
//*/


//        // setting up system
        System system(nrUnitCellsEachDirection);
        system.createFCCLattice(latticeConstant, initialTemperature);
        system.removeTotalMomentum();

//        // input parametres for force and pot.
        system.potential().setEpsilon(UnitConverter::energyFromSI(119.8*UnitConverter::kb));
        system.potential().setSigma(UnitConverter::lengthToAngstroms(3.405));

        // preparing output files
        string movietitle = "../results/movies2/movie_long_T_"+to_string(temperature_current)+".xyz";
        IO movie(movietitle.c_str());
        string txtfilename = "../results/txt2/nearTc_long_T_"+to_string(temperature_current)+".txt";
        StatisticsSampler statisticsSampler(txtfilename.c_str());


//        // integration loop
        system.calculateForces();   // in order to sample both kin and pot energy at t=0
        for(int timestep=0; timestep<timeLimit; timestep++) {

            statisticsSampler.sample(system); // system - same as *this within a object.
//            //write  to file (and print)
            if( timestep % printrate == 0||timestep ==0 ) {
///*              cout << setw(20) << system.steps()<<
//                      setw(20) << system.time() <<
//                      setw(20) << statisticsSampler.temperature() <<
//                      setw(20) << statisticsSampler.kineticEnergy() <<
//                      setw(20) << statisticsSampler.potentialEnergy() <<
//                      setw(20) << statisticsSampler.totalEnergy() << endl;
//*/

                statisticsSampler.saveToFile(system);
                movie.saveState(system);
            }

           system.step(dt);

        } // End integration loop

        movie.close();
        statisticsSampler.closeFile();



    } // end Temperature loop :)











    return 0;
}
