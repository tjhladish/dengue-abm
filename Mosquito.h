// Mosquito.h
// A single infected mosquito.
// This just takes care of mosquito lifespan, infecting serotype
// Movement is handled by other classes (Community)
#ifndef __MOSQUITO_H
#define __MOSQUITO_H
#include "Parameters.h"

class Location;

class Mosquito
{
    public:
        Mosquito();
        Mosquito(Location *p, Serotype serotype, int nInfectedAtID);
        virtual ~Mosquito();
        int getID() { return _nID; }
        Location *getLocation() { return _pLocation; }
        void setLocation(Location *p) { _pLocation = p; }
        Location *getOriginLocation() { return _pOriginLocation; }
        int getAgeInfected() { return _nAgeInfected; }
        int getAgeDeath() { return _nAgeDeath; }
        bool isDead() { return _bDead; }
        Serotype getSerotype() { return _eSerotype; }


    protected:
        int _nID;                                                     // unique identifier
        Location *_pLocation;                                         // pointer to present location
        Location *_pOriginLocation;                                   // pointer to origin (where infected) location
        int _nAgeInfected;                                            // age when infected in days
        int _nAgeDeath;                                               // lifespan in days
        bool _bDead;                                                  // is dead?
        Serotype _eSerotype;                                               // infecting serotype
        int _nInfectedAtID;                                           // location ID where infected
        //  int _nInfectedByID;     // which person infected me?
        static int _nNextID;                                          // unique ID to assign to the next Mosquito allocated
        //  static double _fDeathProbability[MAXMOSQUITOAGE]; // probability of death each day
        static double _fAgeDistribution[MAXMOSQUITOAGE];            // cumulative density of mosquito ages (in days)
        //  static double _fLifespanDistribution[MAXMOSQUITOAGE];  // cumulative density of mosquito lifespans (in days)
};
#endif
