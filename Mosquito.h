// Mosquito.h
// A single infected mosquito.
// This just takes care of mosquito lifespan, infecting serotype
// Movement is handled by other classes (Community)
#ifndef __MOSQUITO_H
#define __MOSQUITO_H
#include "Parameters.h"
#include "Location.h" 

class Location;
//Location::removeInfectedMosquito();
//Location::addInfectedMosquito();

class Mosquito
{
    public:
        Mosquito();
        Mosquito(Location *p, Serotype serotype, int nInfectedAtID, int nExternalIncubationPeriod);
        virtual ~Mosquito();
        int getID() { return _nID; }
        Location *getLocation() { return _pLocation; }
        void setLocation(Location *p) { _pLocation = p; }
        void updateLocation(Location *p) { _pLocation->removeInfectedMosquito(); setLocation(p); _pLocation->addInfectedMosquito(); }
        Location *getOriginLocation() { return _pOriginLocation; }
        int getAgeInfected() { return _nAgeInfected; }
        int getAgeInfectious() { return _nAgeInfectious; }
        int getAgeDeath() { return _nAgeDeath; }
        bool isDead() { return _bDead; }
        Serotype getSerotype() { return _eSerotype; }


    protected:
        int _nID;                                                     // unique identifier
        Location *_pLocation;                                         // pointer to present location
        Location *_pOriginLocation;                                   // pointer to origin (where infected) location
        int _nAgeInfected;                                            // age when infected in days
        int _nAgeInfectious;                                          // age when infectious in days
        int _nAgeDeath;                                               // lifespan in days
        bool _bDead;                                                  // is dead?
        Serotype _eSerotype;                                          // infecting serotype
        int _nInfectedAtID;                                           // location ID where infected
        static int _nNextID;                                          // unique ID to assign to the next Mosquito allocated
};
#endif
