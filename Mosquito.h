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
        Mosquito(Location* p, Serotype s, int nInfectedAtID, int nExternalIncubationPeriod);
        Mosquito(Location* p, Serotype s, int ageInfd, int ageInfs, int ageDead);

        virtual ~Mosquito();
        int getID() const { return _nID; }
        Location* getLocation() const { return _pLocation; }
        void setLocation(Location *p) { _pLocation = p; }
        void updateLocation(Location *p) { _pLocation->removeInfectedMosquito(); setLocation(p); p->addInfectedMosquito(); }
        Location* getOriginLocation() const { return _pOriginLocation; }
        int getAgeInfected() const { return _nAgeInfected; }
        int getAgeInfectious() const { return _nAgeInfectious; }
        int getAgeDeath() const { return _nAgeDeath; }
        bool isDead() const { return _bDead; }
        Serotype getSerotype() const { return _eSerotype; }


    protected:
        int _nID;                                                     // unique identifier
        Location* _pLocation;                                         // pointer to present location
        Location* _pOriginLocation;                                   // pointer to origin (where infected) location
        Serotype _eSerotype;                                          // infecting serotype
        int _nAgeInfected;                                            // age when infected in days
        int _nAgeInfectious;                                          // age when infectious in days
        int _nAgeDeath;                                               // lifespan in days
        bool _bDead;                                                  // is dead?
        int _nInfectedAtID;                                           // location ID where infected
        static int _nNextID;                                          // unique ID to assign to the next Mosquito allocated
};
#endif
