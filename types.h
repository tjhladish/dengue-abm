#ifndef __TYPES_H
#define __TYPES_H
#include <string>
#include <vector>
#include <map>

#ifdef _OPENMP
# include <omp.h>
struct MutexType {
    MutexType() { omp_init_lock(&lock); }
    ~MutexType() { omp_destroy_lock(&lock); }
    void Lock() { omp_set_lock(&lock); }
    void Unlock() { omp_unset_lock(&lock); }

    MutexType(const MutexType& ) { omp_init_lock(&lock); }
    MutexType& operator= (const MutexType& ) { return *this; }
    public:
    omp_lock_t lock;
};
#else
/* A dummy mutex that doesn't actually exclude anything,
 * but as there is no parallelism either, no worries. */
struct MutexType {
    void Lock() {}
    void Unlock() {}
};
#endif

/* An exception-safe scoped lock-keeper. */
struct ScopedLock {
    explicit ScopedLock(MutexType& m) : mut(m), locked(true) { mut.Lock(); }
    ~ScopedLock() { Unlock(); }
    void Unlock() { if(!locked) return; locked=false; mut.Unlock(); }
    void LockAgain() { if(locked) return; mut.Lock(); locked=true; }
    private:
    MutexType& mut;
    bool locked;
    private: // prevent copying the scoped lock.
    void operator=(const ScopedLock&);
    ScopedLock(const ScopedLock&);
};

template <class T>
class VectorMP {
    public:
        VectorMP() {};
        VectorMP(int x):_vec(x) { };
        //VectorMP(int x, T & const t ):_vec(x,t) { }  
        VectorMP(int x, const T & t):_vec(x,t) { }  

        void push_back_mp( T const x ) {
            ScopedLock lck(lock);
            _vec.push_back(x);
        }
        T& back_mp() {
            ScopedLock lck(lock);
            return _vec.back(); 
        }
        int size_mp() {
            ScopedLock lck(lock);
            return (int) _vec.size();
        }
        void clear_mp() {
            ScopedLock lck(lock);
            _vec.clear();
        }
        T& operator[] (const int x) {
            ScopedLock lck(lock);
            return _vec[x];
        }

    private:
        std::vector<T*> _vec;
        MutexType lock;
};

    //map< Location*, map<int,bool> >::iterator it;
    //for ( it=_isHot.begin_mp() ; it != _isHot.end_mp(); it++ ) (*it).second.clear_mp();
//template <class KT, class VT> class MapMP;
//template <class KT, class VT> class MapMP<KT, VT>::iterator;

template <class KT, class VT>
class MapMP {
    //friend class MapMP::iterator;

    public:
        int count_mp( KT x ) {
            ScopedLock lck(lock);
            return _map.count(x);
        }
        void erase_mp(KT key) {
            ScopedLock lck(lock);
            _map.erase(key);
        }
        void clear_mp() {
            ScopedLock lck(lock);
            _map.clear();
        }
        VT& operator[] (const KT key) {
            ScopedLock lck(lock);
            return _map[key];
        }
        typename std::map<KT, VT>::iterator begin() {
            return _map.begin();
        }
        typename std::map<KT, VT>::iterator end() {
            return _map.end();
        }

    private:
        std::map<KT, VT> _map;
        MutexType lock;
};

/*template <class KT, class VT>
class MapMP::iterator {

}*/

/*class LocationMap2D {
    private:
        std::map< MosquitoMap1D > _mosquitoMap2D;  // queue of infectious mosquitoes with n days
        MutexType lock;
    public:
        void push_back_mp( MosquitoMap1D mv ) {
            ScopedLock lck(lock);
            _mosquitoMap2D.push_back(mv);
        }
        int size_mp() {
            ScopedLock lck(lock);
            return (int) _mosquitoMap2D.size();
        }
}*/


//        static std::map< Location*, std::map<int, bool> > _isHot;

 //       std::vector< std::vector<Mosquito*> > _infectiousMosquitoQueue;  // queue of infectious mosquitoes with n days
                                                                         // left to live
  //      std::vector< std::vector<Mosquito*> > _exposedMosquitoQueue;  // queue of exposed mosquitoes with n days of latency left

/*#include <set>
class data {
    private:
        std::set<int> flags;
        MutexType lock;
    public:
        bool set_get(int c)
        {
            ScopedLock lck(lock); // locks the mutex

            if(flags.find(c) != flags.end()) return true; // was found
            flags.insert(c);
            return false; // was not found
        } // automatically releases the lock when lck goes out of scope.
};*/

#endif
