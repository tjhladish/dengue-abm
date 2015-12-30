#ifndef DENGUE_UTILITY_H
#define DENGUE_UTILITY_H

#include <cstdlib>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

namespace dengue {
    namespace util {
        vector<string> split(const string &s, char delim);

        inline vector<string> read_vector_file(string filename, char sep=' ') {
            ifstream myfile(filename.c_str());
            if (!myfile) {
                cerr << "ERROR: " << filename << " not found." << endl;
                exit(116);
            }

            vector<string> V;
            if (myfile.is_open()) {
                string line;
                while ( getline(myfile,line) ) {
                    vector<string> fields = split(line, sep);
                    if (fields.size() == 0) {
                        cerr << "ERROR: Found line with no values in file: " << filename << " at line: " << V.size() << endl;
                        exit(117);
                    } else {
                        V.push_back(fields[0]);
                    }
                }
            }
            return V;
        }

        inline vector<vector<string> > read_2D_vector_file(string filename, char sep=' ') {
            ifstream myfile(filename.c_str());
            if (!myfile) {
                cerr << "ERROR: " << filename << " not found." << endl;
                exit(118);
            }

            vector<vector<string> > M;
            if (myfile.is_open()) {
                string line;

                while ( getline(myfile,line) ) {
                    vector<string> fields = split(line, sep);

                    vector<string> row(fields.size());
                    for( unsigned int i=0; i < fields.size(); i++ ) {
                            row[i] = fields[i];
                    }
                    M.push_back(row);
                }
            }
            return M;
        }

        class Fit {
            public:
                double m;
                double b;
                double rsq;
        };

        Fit* lin_reg(const std::vector<double> &x, const std::vector<double> &y);

        template <typename T>
        inline void cerr_vector(vector<T> & my_vector, string sep = " ") {
            for (int i = 0; i < my_vector.size() - 1; i++ ) cerr << my_vector[i] << sep;
            cerr << my_vector.back();
        }

        template <typename T>
        inline void cout_vector(vector<T> & my_vector, string sep = " ") {
            for (int i = 0; i < my_vector.size() - 1; i++ ) cout << my_vector[i] << sep;
            cout << my_vector.back();
        }

        inline double string2double(const std::string& s){ std::istringstream i(s); double x = 0; i >> x; return x; }

        template <typename T> inline T sum(vector<T> list) { T sum=0; for (unsigned int i=0; i<list.size(); i++) sum += list[i]; return sum;}
        template <typename T> inline double mean(vector<T> list) { return (double) sum(list) / list.size(); }
        template <typename T> inline long double meanl(vector<T> list) { return (long double) sum(list) / list.size(); }

        template <typename T> inline
        double median(vector<T> L) { 
            sort(L.begin(), L.end());
            float idx = (L.size() - 1.0) * 0.5;
            return ( L[ (int) ceil(idx) ] + L[ (int) floor(idx) ] ) /2.0;
        }

        // five number summary (min, 1st quartile, median, 3rd quartile, max)
        template <typename T> inline 
        vector<double> fivenum(vector<T> L) {
            assert(L.size() > 2);
            vector<double> stats(5);
            sort(L.begin(), L.end());
            stats[0] = L[0];        // min
            stats[4] = L.back();    // max

            float idx1 = (L.size() -1) * 0.25;
            float idx2 = (L.size() -1) * 0.5;
            float idx3 = (L.size() -1) * 0.75;
            
            stats[1] = (L[ceil(idx1)] + L[floor(idx1)]) /2.0;
            stats[2] = (L[ceil(idx2)] + L[floor(idx2)]) /2.0;
            stats[3] = (L[ceil(idx3)] + L[floor(idx3)]) /2.0;
         
            return stats;
        }

        template <typename T> inline 
        T min_element(vector<T> list) {
            T element = list[0];
            for (unsigned int i = 0; i < list.size(); i++) {
                element = min(element, list[i]);
            }
            return element;
        }

        template <typename T> inline 
        T max_element(vector<T> list) {
            T element = list[0];
            for (unsigned int i = 0; i < list.size(); i++) {
                element = max(element, list[i]);
            }
            return element;
        }

        template <typename T> inline 
        T range(vector<T> list) {
            return max_element(list) - min_element(list);
        }

        template <typename T>
        vector<double> normalize_dist (vector<T> dist, T sum) {
            vector<double> normed(dist.size());
            for (unsigned int i = 0; i < dist.size(); i++) normed[i] = ((double) dist[i]) / sum;
            return normed;
        }

        template <typename T>
        vector<double> normalize_dist (vector<T> dist) {
            return normalize_dist(dist, sum(dist));
        }

        template <typename T>
        inline std::string to_string (const T& t) {
            std::stringstream ss;
            ss << t;
            return ss.str();
        }

        inline float to_float(const std::string& s){
            std::istringstream i(s);
            float x = 0;
            i >> x;
            return x;
        }

        inline double to_double(const std::string& s){
            std::istringstream i(s);
            double x = 0;
            i >> x;
            return x;
        }

        inline int to_int(const std::string& s){
            std::istringstream i(s);
            int x = 0;
            i >> x;
            return x;
        }

        template <typename T>
        double variance(vector<T> & numbers) {
            double x = mean(numbers);
            double var_num = 0;
            int N = numbers.size();
            if (N == 1) return 0;
            for (int i=0; i<N; i++) var_num += pow(numbers[i] - x, 2);
            double var = var_num/(N-1);
            return var;
        }

        template <typename T>
        double stdev(vector<T> & numbers) { return sqrt( variance(numbers) ); }

        template <typename T>
        long double stdevl(vector<T> & numbers) { return sqrt( variance(numbers) ); }

        template <typename T>
        inline int sign(const T& _a) { return (int)((_a)<0 ? (-1) : (1)); }

        template <typename T>
        inline T MIN(const T& _a, const T& _b) { return ((_a)<(_b)?(_a):(_b));}

        template <typename T>
        inline T MAX(const T& _a,const T& _b) { return ((_a)>(_b)?(_a):(_b));}

        template <typename T>
        inline void delete_element(vector<T> & my_vector, T element) {
            for(int i = 0; i < my_vector.size(); i++) {
                if (my_vector[i] == element) {
                    my_vector.erase(my_vector.begin() + i);
                    break;
                }
            }
        }

        inline vector<int> tabulate_vector( vector<int> & my_vector ) {
            vector<int> tabulated(max_element(my_vector) + 1, 0);
            for (unsigned int i = 0; i<my_vector.size(); i++) tabulated[my_vector[i]]++;
            return tabulated;
        }

        template <typename T>  
        vector<T> shuffle_periods(const gsl_rng* RNG, vector<T> & sequence, int period_length = 365) {  
            vector<int> periods(sequence.size()/period_length);  
            for (unsigned int i = 0; i<periods.size(); i++) periods[i] = i;  
            gsl_ran_shuffle (RNG, periods.data(), periods.size(), sizeof(int));  
         
            vector<T> new_seq(sequence.size());  
            for(unsigned int i=0; i<periods.size(); ++i) {  
                for (int j = 0; j<period_length; j++) new_seq[i*period_length + j] = sequence[periods[i]*period_length + j];  
            }  
         
            if (sequence.size()/period_length < ((float) sequence.size())/period_length) {  
                int r = gsl_rng_uniform_int(RNG, 1+ceil(sequence.size()/period_length));  
                int offset = periods.size()*period_length;  
                for (unsigned int i = 0; i<new_seq.size() - offset; ++i) {  
                    new_seq[offset + i] = sequence[r*period_length + i];  
                }  
            }  
            return new_seq;  
        }

        inline int parseLine(char line[]){
            int i = strlen(line);
            while (*line < '0' || *line > '9') line++;
            line[i-3] = '\0';
            i = atoi(line);
            return i;
        }

        inline int getRamUsage(){ //Note: this value is in KB!
            FILE* file = fopen("/proc/self/status", "r");
            int result = -1;
            char line[128];

            while (fgets(line, 128, file) != NULL){
                if (strncmp(line, "VmSize:", 7) == 0){
                    result = parseLine(line);
                    break;
                }
            }
            fclose(file);
            printf("%d", result);
            return result;
        }
    }
}
#endif
