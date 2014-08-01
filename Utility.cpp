#include "Utility.h"

namespace dengue {
    namespace util {
        Fit* lin_reg(const std::vector<double> &x, const std::vector<double> &y) {
            assert( x.size() == y.size() );
            Fit* fit = new Fit();
            const int n = x.size();
            double sumx = 0.0;                        // sum of x                      
            double sumx2 = 0.0;                       // sum of x**2                   
            double sumxy = 0.0;                       // sum of x * y                  
            double sumy = 0.0;                        // sum of y                      
            double sumy2 = 0.0;                       // sum of y**2                   

            for (int i=0; i<n; i++)   { 
                sumx  += x[i];       
                sumx2 += pow(x[i],2);  
                sumxy += x[i] * y[i];
                sumy  += y[i];      
                sumy2 += pow(y[i],2); 
            } 

            double denom = n * sumx2 - pow(sumx,2);
            if (denom == 0) {
                // singular matrix. can't solve the problem.
                fit->m   = 0;
                fit->b   = 0;
                fit->rsq = 0;
                return fit;
            }

            fit->m = (n * sumxy  -  sumx * sumy) / denom;
            fit->b = (sumy * sumx2  -  sumx * sumxy) / denom;
            // compute correlation coeff
            fit->rsq = pow((sumxy - sumx * sumy / n) / sqrt((sumx2 - pow(sumx,2)/n) * (sumy2 - pow(sumy,2)/n)),2);

            return fit; 
        }

        vector<string> split(const string &s, char delim) {
            vector<string> tokens;
            stringstream ss(s);
            string token;
            while (getline(ss, token, delim)) {
                tokens.push_back(token);
            }
            return tokens;
        }
    }
}
