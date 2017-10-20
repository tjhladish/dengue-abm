#include <iostream>
#include "TestPar.h"

using namespace std;

int main() {
  for(auto p:MOSQUITO_DAILY_DEATH_PROBABILITY) cout << p << " ";

  cout << "\n\n" << "SURVIVAL" << "\n\n";

  for(auto p:MOSQUITO_DAILY_SURVIVAL_PROBABILITY) cout << p << " ";

  cout << "\n\n" << "SURVIVAL CUMULATIVE" << "\n\n";

  for(auto p:MOSQUITO_SURVIVE_CUMPROB) cout << p << " ";

  cout << "\n\n" << "RELATIVE AGE FRACTION" << "\n\n";

  for(auto p:MOSQUITO_AGE_RELFRAC) cout << p << " ";

  cout << "\n\n" << "AGE CDF" << "\n\n";

  for(auto p:MOSQUITO_AGE_CDF) cout << p << " ";

  cout << "\n\n" << "DEATH AGE CDF" << "\n\n";

  for(auto p:MOSQUITO_DEATHAGE_CDF) cout << p << " ";

  cout << "\n\n";

  return 0;
}
