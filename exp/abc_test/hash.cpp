//unordered_set::hash_function
#include <iostream>
#include <string>
#include <unordered_set>

typedef std::unordered_set<std::string> stringset;

int main ()
{
    stringset myset;

    stringset::hasher fn = myset.hash_function();

    std::cout << "that: " << fn ("that") << std::endl;
    std::cout << "than: " << fn ("than") << std::endl;

    return 0;
}
