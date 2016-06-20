#include <iostream>
#include "StructuredJacobian.hh"
#include <cstdlib>


using namespace std;

int main(int argc, char *argv[])
{

    for(int i = 1; i < 20; i++)
    {
        StructuredJacobian structuredJacobian(i);
        // Get the Sparsity Pattern from the Object.
        auto_ptr<Sparsity> sparsity = structuredJacobian.getSparsityPattern();
        cout << "For n = " << i << " " << "Size: " << sparsity->size() - 1 << " 8n-16= " << 8*i -16 << endl;
    }

    // TODO: The size can be replaced by 8n-16 also for n > 6


    return 0;
}
