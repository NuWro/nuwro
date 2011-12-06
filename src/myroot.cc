//--------------------------------------------------
#include "TROOT.h"
#include "TRint.h"


int main(int argc, char **argv)
{
    // Create interactive interface
    TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0);

    // Run interactive interface
    theApp->Run();

    return(0);
}
//--------------------------------------------------
