#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "bdd/extrab/extraBdd.h"
#include "sat/cnf/cnf.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <string>
#include <unordered_map>

// #include <boost/filesystem.hpp>
// #include "abc_plus.h"
// #include "sta.h"
#include "dals.h"
// #include "playground.h"

// using namespace boost::filesystem;
// using namespace abc_plus;

static int Lsv_CommandLsvFinal(Abc_Frame_t* pAbc, int argc, char** argv);

void init(Abc_Frame_t* pAbc) {
    Cmd_CommandAdd(pAbc, "LSV", "lsv_final", Lsv_CommandLsvFinal, 0);
}

void destroy(Abc_Frame_t* pAbc) {}

Abc_FrameInitializer_t frame_initializer = { init, destroy };

struct PackageRegistrationManager {
    PackageRegistrationManager() { Abc_FrameAddInitializer(&frame_initializer); }
} lsvPackageRegistrationManager;

int Lsv_CommandLsvFinal(Abc_Frame_t* pAbc, int argc, char** argv)
{
    Abc_Ntk_t *pNtk = Abc_FrameReadNtk(pAbc);
    auto dals = DALS::GetDALS();
    std::cout << ":)" << std::endl;
    dals->SetTargetNtk(pNtk);
    dals->SetSim64Cycles(10000);
    dals->Run(0.15);

    auto approx_ntk = dals->GetApproxNtk();
    abc_plus::NtkWriteBlif(approx_ntk, "output.blif");

    return 0;
}