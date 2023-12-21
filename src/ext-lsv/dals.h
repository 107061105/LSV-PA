/**
 * @file dals.h
 * @brief
 * @author Nathan Zhou
 * @date 2019-01-20
 * @bug No known bugs.
 */

#ifndef DALS_DALS_H
#define DALS_DALS_H

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "bdd/extrab/extraBdd.h"
#include "sat/cnf/cnf.h"
#include <map>
#include <memory>
#include <vector>
#include <unordered_map>
// #include "network.h"
// #include "framework.h"
// #include "sta.h"
// #include "dinic.h"
#include "abc_plus.h"

// using namespace abc_plus;
using NtkPtr = Abc_Ntk_t *;
using ObjPtr = Abc_Obj_t *;
/////////////////////////////////////////////////////////////////////////////
/// Class ALC, Approximate Local Change
/////////////////////////////////////////////////////////////////////////////

class ALC {
public:
    //---------------------------------------------------------------------------
    // Getters & Setters
    //---------------------------------------------------------------------------
    double GetError() const;

    Abc_Obj_t * GetTarget() const;

    Abc_Obj_t * GetSubstitute() const;

    bool IsComplemented() const;

    void SetError(double err);

    void SetTarget(Abc_Obj_t * t);

    void SetSubstitute(Abc_Obj_t * sub);

    void SetComplemented(bool is_complemented);

    //---------------------------------------------------------------------------
    // ALC Methods
    //---------------------------------------------------------------------------
    void Do();

    void Recover();

    //---------------------------------------------------------------------------
    // Operator Methods, Constructors & Destructors
    //---------------------------------------------------------------------------
//    ALC &operator=(const ALC &other);

    ALC();

    ALC(Abc_Obj_t * t, Abc_Obj_t * s, bool is_complemented, double error = 1);

    ~ALC();

private:
    double error_;
    bool is_complemented_;
    Abc_Obj_t * target_;
    Abc_Obj_t * substitute_;
    Abc_Obj_t * inv_;
    std::vector<Abc_Obj_t *> target_fan_outs_;
    std::unordered_map<Abc_Obj_t *, std::vector<Abc_Obj_t *>> target_fan_out_fan_ins_;
};

/////////////////////////////////////////////////////////////////////////////
/// Singleton Class DALS, Delay-Driven Approximate Logic Synthesis
/////////////////////////////////////////////////////////////////////////////


class DALS {
public:
    //---------------------------------------------------------------------------
    // Getters & Setters
    //---------------------------------------------------------------------------
    static std::shared_ptr<DALS> GetDALS();

    Abc_Ntk_t * GetApproxNtk();

    void SetTargetNtk(Abc_Ntk_t * ntk);

    void SetSim64Cycles(int sim_64_cycles);

    //---------------------------------------------------------------------------
    // DALS Methods
    //---------------------------------------------------------------------------
    void CalcTruthVec(bool show_progress_bar = false);

    void CalcALCs(const std::vector<Abc_Obj_t *> &target_nodes, bool show_progress = false, int top_k = 3);

    double EstSubPairError(Abc_Obj_t * target, Abc_Obj_t * substitute);

    void Run(double err_constraint = 0.15);

    //---------------------------------------------------------------------------
    // Operator Methods, Constructors & Destructors
    //---------------------------------------------------------------------------
    void operator=(DALS const &) = delete;

    DALS(abc_plus::Framework const &) = delete;

    ~DALS();

private:
    Abc_Ntk_t * target_ntk_;
    Abc_Ntk_t * approx_ntk_;
    int sim_64_cycles_;
    std::unordered_map<Abc_Obj_t *, std::vector<uint64_t>> truth_vec_;
    std::unordered_map<Abc_Obj_t *, std::vector<ALC>> cand_alcs_;
    std::unordered_map<Abc_Obj_t *, ALC> opt_alc_;

    DALS();
};

#endif
