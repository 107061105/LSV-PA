/**
 * @file dals.cpp
 * @brief
 * @author Nathan Zhou
 * @date 2019-01-20
 * @bug In some cases, the depth of the circuit stays the same for multiple rounds.
 */

#include <iostream>
#include <algorithm>
#include <bitset>
#include <iomanip>

// resolve conflict between cpu timers and original timers (deprecated) in boost library
// #define timer timer_deprecated


// #undef timer

#include "dals.h"
#include "sta.h"
#include "dinic.h"


/////////////////////////////////////////////////////////////////////////////
/// Class ALC, Approximate Local Change
/////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Getters & Setters
//---------------------------------------------------------------------------
double ALC::GetError() const { return error_; }

Abc_Obj_t * ALC::GetTarget() const { return target_; }

Abc_Obj_t * ALC::GetSubstitute() const { return substitute_; }

bool ALC::IsComplemented() const { return is_complemented_; }

void ALC::SetError(double err) { error_ = err; }

void ALC::SetTarget(Abc_Obj_t * t) { target_ = t; }

void ALC::SetSubstitute(Abc_Obj_t * sub) { substitute_ = sub; }

void ALC::SetComplemented(bool is_complemented) { is_complemented_ = is_complemented; }

//---------------------------------------------------------------------------
// ALC Methods
//---------------------------------------------------------------------------
void ALC::Do() {
    if (is_complemented_) {
        inv_ = ObjCreateInv(substitute_);
        ObjReplace(target_, inv_);
    } else
        ObjReplace(target_, substitute_);
}

void ALC::Recover() {

    if (is_complemented_)
        ObjDelete(inv_);
    for (auto const &fan_out : target_fan_outs_) {
        Abc_ObjRemoveFanins(fan_out);
        for (auto const &fan_in : target_fan_out_fan_ins_.at(fan_out))
            Abc_ObjAddFanin(fan_out, fan_in);
    }
}

//---------------------------------------------------------------------------
// Operator Methods, Constructors & Destructors
//---------------------------------------------------------------------------


ALC::ALC() = default;

ALC::ALC(Abc_Obj_t * t, Abc_Obj_t * s, bool is_complemented, double error) : error_(error), is_complemented_(is_complemented),
                                                                   target_(t), substitute_(s), inv_(nullptr),
                                                                   target_fan_outs_(std::vector<Abc_Obj_t *>()) {
    for (auto const &fan_out : abc_plus::ObjFanouts(target_)) {
        target_fan_outs_.push_back(fan_out);
        for (auto const &fan_in : abc_plus::ObjFanins(fan_out))
            target_fan_out_fan_ins_[fan_out].push_back(fan_in);
    }
}

ALC::~ALC() = default;

/////////////////////////////////////////////////////////////////////////////
/// Singleton Class DALS, Delay-Driven Approximate Logic Synthesis
/////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Getters & Setters
//---------------------------------------------------------------------------
std::shared_ptr<DALS> DALS::GetDALS() {
    static std::shared_ptr<DALS> dals(new DALS);
    return dals;
}

NtkPtr DALS::GetApproxNtk() { return approx_ntk_; }

void DALS::SetTargetNtk(NtkPtr ntk) {
    target_ntk_ = abc_plus::NtkDuplicate(ntk);
    approx_ntk_ = abc_plus::NtkDuplicate(target_ntk_);
}

void DALS::SetSim64Cycles(int sim_64_cycles) { sim_64_cycles_ = sim_64_cycles; }

//---------------------------------------------------------------------------
// DALS Methods
//---------------------------------------------------------------------------
void DALS::CalcTruthVec(bool show_progress_bar) { truth_vec_ = SimTruthVec(approx_ntk_, show_progress_bar, sim_64_cycles_); }

void DALS::CalcALCs(const std::vector<Abc_Obj_t *> &target_nodes, bool show_progress, int top_k) {
    cand_alcs_.clear();

    CalcTruthVec(show_progress);

    for (auto const &t_node : target_nodes)
        cand_alcs_.emplace(t_node, std::vector<ALC>());
    auto time_info = CalcSlack(approx_ntk_);
    auto s_nodes = NtkTopoSortPINode(approx_ntk_);

    // calculate candidate ALCs for each target node
    for (auto const &t_node : target_nodes) {
        // if (show_progress) ++(*pd);
        for (auto const &s_node : s_nodes)
            if (t_node != s_node && time_info.at(s_node).arrival_time < time_info.at(t_node).arrival_time) {
                double est_error = EstSubPairError(t_node, s_node);
                if (time_info.at(s_node).arrival_time < time_info.at(t_node).arrival_time - 1)
                    cand_alcs_[t_node].emplace_back(t_node, s_node, est_error > 0.5, std::min(est_error, 1 - est_error));
                else
                    cand_alcs_[t_node].emplace_back(t_node, s_node, false, est_error);
            }
        if (cand_alcs_[t_node].size() > top_k)
            std::partial_sort(cand_alcs_[t_node].begin(), cand_alcs_[t_node].begin() + top_k, cand_alcs_[t_node].end(),
                              [](const auto &a, const auto &b) {
                                  return a.GetError() < b.GetError();
                              });
        else
            std::sort(cand_alcs_[t_node].begin(), cand_alcs_[t_node].end(),
                      [](const auto &a, const auto &b) {
                          return a.GetError() < b.GetError();
                      });
    }
    // if (show_progress)
    //     std::cout << "Calc Candidate ALCs Finished" << timer.format() << std::endl;

    // calculate the most optimal ALC for each target node
    std::vector<ALC> k_alcs;
    for (auto const &t_node : target_nodes) {
        // if (show_progress) ++(*pd);
        k_alcs.clear();
        for (auto alc: cand_alcs_[t_node]) {
            if (k_alcs.size() == top_k) break;
            alc.Do();
            alc.SetError(SimER(target_ntk_, approx_ntk_));
            alc.Recover();
            k_alcs.push_back(alc);
        }
        std::sort(k_alcs.begin(), k_alcs.end(),
                  [](const auto &a, const auto &b) {
                      return a.GetError() < b.GetError();
                  });
        opt_alc_.emplace(t_node, k_alcs.front());
    }
    // if (show_progress)
    //     std::cout << "Calc Optimal ALC Finished" << timer.format() << std::endl;
}

double DALS::EstSubPairError(Abc_Obj_t * target, Abc_Obj_t * substitute) {
    int err_cnt = 0;
    for (int i = 0; i < sim_64_cycles_; i++)
        err_cnt += std::bitset<64>(truth_vec_.at(target)[i] ^ truth_vec_.at(substitute)[i]).count();
    return (double) err_cnt / (double) (64 * sim_64_cycles_);
}

void DALS::Run(double err_constraint) {
    std::cout << "h1" << std::endl;
    if (approx_ntk_ == NULL) std::cout<<"empty ckt"<<std::endl;
    double err = 0;
    int round = 0;
    while (err < err_constraint) {
        round++;
        auto time_info = CalcSlack(approx_ntk_);
        std::cout << "h2" << std::endl;

        std::vector<Abc_Obj_t *> pis_nodes_0, nodes_0;
        for (auto const &obj : NtkTopoSortPINode(approx_ntk_))
            if (time_info.at(obj).slack == 0) {
                pis_nodes_0.push_back(obj);
                if (ObjIsNode(obj)) nodes_0.push_back(obj);
            }

        CalcALCs(nodes_0, false, 3);
        std::cout << "h3" << std::endl;

        int N = Abc_NtkObjNumMax(approx_ntk_) + 1;
        int source = 0, sink = N - 1;
        Dinic dinic(N * 2);

        for (const auto &obj_0 : pis_nodes_0) {
            int u = ObjID(obj_0);
            if (ObjIsPI(obj_0))
                dinic.AddEdge(source, u, std::numeric_limits<double>::max());
            else {
                if (opt_alc_.at(obj_0).GetError() == 0)
                    dinic.AddEdge(u, u + N, std::numeric_limits<double>::min());
                else
                    dinic.AddEdge(u, u + N, opt_alc_.at(obj_0).GetError());
                if (ObjIsPONode(obj_0))
                    dinic.AddEdge(u + N, sink, std::numeric_limits<double>::max());
            }
        }

        for (auto &[u, vs] : GetCriticalGraph(approx_ntk_)) {
            for (auto &v : vs) {
                if (ObjIsPI(NtkObjbyID(approx_ntk_, u)))
                    dinic.AddEdge(u, v, std::numeric_limits<double>::max());
                else
                    dinic.AddEdge(u + N, v, std::numeric_limits<double>::max());
            }
        }

        std::cout << "---------------------------------------------------------------------------" << std::endl;
        std::cout << "> Round " << round << std::endl;
        std::cout << "---------------------------------------------------------------------------" << std::endl;
        std::cout << "MinCut: " << std::endl;
        for (const auto edge : dinic.MinCut(source, sink)) {
            auto obj = NtkObjbyID(approx_ntk_, edge.u);
            std::cout << ObjName(obj) << "--->" << ObjName(opt_alc_.at(obj).GetSubstitute())
                      << " : " << opt_alc_.at(obj).IsComplemented()
                      << " : " << opt_alc_.at(obj).GetError()
                      << std::endl;
            opt_alc_.at(obj).Do();
        }

        err = SimER(target_ntk_, approx_ntk_);
        std::cout << "Error Rate: " << err << std::endl;
        std::cout << "Delay: "
                  << GetKMostCriticalPaths(target_ntk_, 1)[0].max_delay << "--->"
                  << GetKMostCriticalPaths(approx_ntk_, 1)[0].max_delay << std::endl;
    }
}

//---------------------------------------------------------------------------
// Operator Methods, Constructors & Destructors
//---------------------------------------------------------------------------
DALS::~DALS() {
    Abc_NtkDelete(target_ntk_);
    Abc_NtkDelete(approx_ntk_);
}

DALS::DALS() = default;
