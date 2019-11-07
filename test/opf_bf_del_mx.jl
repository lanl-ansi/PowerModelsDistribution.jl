# @info "running load model relaxation and matrix KCL tests"
#
# @testset "test distflow formulations" begin
#     @testset "test linearised distflow opf_bf" begin
#         @testset "5-bus lplinubf opf_bf" begin
#             mp_data = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
#             pm = PMs.build_model(mp_data, PMD.SDPUBFPowerModel, PMD.post_mc_opf_bf, ref_extensions=[ref_add_arcs_trans!], multiconductor=true)
#             sol = PMs.optimize_model!(pm, scs_solver)
#             #@test sol["termination_status"] == PMs.OPTIMAL
#         end
# end

using PowerModelsDistribution
const PMD = PowerModelsDistribution

import Memento

import InfrastructureModels

import PowerModels
const PMs = PowerModels

# Suppress warnings during testing.
const TESTLOG = Memento.getlogger(PowerModels)
Memento.setlevel!(TESTLOG, "error")

import JuMP
import Ipopt
import Cbc
import Juniper
import SCS
import Mosek, MosekTools

import JSON

using Test
using LinearAlgebra

pms_path = joinpath(dirname(pathof(PowerModels)), "..")
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6)
cbc_solver = with_optimizer(Cbc.Optimizer, logLevel=0)
scs_solver = with_optimizer(SCS.Optimizer, max_iters=10000)
mosek_solver = with_optimizer(Mosek.Optimizer)
juniper_solver = with_optimizer(Juniper.Optimizer, nl_solver=with_optimizer(Ipopt.Optimizer, tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])

include("common.jl")

data = PMD.parse_file("test_mx.dss")
data["gen"]["1"]["cost"] *= 100

pm_nlp = PMs.build_model(data, PMs.ACPPowerModel, PMD.post_mc_opf; ref_extensions=[ref_add_arcs_trans!], multiconductor=true)
sol_nlp = PMs.optimize_model!(pm_nlp, ipopt_solver)

##
pm_sdp_wye_vec = PMs.build_model(data, PMD.SDPUBFPowerModel, PMD.post_mc_opf_bf; ref_extensions=[ref_add_arcs_trans!], multiconductor=true)
sol_sdp_wye_vec = PMs.optimize_model!(pm_sdp_wye_vec, mosek_solver)
sol_sdp_wye_vec["objective"]-sol_nlp["objective"]

##
pm_sdp_del_vec = PMs.build_model(data, PMD.SDPUBFPowerModel, PMD.post_mc_opf_bf_del; ref_extensions=[ref_add_arcs_trans!], multiconductor=true)
sol_sdp_del_vec = PMs.optimize_model!(pm_sdp_del_vec, mosek_solver)
sol_sdp_del_vec["objective"]-sol_nlp["objective"]
##
pm_sdp_del_mx = PMs.build_model(data, PMD.SDPUBFPowerModel, PMD.post_mc_opf_bf_del_mx; ref_extensions=[ref_add_arcs_trans!], multiconductor=true)
sol_sdp_del_mx = PMs.optimize_model!(pm_sdp_del_mx, mosek_solver)
sol_sdp_del_mx["objective"]-sol_nlp["objective"]
