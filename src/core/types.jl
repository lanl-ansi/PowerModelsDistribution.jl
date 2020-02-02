""
abstract type AbstractNLPUBFModel <: _PMs.AbstractBFQPModel end


""
abstract type AbstractConicUBFModel <: _PMs.AbstractBFConicModel end


AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFModel <: AbstractConicUBFModel end


"SDP BFM with KCL as matrix equation, Geth 2020 (under review)"
abstract type SDPUBFKCLMXModel <: SDPUBFModel end


KCLMXModels = Union{SDPUBFKCLMXModel}


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFModel <: AbstractNLPUBFModel end


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFModel <: AbstractConicUBFModel end


SOCUBFModels = Union{SOCNLPUBFModel, SOCConicUBFModel}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFModel <: AbstractNLPUBFModel end


"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPUBFDiagModel <: AbstractLPUBFModel end
const LinDist3FlowModel = LPUBFDiagModel # more popular name for it


"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
# This formulation is equivalent to LPUBFDiagModel, with redundant variables and constraints
const LPUBFFullModel = LPUBFDiagModel


"default SDP unbalanced DistFlow constructor"
mutable struct SDPUBFPowerModel <: SDPUBFModel _PMs.@pm_fields end


"default SDP unbalanced DistFlow with matrix KCL constructor"
mutable struct SDPUBFKCLMXPowerModel <: SDPUBFKCLMXModel _PMs.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCNLPUBFPowerModel <: SOCNLPUBFModel _PMs.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCConicUBFPowerModel <: SOCConicUBFModel _PMs.@pm_fields end


"default LP unbalanced DistFlow constructor"
mutable struct LPUBFDiagPowerModel <: LPUBFDiagModel _PMs.@pm_fields end
const LinDist3FlowPowerModel = LPUBFDiagPowerModel # more popular name


"Legacy name for LPUBFDiagModel; print a warning that this will be dropped soon."
const LPUBFFullPowerModel = LPUBFDiagPowerModel


"This is completely equivalent to LPUBFDiagModel."
#const LPLinUBFPowerModel = LPdiagUBFPowerModel
mutable struct LPLinUBFPowerModel <: LPUBFDiagModel _PMs.@pm_fields end
