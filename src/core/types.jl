""
abstract type AbstractNLPUBFModel <: _PMs.AbstractBFQPModel end


""
abstract type AbstractConicUBFModel <: _PMs.AbstractBFConicModel end


AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFModel <: AbstractConicUBFModel end


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFModel <: AbstractNLPUBFModel end


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFModel <: AbstractConicUBFModel end

SOCUBFModels = Union{SOCNLPUBFModel, SOCConicUBFModel}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFModel <: AbstractNLPUBFModel end


"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
abstract type LPfullUBFModel <: AbstractLPUBFModel end


"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPdiagUBFModel <: AbstractLPUBFModel end


"default SDP unbalanced DistFlow constructor"
mutable struct SDPUBFPowerModel <: SDPUBFModel _PMs.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCNLPUBFPowerModel <: SOCNLPUBFModel _PMs.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCConicUBFPowerModel <: SOCConicUBFModel _PMs.@pm_fields end


"default LP unbalanced DistFlow constructor"
mutable struct LPfullUBFPowerModel <: LPfullUBFModel _PMs.@pm_fields end


"default LP unbalanced DistFlow constructor"
mutable struct LPdiagUBFPowerModel <: LPdiagUBFModel _PMs.@pm_fields end


"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current in scalar form"
abstract type LPLinUBFModel <: _PMs.AbstractBFModel end


"default Lin3Distflow constructor for scalar form"
mutable struct LPLinUBFPowerModel <: LPLinUBFModel _PMs.@pm_fields end
