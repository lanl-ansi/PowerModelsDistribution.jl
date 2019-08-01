""
abstract type AbstractNLPUBFForm <: _PMs.AbstractBFQPForm end


""
abstract type AbstractConicUBFForm <: _PMs.AbstractBFConicForm end


AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFForm <: AbstractConicUBFForm end


"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFForm <: AbstractNLPUBFForm end


"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFForm <: AbstractConicUBFForm end

SOCUBFForm = Union{SOCNLPUBFForm, SOCConicUBFForm}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFForm <: AbstractNLPUBFForm end


"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
abstract type LPfullUBFForm <: AbstractLPUBFForm end


"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPdiagUBFForm <: AbstractLPUBFForm end


""
const SDPUBFPowerModel = _PMs.GenericPowerModel{SDPUBFForm}


"default SDP unbalanced DistFlow constructor"
SDPUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, SDPUBFForm; kwargs...)


""
const SOCNLPUBFPowerModel = _PMs.GenericPowerModel{SOCNLPUBFForm}


"default SOC unbalanced DistFlow constructor"
SOCNLPUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, SOCNLPUBFForm; kwargs...)


""
const SOCConicUBFPowerModel = _PMs.GenericPowerModel{SOCConicUBFForm}


"default SOC unbalanced DistFlow constructor"
SOCConicUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, SOCConicUBFForm; kwargs...)


""
const LPfullUBFPowerModel = _PMs.GenericPowerModel{LPfullUBFForm}


"default LP unbalanced DistFlow constructor"
LPfullUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, LPfullUBFForm; kwargs...)


""
const LPdiagUBFPowerModel = _PMs.GenericPowerModel{LPdiagUBFForm}


"default LP unbalanced DistFlow constructor"
LPdiagUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, LPdiagUBFForm; kwargs...)


"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current in scalar form"
abstract type LPLinUBFForm <: _PMs.AbstractBFForm end


""
const LPLinUBFPowerModel = _PMs.GenericPowerModel{LPLinUBFForm}


"default Lin3Distflow constructor for scalar form"
LPLinUBFPowerModel(data::Dict{String,Any}; kwargs...) = _PMs.GenericPowerModel(data, LPLinUBFForm; kwargs...)
