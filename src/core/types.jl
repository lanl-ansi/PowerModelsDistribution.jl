"Supported data model types"
@enum DataModel ENGINEERING MATHEMATICAL DSS MATPOWER

"Load Models"
@enum LoadModel POWER CURRENT IMPEDANCE EXPONENTIAL ZIP

"Shunt Models"
@enum ShuntModel GENERIC CAPACITOR REACTOR

"Switch States"
@enum SwitchState OPEN CLOSED

"Generator, Solar, Storage, Wind Control Modes"
@enum ControlMode FREQUENCYDROOP ISOCHRONOUS

"Configurations"
@enum ConnConfig WYE DELTA

"Dispatchable"
@enum Dispatchable NO YES

"Status"
@enum Status DISABLED ENABLED

PowerModelsDistributionEnums = Union{DataModel,LoadModel,ShuntModel,SwitchState,ControlMode,ConnConfig,Dispatchable,Status}

"Base Abstract NLP Unbalanced Branch Flow Model"
abstract type AbstractNLPUBFModel <: _PM.AbstractBFQPModel end


"Base Abstract Conic Unbalanced Branch Flow Model"
abstract type AbstractConicUBFModel <: _PM.AbstractBFConicModel end


"Collection of Unbalanced Branch Flow Models"
AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFModel <: AbstractConicUBFModel end


"SDP BFM with KCL as matrix equation, Geth 2020 (under review)"
abstract type SDPUBFKCLMXModel <: SDPUBFModel end


"Collection of Semidefinite Models"  # TODO Better documentation, name?
KCLMXModels = Union{SDPUBFKCLMXModel}


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFModel <: AbstractNLPUBFModel end


"SOC relaxation of SDPUBFModel per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFModel <: AbstractConicUBFModel end


"Collection of Second Order Cone Models"
SOCUBFModels = Union{SOCNLPUBFModel, SOCConicUBFModel}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFModel <: AbstractNLPUBFModel end


"""
LinDist3Flow per Arnold et al. (2016), using vector variables for power, voltage and current

D. B. Arnold, M. Sankur, R. Dobbe, K. Brady, D. S. Callaway and A. Von Meier, "Optimal dispatch of reactive power for voltage regulation and balancing in unbalanced distribution systems," 2016 IEEE Power and Energy Society General Meeting (PESGM), Boston, MA, 2016, pp. 1-5, doi: 10.1109/PESGM.2016.7741261.
"""
abstract type LPUBFDiagModel <: AbstractLPUBFModel end
const LinDist3FlowModel = LPUBFDiagModel # more popular name for it


"default SDP unbalanced DistFlow constructor"
mutable struct SDPUBFPowerModel <: SDPUBFModel _PM.@pm_fields end


"default SDP unbalanced DistFlow with matrix KCL constructor"
mutable struct SDPUBFKCLMXPowerModel <: SDPUBFKCLMXModel _PM.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCNLPUBFPowerModel <: SOCNLPUBFModel _PM.@pm_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCConicUBFPowerModel <: SOCConicUBFModel _PM.@pm_fields end


"default LP unbalanced DistFlow constructor"
mutable struct LPUBFDiagPowerModel <: LPUBFDiagModel _PM.@pm_fields end
const LinDist3FlowPowerModel = LPUBFDiagPowerModel # more popular name
