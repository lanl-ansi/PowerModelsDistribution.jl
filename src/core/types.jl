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

#================================================
    # exact non-convex models
    ACPUPowerModel, ACRUPowerModel, IVRUPowerModel
    # linear approximations
    DCPUPowerModel, DCMPPowerModel, BFAPowerModel, NFAUPowerModel
    # quadratic approximations
    DCPLLPowerModel, LPACCPowerModel
    # quadratic relaxations
    SOCWRPowerModel, SOCWRConicPowerModel,
    SOCBFPowerModel, SOCBFConicPowerModel,
    QCRMPowerModel, QCLSPowerModel,
    # sdp relaxations
    SDPWRMPowerModel, SparseSDPWRMPowerModel
================================================#

##### Top Level Abstract Types #####

"active power only models"
abstract type AbstractUnbalancedActivePowerModel <: AbstractUnbalancedPowerModel end

"variants that target conic solvers"
abstract type AbstractUnbalancedConicModel <: AbstractUnbalancedPowerModel end

"for branch flow models"
abstract type AbstractUBFModel <: AbstractUnbalancedPowerModel end


"for variants of branch flow models that target LP solvers"
abstract type AbstractUBFAModel <: AbstractUBFModel end

"for variants of branch flow models that target QP or NLP solvers"
abstract type AbstractUBFQPModel <: AbstractUBFModel end

"for variants of branch flow models that target conic solvers"
abstract type AbstractUBFConicModel <: AbstractUBFModel end

""
abstract type AbstractUnbalancedACPModel <: AbstractUnbalancedPowerModel end

"""
AC power flow Model with polar bus voltage variables.
The seminal reference of AC OPF:
```
@article{carpentier1962contribution,
  title={Contribution to the economic dispatch problem},
  author={Carpentier, J},
  journal={Bulletin de la Societe Francoise des Electriciens},
  volume={3},
  number={8},
  pages={431--447},
  year={1962}
}
```
History and discussion:
```
@techreport{Cain2012,
  author = {Cain, Mary B and {O' Neill}, Richard P and Castillo, Anya},
  title = {{History of optimal power flow and Models}},
  year = {2012}
  pages = {1--36},
  url = {https://www.ferc.gov/industries/electric/indus-act/market-planning/opf-papers/acopf-1-history-Model-testing.pdf}
}
```
"""
mutable struct ACPUPowerModel <: AbstractUnbalancedACPModel @pmd_fields end

""
abstract type AbstractUnbalancedACRModel <: AbstractUnbalancedPowerModel end


"""
AC power flow Model with rectangular bus voltage variables.
```
@techreport{Cain2012,
  author = {Cain, Mary B and {O' Neill}, Richard P and Castillo, Anya},
  pages = {1--36},
  title = {{History of optimal power flow and Models}},
  url = {https://www.ferc.gov/industries/electric/indus-act/market-planning/opf-papers/acopf-1-history-Model-testing.pdf}
  year = {2012}
}
```
"""
mutable struct ACRUPowerModel <: AbstractUnbalancedACRModel @pmd_fields end

""
abstract type AbstractUnbalancedIVRModel <: AbstractUnbalancedACRModel end

"""
Current voltage formulation of AC OPF. The formulation uses rectangular
coordinates for both current and voltage.  Note that, even though Kirchhoff's
circuit laws are linear in current and voltage, this formulation is nonconvex
due to constants power loads/generators and apparent power limits.
```
@techreport{ONeill2012,
    author = {{O' Neill}, Richard P and Castillo, Anya and Cain, Mary B},
    pages = {1--18},
    title = {{The IV formulation and linear approximations of the ac optimal power flow problem}},
    year = {2012}
}
```
Applicable to problem formulations with `_iv` in the name.
"""
mutable struct IVRUPowerModel <: AbstractUnbalancedIVRModel @pmd_fields end


##### Linear Approximations #####



abstract type AbstractUnbalancedDCPModel <: AbstractUnbalancedActivePowerModel end


"""
Linearized 'DC' power flow Model with polar voltage variables.
This model is a basic linear active-power-only approximation, which uses branch susceptance values
`br_b = -br_x / (br_x^2 + br_x^2)` for determining the network phase angles.  Furthermore, transformer
parameters such as tap ratios and phase shifts are not considered as part of this model.
It is important to note that it is also common for active-power-only approximations to use `1/br_x` for
determining the network phase angles, instead of the `br_b` value that is used here.  Small discrepancies
in solutions should be expected when comparing active-power-only approximations across multiple tools.
```
@ARTICLE{4956966,
  author={B. Stott and J. Jardim and O. Alsac},
  journal={IEEE Transactions on Power Systems},
  title={DC Power Flow Revisited},
  year={2009},
  month={Aug},
  volume={24},
  number={3},
  pages={1290-1300},
  doi={10.1109/TPWRS.2009.2021235},
  ISSN={0885-8950}
}
```
"""
mutable struct DCPUPowerModel <: AbstractUnbalancedDCPModel @pmd_fields end


abstract type AbstractUnbalancedNFAModel <: AbstractUnbalancedDCPModel end

"""
The an active power only network flow approximation, also known as the transportation model.
"""
mutable struct NFAUPowerModel <: AbstractUnbalancedNFAModel @pmd_fields end




"Base Abstract NLP Unbalanced Branch Flow Model"
abstract type AbstractNLPUBFModel <: AbstractUBFQPModel end


"Base Abstract Conic Unbalanced Branch Flow Model"
abstract type AbstractConicUBFModel <: AbstractUBFConicModel end


"Collection of Unbalanced Branch Flow Models"
AbstractUBFModels = Union{AbstractNLPUBFModel, AbstractConicUBFModel}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFModel <: AbstractConicUBFModel end


"SDP BFM with KCL as matrix equation, Geth 2020 (under review)"
abstract type SDPUBFKCLMXModel <: SDPUBFModel end


"Collection of Semidefinite Models"
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
mutable struct SDPUBFPowerModel <: SDPUBFModel @pmd_fields end


"default SDP unbalanced DistFlow with matrix KCL constructor"
mutable struct SDPUBFKCLMXPowerModel <: SDPUBFKCLMXModel @pmd_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCNLPUBFPowerModel <: SOCNLPUBFModel @pmd_fields end


"default SOC unbalanced DistFlow constructor"
mutable struct SOCConicUBFPowerModel <: SOCConicUBFModel @pmd_fields end


"default LP unbalanced DistFlow constructor"
mutable struct LPUBFDiagPowerModel <: LPUBFDiagModel @pmd_fields end
const LinDist3FlowPowerModel = LPUBFDiagPowerModel # more popular name


AbstractUnbalancedWModels = Union{AbstractUBFModel}
AbstractUnbalancedAPLossLessModels = Union{DCPUPowerModel, AbstractUnbalancedNFAModel}
AbstractUnbalancedPolarModels = Union{AbstractUnbalancedACPModel, AbstractUnbalancedDCPModel}
AbstractUnbalancedWConvexModels = Union{AbstractUBFModel}