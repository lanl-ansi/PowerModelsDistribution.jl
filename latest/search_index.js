var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ThreePhasePowerModels.jl-Documentation-1",
    "page": "Home",
    "title": "ThreePhasePowerModels.jl Documentation",
    "category": "section",
    "text": "CurrentModule = ThreePhasePowerModels"
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "ThreePowerModels.jl is a Julia/JuMP extension package to PowerModels.jl for modeling of Multi-Phase (with a focus on three-phase) power grids. "
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The latest stable release of PowerModels can be installed using the Julia package manager withPkg.add(\"ThreePhasePowerModels\")For the current development version, \"checkout\" this package withPkg.checkout(\"ThreePhasePowerModels\")At least one solver is required for running ThreePhasePowerModels.  The open-source solver Ipopt is recommended, as it is extremely fast, and can be used to solve a wide variety of the problems and network formulations provided in ThreePhasePowerModels.  The Ipopt solver can be installed via the package manager withPkg.add(\"Ipopt\")Test that the package works by runningPkg.test(\"ThreePhasePowerModels\")"
},

{
    "location": "quickguide.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "quickguide.html#Quick-Start-Guide-1",
    "page": "Getting Started",
    "title": "Quick Start Guide",
    "category": "section",
    "text": "Once ThreePhasePowerModels is installed, Ipopt is installed, and a network data file (e.g. \"case5_c_r_a.m\" or \"case3_unbalanced.dss\" in the package folder under ./test/data) has been acquired, an unbalanced AC Optimal Power Flow can be executed with,using ThreePhasePowerModels\nusing Ipopt\n\nrun_ac_tp_opf(\"case3_unbalanced.dss\", IpoptSolver())"
},

{
    "location": "quickguide.html#Getting-Results-1",
    "page": "Getting Started",
    "title": "Getting Results",
    "category": "section",
    "text": "The run commands in ThreePhasePowerModels return detailed results data in the form of a dictionary. Results dictionaries from either Matpower-style .m or OpenDSS\' .dss files will be identical in format. This dictionary can be saved for further processing as follows,result = run_ac_tp_opf(\"case3_unbalanced.dss\", IpoptSolver())"
},

{
    "location": "quickguide.html#Accessing-Different-Formulations-1",
    "page": "Getting Started",
    "title": "Accessing Different Formulations",
    "category": "section",
    "text": "The function \"run_ac_tp_opf\" is a shorthands for a more general formulation-independent OPF execution, \"run_tp_opf\". For example, run_ac_tp_opf is equivalent to,using PowerModels\nrun_tp_opf(\"case3_unbalanced.dss\", ACPPowerModel, IpoptSolver())Note that PowerModels needs to be loaded to access formulations which are extended by ThreePhasePowerModels, here \"ACPPowerModel\". The PowerModel \"ACPPowerModel\" indicates an AC formulation in polar coordinates.  This more generic run_tp_opf() allows one to solve an OPF problem with any power network formulation implemented in PowerModels or ThreePhasePowerModels.  For example, the SDP relaxation of unbalanced Optimal Power Flow (branch flow model) can be run with,using SCS\nrun_tp_opf_bf(\"case3_unbalanced.dss\", SDPUBFPowerModel, SCSSolver())Note that you have to use a SDP-capable solver, e.g. the open-source solver SCS, to solve SDP models."
},

{
    "location": "quickguide.html#Inspecting-the-Formulation-1",
    "page": "Getting Started",
    "title": "Inspecting the Formulation",
    "category": "section",
    "text": "The following example demonstrates how to break a run_tp_opf call into seperate model building and solving steps.  This allows inspection of the JuMP model created by ThreePhasePowerModels for the AC-OPF problem,data = ThreePhasePowerModels.parse_file(\"case3_unbalanced.dss\")\npm = PowerModels.build_generic_model(data, ACPPowerModel, ThreePhasePowerModels.post_tp_opf; multiconductor=true)\nprint(pm.model)\nsolve_generic_model(pm, IpoptSolver())"
},

{
    "location": "math-model.html#",
    "page": "Mathematical Model",
    "title": "Mathematical Model",
    "category": "page",
    "text": ""
},

{
    "location": "math-model.html#The-ThreePhasePowerModels-Mathematical-Model-1",
    "page": "Mathematical Model",
    "title": "The ThreePhasePowerModels Mathematical Model",
    "category": "section",
    "text": "As ThreePhasePowerModels implements a variety of power network optimization problems, the implementation is the best reference for precise mathematical formulations.  This section provides a complex number based mathematical specification for a prototypical unbalanced AC Optimal Power Flow problem, to provide an overview of the typical mathematical models in ThreePhasePowerModels."
},

{
    "location": "math-model.html#Unbalanced-AC-Optimal-Power-Flow-1",
    "page": "Mathematical Model",
    "title": "Unbalanced AC Optimal Power Flow",
    "category": "section",
    "text": "ThreePhasePowerModels implements a  generalized version of the AC Optimal Power Flow problem, from Matpower but extended to take into account phase unbalance [1].  These generalizations make it possible for ThreePhasePowerModels to more accurately capture real-world distribution network datasets.  The core generalizations are,Support for multiple load and shunt components on each bus\nLine charging (shunt) that supports a conductance and asymmetrical valuesIn the mathematical description below,Bold typeface indicates a vector (in mathbbC^c) or matrix (in mathbbC^ctimes c)\nOperator diag takes the diagonal (vector) from a square matrix\nThe set of complex numbers is mathbbC and real numbers is mathbbR\nSuperscript H indicates complex conjugate transpose (Hermitian adjoint)\nNote that complex power is defined as mathbfS_ij = mathbfV_i mathbfI_ij^H and is therefore a complex matrix of dimension c times c\nThe line mathbfY^c_ij mathbfY^c_ji and bus mathbfY^s_k shunt matrices do not need to be diagonal"
},

{
    "location": "math-model.html#Sets-1",
    "page": "Mathematical Model",
    "title": "Sets",
    "category": "section",
    "text": "The definitions of the sets involved remain unchanged w.r.t. the balanced OPF problem definition, except for the addition of the conductor set:beginalign\n\nmboxsets  nonumber \n N mbox - busesnonumber \n R mbox - references busesnonumber \n E E^R mbox - branches forward and reverse orientation nonumber \n G G_i mbox - generators and generators at bus i nonumber \n L L_i mbox - loads and loads at bus i nonumber \n S S_i mbox - shunts and shunts at bus i nonumber \n C mbox - conductors nonumber \n\nendalignwhere the set of conductors C typically equals  abc."
},

{
    "location": "math-model.html#Data-1",
    "page": "Mathematical Model",
    "title": "Data",
    "category": "section",
    "text": "beginalign\nmboxdata  nonumber \n S^gl_kc S^gu_kc in mathbbC  forall k in G forall c in C nonumber mathbfS^gl_k= S^gl_kc_c in C mathbfS^gu_k = S^gu_kc_c in C  \n c_2k c_1k c_0k in mathbbR  forall k in G nonumber \n v^l_ic v^u_ic in mathbbR  forall i in N forall c in C nonumber mathbfv^l_i = v^l_ic_c in C mathbfv^u_i = v^u_ic_c in C \n S^d_kcin mathbbC  forall k in L forall c in C nonumber mathbfS^d_k = S^d_kc_c in C \n mathbfY^s_kin mathbbC^ctimes c  forall k in S nonumber \n mathbfY_ij mathbfY^c_ij mathbfY^c_jiin mathbbC^ctimes c  forall (ij) in E nonumber \n s^u_ijc theta^Delta l_ijc theta^Delta u_ijc in mathbbR forall (ij) in E forall c in C nonumber mathbfs^u_ij = s^u_ijc_c in C \n V^textref_ic  in mathbbC  forall r in R  mathbfV^textref_i =  V^textref_ic_c in C \n\nendalignwhere the notation mathbfv^l_i = v^l_ic_c in C reflects that the vector mathbfv^l_i is constructed by putting the individual phase values v^l_ic in a vector (in order abc).Alternatively, the series impedance of a line can be written in impedance form:mathbfZ_ij in mathbbC^ctimes c  forall (ij) in E nonumber mathbfY_ij = ( mathbfZ_ij)^-1where superscript -1 indicates the matrix inverse. Note that mathbfY_ij or mathbfZ_ij may not be invertible, e.g. in case of single-phase branches in a three-phase grid. In this case the pseudo-inverse can be used."
},

{
    "location": "math-model.html#Variables-for-a-Bus-Injection-Model-1",
    "page": "Mathematical Model",
    "title": "Variables for a Bus Injection Model",
    "category": "section",
    "text": "beginalign\n S^g_kc  in mathbbC  forall kin G forall c in C nonumber mathbfS^g_k = S^g_kc_c in C \n V_ic  in mathbbC  forall iin N forall c in C nonumber mathbfV_i = V_ic_c in C \n mathbfS_ij  in mathbbC^ctimes c  forall (ij) in E cup E^R \n\nendalign"
},

{
    "location": "math-model.html#Mathematical-Formulation-of-a-Bus-Injection-Model-1",
    "page": "Mathematical Model",
    "title": "Mathematical Formulation of a Bus Injection Model",
    "category": "section",
    "text": "A complete mathematical model is as follows,\nbeginalign\nmboxminimize   sum_k in G c_2k left( sum_c in C Re(S^g_kc) right)^2 + c_1k  sum_c in C Re(S^g_kc) + c_0k \n\nmboxsubject to   nonumber \n mathbfV_i = mathbfV^textref_i    forall r in R \n S^gl_kc leq S^g_kc leq S^gu_kc  forall k in G forall c in C  \n v^l_ic leq V_ic leq v^u_ic  forall i in N forall c in C \n sum_substackk in G_i mathbfS^g_k - sum_substackk in L_i mathbfS^d_k - sum_substackk in S_i  mathbfV_i mathbfV^H_i (mathbfY^s_k)^H = sum_substack(ij)in E_i cup E_i^R diag(mathbfS_ij)  forall iin N \n mathbfS_ij =  mathbfV_i mathbfV_i^H left( mathbfY_ij + mathbfY^c_ijright)^H - mathbfV_i mathbfV^H_j mathbfY^H_ij   forall (ij)in E \n mathbfS_ji = mathbfV_j mathbfV_j^H left( mathbfY_ij + mathbfY^c_ji right)^H - mathbfV^H_i mathbfV_j mathbfY^H_ij  forall (ij)in E \n diag(mathbfS_ij) leq mathbfs^u_ij  forall (ij) in E cup E^R \n theta^Delta l_ijc leq angle (V_ic V^*_jc) leq theta^Delta u_ijc  forall (ij) in E forall c in C\n\nendalign"
},

{
    "location": "math-model.html#Variables-for-a-Branch-Flow-Model-1",
    "page": "Mathematical Model",
    "title": "Variables for a Branch Flow Model",
    "category": "section",
    "text": "beginalign\n S^g_kc  in mathbbC forall kin G forall c in C nonumber mathbfS^g_k = S^g_kc_c in C \n V_ic in mathbbC  forall iin N forall c in C nonumber mathbfV_i = V_ic_c in C \n I^s_ijc  in mathbbC forall e in E forall c in C nonumber mathbfI^s_ij = I^s_ijc_c in C \n mathbfS_ij  in mathbbC^ctimes c  forall (ij) in E cup E^R \n\nendalign"
},

{
    "location": "math-model.html#Mathematical-Formulation-of-a-Branch-Flow-Model-1",
    "page": "Mathematical Model",
    "title": "Mathematical Formulation of a Branch Flow Model",
    "category": "section",
    "text": "A complete mathematical model is as follows,beginalign\nmboxminimize   sum_k in G c_2k left( sum_c in C Re(S^g_kc) right)^2 + c_1k  sum_c in C Re(S^g_kc) + c_0k \n\nmboxsubject to   nonumber \n mathbfV_i = mathbfV^textref_i    forall r in R \n S^gl_kc leq S^g_kc leq S^gu_kc  forall k in G forall c in C  \n v^l_ic leq V_ic leq v^u_ic  forall i in N forall c in C \n sum_substackk in G_i mathbfS^g_k - sum_substackk in L_i mathbfS^d_k - sum_substackk in S_i  mathbfV_i mathbfV^H_i (mathbfY^s_k)^H = sum_substack(ij)in E_i cup E_i^R diag(mathbfS_ij)  forall iin N \n mathbfS_ij + mathbfS_ji =  mathbfV_i mathbfV_i^H (mathbfY^c_ij)^H + mathbfZ_ij mathbfI^s_ij(mathbfI^s_ij)^H + mathbfV_j mathbfV_j^H (mathbfY^c_ji)^H   forall (ij)in E \n mathbfS^s_ij = mathbfS_ij + mathbfV_i mathbfV_i^H (mathbfY^c_ij)^H   forall (ij) in E cup E^R \n mathbfS^s_ij = mathbfV_i (mathbfI^s_ij)^H   forall (ij) in E cup E^R\n mathbfV_i = mathbfV_j - mathbfZ_ij mathbfI^s_ij forall (ij)in E \n diag(mathbfS_ij) leq mathbfs^u_ij  forall (ij) in E cup E^R \n theta^Delta l_ijc leq angle (V_ic V^*_jc) leq theta^Delta u_ijc  forall (ij) in E forall c in C\n\nendalign[1] Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399"
},

{
    "location": "formulations.html#",
    "page": "Network Formulations",
    "title": "Network Formulations",
    "category": "page",
    "text": ""
},

{
    "location": "formulations.html#Network-Formulations-1",
    "page": "Network Formulations",
    "title": "Network Formulations",
    "category": "section",
    "text": ""
},

{
    "location": "formulations.html#Type-Hierarchy-1",
    "page": "Network Formulations",
    "title": "Type Hierarchy",
    "category": "section",
    "text": "We begin with the top of the hierarchy, where we can distinguish between conic and non-conic power flow models.PowerModels.AbstractConicForms = Union{PowerModels.AbstractConicPowerFormulation, PowerModels.AbstractBFConicForm}\nPowerModels.AbstractConicPowerFormulation <: PowerModels.AbstractPowerFormulation\nPowerModels.AbstractBFForm <: PowerModels.AbstractPowerFormulation\nPowerModels.AbstractBFQPForm <: PowerModels.AbstractBFForm\nPowerModels.AbstractBFConicForm <: PowerModels.AbstractBFFormWe begin with the top of the hierarchy, where we can distinguish between AC and DC power flow models.PowerModels.AbstractACPForm <: PowerModels.AbstractPowerFormulation\nPowerModels.AbstractDCPForm <: PowerModels.AbstractPowerFormulation\nPowerModels.AbstractWRForm <: PowerModels.AbstractPowerFormulation\nThreePhasePowerModels.AbstractNLPUBFForm <: PowerModels.AbstractBFQPForm\nThreePhasePowerModels.AbstractConicUBFForm <: PowerModels.AbstractBFConicForm\nThreePhasePowerModels.AbstractLPUBFForm <: ThreePhasePowerModels.AbstractNLPUBFFormFrom there, different forms are possible:#Bus injection models:\nPowerModels.StandardACPForm <: PowerModels.AbstractACPForm\nPowerModels.StandardDCPForm <: PowerModels.AbstractDCPForm\nPowerModels.SOCWRForm <: PowerModels.AbstractWRForm\n\n#Branch flow models:\nThreePhasePowerModels.SDPUBFForm <: ThreePhasePowerModels.AbstractConicUBFForm\nThreePhasePowerModels.SOCNLPUBFForm <: ThreePhasePowerModels.AbstractNLPUBFForm\nThreePhasePowerModels.SOCConicUBFForm <: ThreePhasePowerModels.AbstractConicUBFForm\n\nThreePhasePowerModels.LPLinUBFForm <: PowerModels.AbstractBFForm\nThreePhasePowerModels.LPfullUBFForm <: ThreePhasePowerModels.AbstractLPUBFForm\nThreePhasePowerModels.LPdiagUBFForm <: ThreePhasePowerModels.AbstractLPUBFForm"
},

{
    "location": "formulations.html#Power-Models-1",
    "page": "Network Formulations",
    "title": "Power Models",
    "category": "section",
    "text": "Each of these forms can be used as the type parameter for a PowerModel:PowerModels.ACPPowerModel = GenericPowerModel{PowerModels.StandardACPForm}\nPowerModels.DCPPowerModel = GenericPowerModel{PowerModels.StandardDCPForm}\n\nPowerModels.SOCWRPowerModel = GenericPowerModel{PowerModels.SOCWRForm}\n\nThreePhasePowerModels.SDPUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SDPUBFForm}\nThreePhasePowerModels.SOCNLPUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SOCNLPUBFForm}\nThreePhasePowerModels.SOCConicUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.SOCConicUBFForm}\n\nThreePhasePowerModels.LPfullUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPfullUBFForm}\nThreePhasePowerModels.LPdiagUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPdiagUBFForm}\nThreePhasePowerModels.LPLinUBFPowerModel = GenericPowerModel{ThreePhasePowerModels.LPLinUBFForm}"
},

{
    "location": "formulations.html#Union-Types-1",
    "page": "Network Formulations",
    "title": "Union Types",
    "category": "section",
    "text": "To support both conic and quadratically-constrained formulation variants for the unbalanced branch flow model, the union type AbstractUBFForm is defined. These formulations extend AbstractBFForm and are therefore also AbstractWForms (as defined in PowerModels proper).AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}"
},

{
    "location": "formulations.html#Optimization-problem-classes-1",
    "page": "Network Formulations",
    "title": "Optimization problem classes",
    "category": "section",
    "text": "NLP (nonconvex): ACPPowerModel\nSDP: SDPUBFPowerModel\nSOC(-representable): SOCWRPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel\nLinear: LPfullUBFPowerModel, LPdiagUBFPowerModel, LPLinUBFPowerModel, DCPPowerModel"
},

{
    "location": "formulations.html#Matrix-equations-versus-scalar-equations-1",
    "page": "Network Formulations",
    "title": "Matrix equations versus scalar equations",
    "category": "section",
    "text": "JuMP supports vectorized syntax, but not for nonlinear constraints. Therefore, certain formulations must be implemented in a scalar fashion. Other formulations can be written as matrix (in)equalities. The current implementations are categorized as follows:Scalar: ACPPowerModel, DCPPowerModel, LPLinUBFPowerModel, SOCWRPowerModel\nMatrix: SDPUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFPowerModel, LPfullUBFPowerModel, LPdiagUBFPowerModel"
},

{
    "location": "library.html#",
    "page": "Modeling Components",
    "title": "Modeling Components",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPLinUBFForm",
    "category": "type",
    "text": "LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current in scalar form\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPLinUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPLinUBFPowerModel",
    "category": "method",
    "text": "default Lin3Distflow constructor for scalar form\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPdiagUBFForm",
    "category": "type",
    "text": "LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPdiagUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPdiagUBFPowerModel",
    "category": "method",
    "text": "default LP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPfullUBFForm",
    "category": "type",
    "text": "Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPfullUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.LPfullUBFPowerModel",
    "category": "method",
    "text": "default LP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SDPUBFForm",
    "category": "type",
    "text": "SDP BFM per Gan and Low 2014, PSCC\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SDPUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SDPUBFPowerModel",
    "category": "method",
    "text": "default SDP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCConicUBFForm",
    "category": "type",
    "text": "SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCConicUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCConicUBFPowerModel",
    "category": "method",
    "text": "default SOC unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCNLPUBFForm",
    "category": "type",
    "text": "SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an QCP\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "category": "method",
    "text": "default SOC unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_ac_tp_opf-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_ac_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_ac_tp_pf-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_ac_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_bf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_bf-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_pbs-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_pbs-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_ots-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_ots-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_bf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_bf-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_pbs-Tuple{Dict{String,Any},Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_pbs-Tuple{String,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractConicUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.AbstractConicUBFForm",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractLPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.AbstractLPUBFForm",
    "category": "type",
    "text": "Abstract form for linear unbalanced power flow models\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractNLPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.AbstractNLPUBFForm",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#PowerModels.constraint_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "PowerModels.constraint_branch_current",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_branch_current_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_branch_current_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_branch_flow_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_branch_flow_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_bus_slack_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_bus_slack_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_bus_voltage_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_bus_voltage_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_component!-Tuple{Dict,AbstractString,Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_component!",
    "category": "method",
    "text": "add_component!(dss_data, ctype_name, compDict)\n\nAdds a component of type ctype_name with properties given by compDict to the existing dss_data structure. If a component of the same type has already been added to dss_data, the new component is appeneded to the existing array of components of that type, otherwise a new array is created.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_original_variables-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_original_variables",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_property-Tuple{Dict,AbstractString,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.add_property",
    "category": "method",
    "text": "add_property(compDict, key, value)\n\nAdds a property to an existing component properties dictionary compDict given the key and value of the property. If a property of the same name already exists inside compDict, the original value is converted to an array, and the new value is appended to the end.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.adjust_sourcegen_bounds!-Tuple{Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.adjust_sourcegen_bounds!",
    "category": "method",
    "text": "adjust_sourcegen_bounds!(tppm_data)\n\nChanges the bounds for the sourcebus generator by checking the emergamps of all of the branches attached to the sourcebus and taking the sum of non-infinite values. Defaults to Inf if all emergamps connected to sourcebus are also Inf.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.assign_property!-Tuple{Dict,AbstractString,AbstractString,AbstractString,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.assign_property!",
    "category": "method",
    "text": "assign_property!(dss_data, cType, cName, propName, propValue)\n\nAssigns a property with name propName and value propValue to the component of type cType named cName in dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.calc_tp_voltage_product_bounds-Tuple{PowerModels.GenericPowerModel,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.calc_tp_voltage_product_bounds",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.check_duplicate_components!-Tuple{Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.check_duplicate_components!",
    "category": "method",
    "text": "check_duplicate_components!(dss_data)\n\nFinds duplicate components in dss_data and merges up, meaning that older data (lower indices) is always overwritten by newer data (higher indices).\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.check_network_data-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.check_network_data",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractBFForm, PowerModels.AbstractWRConicForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Int64,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))\nq[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] == -b*(t[f_bus] - t[t_bus])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractWRConicForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))\nq[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])\nq[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))\nq[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Do nothing, this model is symmetric\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractWRConicForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))\nq[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "Do nothing, this model is symmetric\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])\nq[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.AbstractLPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.SDPUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.SOCConicUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.SOCConicUBFForm, ThreePhasePowerModels.SOCNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow loss equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACPForm, PowerModels.AbstractACTForm, PowerModels.AbstractDCPForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "Creates phase angle constraints at reference buses\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractBFForm, PowerModels.AbstractWRConicForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "do nothing, no way to represent this in these variables\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "do nothing, this model does not have complex voltage constraints\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_voltage_magnitude_difference-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_voltage_magnitude_difference",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.constraint_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createCapacitor",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createCapacitor",
    "category": "function",
    "text": "createCapacitor(bus1, name, bus2=0; kwargs)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Capacitor. If bus2 is not specified, the capacitor will be treated as a shunt. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createGenerator-Tuple{Any,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createGenerator",
    "category": "method",
    "text": "createGenerator(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Generator. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLine-Tuple{Any,Any,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createLine",
    "category": "method",
    "text": "createLine(bus1, bus2, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the properties for a Line. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLinecode-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createLinecode",
    "category": "method",
    "text": "createLinecode(name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the properties of a Linecode. See OpenDSS documentation for valid fields and ways to specify the different properties. DEPRECIATED: Calculation all done inside of createLine() due to Rg, Xg. Merge linecode values into line kwargs values BEFORE calling createLine(). This is now mainly used for parsing linecode dicts into correct data types.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLoad-Tuple{Any,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createLoad",
    "category": "method",
    "text": "createLoad(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Load. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createPVSystem-Tuple{Any,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createPVSystem",
    "category": "method",
    "text": "createPVSystem(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a PVSystem. See OpenDSS document https://github.com/tshort/OpenDSS/blob/master/Doc/OpenDSS%20PVSystem%20Model.doc for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createReactor",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createReactor",
    "category": "function",
    "text": "createReactor(bus1, name, bus2=0; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Reactor. If bus2 is not specified Reactor is treated like a shunt. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createTransformer-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createTransformer",
    "category": "method",
    "text": "createTransformer(name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Transformer. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createVSource",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.createVSource",
    "category": "function",
    "text": "createVSource(bus1, name, bus2=0; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Voltage Source. If bus2 is not specified, VSource will be treated like a generator. Mostly used as sourcebus which represents the circuit. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.create_starbus-Tuple{Dict,Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.create_starbus",
    "category": "method",
    "text": "creates a starbus from a 3-winding transformer\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.discover_buses-Tuple{Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.discover_buses",
    "category": "method",
    "text": "discover_buses(dss_data)\n\nDiscovers all of the buses (not separately defined in OpenDSS), from \"lines\".\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_branch!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_branch!",
    "category": "method",
    "text": "dss2tppm_branch!(tppm_data, dss_data)\n\nAdds PowerModels-style branches to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_bus!",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_bus!",
    "category": "function",
    "text": "dss2tppm_bus!(tppm_data, dss_data)\n\nAdds PowerModels-style buses to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_gen!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_gen!",
    "category": "method",
    "text": "dss2tppm_gen!(tppm_data, dss_data)\n\nAdds PowerModels-style generators to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_load!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_load!",
    "category": "method",
    "text": "dss2tppm_load!(tppm_data, dss_data)\n\nAdds PowerModels-style loads to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_pvsystem!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_pvsystem!",
    "category": "method",
    "text": "dss2tppm_pvsystem!(tppm_data, dss_data)\n\nAdds PowerModels-style pvsystems to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_shunt!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_shunt!",
    "category": "method",
    "text": "dss2tppm_shunt!(tppm_data, dss_data)\n\nAdds PowerModels-style shunts to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_transformer!-Tuple{Dict,Dict,Bool}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.dss2tppm_transformer!",
    "category": "method",
    "text": "dss2tppm_transformer!(tppm_data, dss_data, import_all)\n\nAdds PowerModels-style transformers (branches) to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.find_bus-Tuple{AbstractString,Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.find_bus",
    "category": "method",
    "text": "find_bus(busname, tppm_data)\n\nFinds the index number of the bus in existing data from the given busname.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.find_component-Tuple{Dict,AbstractString,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.find_component",
    "category": "method",
    "text": "find_component(tppm_data, name, compType)\n\nReturns the component of compType with name from data of type Dict{String,Array}.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_dtypes-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_dtypes",
    "category": "method",
    "text": "Returns a Dict{String,Type} for the desired component comp, giving all of the expected data types\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_linecode-Tuple{Dict,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_linecode",
    "category": "method",
    "text": "returns the linecode with name id\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_pbs_solution-Tuple{PowerModels.GenericPowerModel,Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_pbs_solution",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_prop_name-Tuple{AbstractString,Int64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_prop_name",
    "category": "method",
    "text": "get_prop_name(ctype, i)\n\nReturns the ith property name for a given component type ctype.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_prop_name-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_prop_name",
    "category": "method",
    "text": "get_prop_name(ctype)\n\nReturns the property names in order for a given component type ctype.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_solution_tp-Tuple{PowerModels.GenericPowerModel,Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.get_solution_tp",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_array-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.isa_array",
    "category": "method",
    "text": "checks if data is an opendss-style array string\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_conn-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.isa_conn",
    "category": "method",
    "text": "checks is a string is a connection by checking the values\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_matrix-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.isa_matrix",
    "category": "method",
    "text": "checks if data is an opendss-style matrix string\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_rpn-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.isa_rpn",
    "category": "method",
    "text": "detects if expr is Reverse Polish Notation expression\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.make_mpv!-Tuple{Dict{String,Any},String,Array{String,1}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.make_mpv!",
    "category": "method",
    "text": "collects several from_keys in an array and sets it to the to_key, removes from_keys\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.matlab_to_tppm-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.matlab_to_tppm",
    "category": "method",
    "text": "Converts a Matlab dict into a ThreePhasePowerModels dict\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.merge_dss!-Tuple{Dict{String,Array},Dict{String,Array}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.merge_dss!",
    "category": "method",
    "text": "merge_dss!(dss_prime, dss_to_add)\n\nMerges two (partially) parsed OpenDSS files to the same dictionary dss_prime. Used in cases where files are referenced via the \"compile\" or \"redirect\" OpenDSS commands inside the originating file.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_branch-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.ml2pm_branch",
    "category": "method",
    "text": "convert raw branch data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_bus-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.ml2pm_bus",
    "category": "method",
    "text": "convert raw bus data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_gen-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.ml2pm_gen",
    "category": "method",
    "text": "convert raw generator data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_load-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.ml2pm_load",
    "category": "method",
    "text": "convert raw load data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_shunt-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.ml2pm_shunt",
    "category": "method",
    "text": "convert raw shunt data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.objective_min_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.objective_min_bus_power_slack",
    "category": "method",
    "text": "a quadratic penalty for bus power slack variables\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_array",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_array",
    "category": "function",
    "text": "parse matrices according to active nodes\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_array-Tuple{Type,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_array",
    "category": "method",
    "text": "parse_array(dtype, data)\n\nParses a OpenDSS style array string data into a one dimensional array of type dtype. Array strings are capped by either brackets, single quotes, or double quotes, and elements are separated by spaces.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_buscoords-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_buscoords",
    "category": "method",
    "text": "parse_buscoords(file)\n\nParses a Bus Coordinate file, in either \"dat\" or \"csv\" formats, where in \"dat\", columns are separated by spaces, and in \"csv\" by commas. File expected to contain \"bus,x,y\" on each line.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_busname-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_busname",
    "category": "method",
    "text": "parse_busname(busname)\n\nParses busnames as defined in OpenDSS, e.g. \"primary.1.2.3.0\".\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_component",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_component",
    "category": "function",
    "text": "parse_component(component, properies, compDict=Dict{String,Any}())\n\nParses a component with properties into a compDict. If compDict is not defined, an empty dictionary will be used. Assumes that unnamed properties are given in order, but named properties can be given anywhere.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_conn-Tuple{String}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_conn",
    "category": "method",
    "text": "parses connection \"conn\" specification reducing to wye or delta\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_dss",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss-Tuple{IOStream}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_dss",
    "category": "method",
    "text": "parse_dss(filename)\n\nParses a OpenDSS file given by filename into a Dict{Array{Dict}}. Only supports components and options, but not commands, e.g. \"plot\" or \"solve\". Will also parse files defined inside of the originating DSS file via the \"compile\", \"redirect\" or \"buscoords\" commands.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss_with_dtypes!",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_dss_with_dtypes!",
    "category": "function",
    "text": "parse_dss_with_dtypes!(dss_data, toParse)\n\nParses the data in keys defined by toParse in dss_data using types given by the default properties from the get_prop_default function.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_element_with_dtype-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_element_with_dtype",
    "category": "method",
    "text": "parses the raw dss values into their expected data types\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_file-Tuple{IOStream}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_file",
    "category": "method",
    "text": "parse_file(io)\n\nParses the IOStream of a file into a Three-Phase PowerModels data structure.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_file-Tuple{String}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_file",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_line",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_line",
    "category": "function",
    "text": "parse_line(elements, curCompDict=Dict{String,Any}())\n\nParses an already separated line given by elements (an array) of an OpenDSS file into curCompDict. If not defined, curCompDict is an empty dictionary.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab-Tuple{IOStream}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_matlab",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab-Tuple{String}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_matlab",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab_string-Tuple{String}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_matlab_string",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matrix-Tuple{Type,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_matrix",
    "category": "method",
    "text": "parse_matrix(dtype, data)\n\nParses a OpenDSS style triangular matrix string data into a two dimensional array of type dtype. Matrix strings are capped by either parenthesis or brackets, rows are separated by \"|\", and columns are separated by spaces.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matrix-Union{Tuple{Array{T,N} where N,Array{Bool,N} where N,Int64,Any}, Tuple{Array{T,N} where N,Array{Bool,N} where N,Int64}, Tuple{Array{T,N} where N,Array{Bool,N} where N}, Tuple{T}} where T",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_matrix",
    "category": "method",
    "text": "parse matrices according to active nodes\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_opendss-Tuple{Dict}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_opendss",
    "category": "method",
    "text": "Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_opendss-Tuple{IOStream}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_opendss",
    "category": "method",
    "text": "Parses a DSS file into a PowerModels usable format.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_options-Tuple{Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_options",
    "category": "method",
    "text": "parse_options(options)\n\nParses options defined with the set command in OpenDSS.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_properties-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_properties",
    "category": "method",
    "text": "parse_properties(properties)\n\nParses a string of properties of a component type, character by character into an array with each element containing (if present) the property name, \"=\", and the property value.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_rpn",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.parse_rpn",
    "category": "function",
    "text": "parses Reverse Polish Notation expr\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf_bf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf_pbs-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_ots-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf_bf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf_pbs-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.post_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_psd_real-Tuple{Any,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_psd_real",
    "category": "method",
    "text": "real-valued SDP to SDP relaxation based on PSDness of principal minors, default is 3x3 SDP relaxation\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc-Tuple{Any,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc",
    "category": "method",
    "text": "See section 4.3 for complex to real PSD constraint transformation: @article{Fazel2001, author = {Fazel, M. and Hindi, H. and Boyd, S.P.}, title = {{A rank minimization heuristic with application to minimum order system approximation}}, doi = {10.1109/ACC.2001.945730}, journal = {Proc. American Control Conf.}, number = {2}, pages = {4734–4739}, url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730}, volume = {6}, year = {2001} }\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc_complex-Tuple{Any,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc_complex",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:\n\n@article{Kim2003,\nauthor = {Kim, S and Kojima, M and Yamashita, M},\ntitle = {{Second order cone programming relaxation of a positive semidefinite constraint}},\ndoi = {10.1080/1055678031000148696},\njournal = {Optimization Methods and Software},\nnumber = {5},\npages = {535--541},\nvolume = {18},\nyear = {2003}\n}\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc_complex_conic-Tuple{Any,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc_complex_conic",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2, applied to complex-value matrix,  as described in:\n\n@article{Kim2003,\nauthor = {Kim, S and Kojima, M and Yamashita, M},\ntitle = {{Second order cone programming relaxation of a positive semidefinite constraint}},\ndoi = {10.1080/1055678031000148696},\njournal = {Optimization Methods and Software},\nnumber = {5},\npages = {535--541},\nvolume = {18},\nyear = {2003}\n}\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc_conic-Tuple{Any,Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc_conic",
    "category": "method",
    "text": "See section 4.3 for complex to real PSD constraint transformation: @article{Fazel2001, author = {Fazel, M. and Hindi, H. and Boyd, S.P.}, title = {{A rank minimization heuristic with application to minimum order system approximation}}, doi = {10.1109/ACC.2001.945730}, journal = {Proc. American Control Conf.}, number = {2}, pages = {4734–4739}, url = {http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=945730}, volume = {6}, year = {2001} }\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc_real-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc_real",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:\n\n@article{Kim2003,\nauthor = {Kim, S and Kojima, M and Yamashita, M},\ntitle = {{Second order cone programming relaxation of a positive semidefinite constraint}},\ndoi = {10.1080/1055678031000148696},\njournal = {Optimization Methods and Software},\nnumber = {5},\npages = {535--541},\nvolume = {18},\nyear = {2003}\n}\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.relaxation_psd_to_soc_real_conic-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.relaxation_psd_to_soc_real_conic",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2, applied to real-value matrix,  as described in:\n\n@article{Kim2003,\nauthor = {Kim, S and Kojima, M and Yamashita, M},\ntitle = {{Second order cone programming relaxation of a positive semidefinite constraint}},\ndoi = {10.1080/1055678031000148696},\njournal = {Optimization Methods and Software},\nnumber = {5},\npages = {535--541},\nvolume = {18},\nyear = {2003}\n}\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.roll-Union{Tuple{Array{T,1},Int64}, Tuple{T}} where T<:Number",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.roll",
    "category": "method",
    "text": "rolls a 1d array left or right by idx\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_dc_tp_pf-Tuple{Any,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.run_dc_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.set_default-Tuple{Dict{String,Any},String,Any}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.set_default",
    "category": "method",
    "text": "checks if the given dict has a value, if not, sets a default value\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.sqr-Tuple{Float64}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.sqr",
    "category": "method",
    "text": "Squares x, for parsing Reverse Polish Notation\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.strip_comments-Tuple{AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.strip_comments",
    "category": "method",
    "text": "Strips comments, defined by \"!\" from the ends of lines\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.strip_lines-Tuple{Array}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.strip_lines",
    "category": "method",
    "text": "strips lines that are either commented (block or single) or empty\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.to_sym_keys-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.to_sym_keys",
    "category": "method",
    "text": "converts Dict{String,Any} to Dict{Symbol,Any} for passing as kwargs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.translate_version!-Tuple{Dict{String,Any}}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.translate_version!",
    "category": "method",
    "text": "Translates legacy versions into current version format\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_active_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_active_bus_power_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_bus_power_slack",
    "category": "method",
    "text": "generates variables for both active and reactive slack at each bus\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_lower_triangle_active_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_lower_triangle_active_branch_flow",
    "category": "method",
    "text": "variable: p_lt[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_lower_triangle_reactive_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_lower_triangle_reactive_branch_flow",
    "category": "method",
    "text": "variable: q_lt[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_reactive_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_reactive_bus_power_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_magnitude_sqr-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage_magnitude_sqr",
    "category": "method",
    "text": "variable: w[i] >= 0 for i in buses\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_prod_hermitian-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage_prod_hermitian",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_prod_hermitian-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage_prod_hermitian",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_product-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_tp_voltage_product",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_upper_triangle_active_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_upper_triangle_active_branch_flow",
    "category": "method",
    "text": "variable: p_ut[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_upper_triangle_reactive_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.variable_upper_triangle_reactive_branch_flow",
    "category": "method",
    "text": "variable: q_ut[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.where_is_comp-Tuple{Array,AbstractString}",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.where_is_comp",
    "category": "method",
    "text": "where_is_comp(data, comp_id)\n\nFinds existing component of id comp_id in array of data and returns index. Assumes all components in data are unique.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.jl-Library-1",
    "page": "Modeling Components",
    "title": "ThreePhasePowerModels.jl Library",
    "category": "section",
    "text": "Modules = [ThreePhasePowerModels]"
},

{
    "location": "developer.html#",
    "page": "Developer",
    "title": "Developer",
    "category": "page",
    "text": ""
},

{
    "location": "developer.html#Developer-Documentation-1",
    "page": "Developer",
    "title": "Developer Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "formulation-details.html#",
    "page": "Formulation Details",
    "title": "Formulation Details",
    "category": "page",
    "text": ""
},

{
    "location": "formulation-details.html#Three-phase-formulation-details-1",
    "page": "Formulation Details",
    "title": "Three-phase formulation details",
    "category": "section",
    "text": ""
},

{
    "location": "formulation-details.html#StandardACPForm-1",
    "page": "Formulation Details",
    "title": "StandardACPForm",
    "category": "section",
    "text": "Real-valued formulation from:Formulation without shunts: Mahdad, B., Bouktir, T., & Srairi, K. (2006). A three-phase power flow modelization: a tool for optimal location and control of FACTS devices in unbalanced power systems. In IEEE Industrial Electronics IECON (pp. 2238–2243)."
},

{
    "location": "formulation-details.html#StandardDCPForm-1",
    "page": "Formulation Details",
    "title": "StandardDCPForm",
    "category": "section",
    "text": "Applying all of the standard DC linearization tricks to the StandardACPForm"
},

{
    "location": "formulation-details.html#SOCWRForm-1",
    "page": "Formulation Details",
    "title": "SOCWRForm",
    "category": "section",
    "text": "Applying the standard BIM voltage cross-product (sine and cosine) substitution tricks to StandardACPForm results immediately in a SOC formulation."
},

{
    "location": "formulation-details.html#SDPUBFForm-1",
    "page": "Formulation Details",
    "title": "SDPUBFForm",
    "category": "section",
    "text": "The BFM SDP relaxation as described in:Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399Note that this formulation is complex-valued and additional steps are needed to implement this in JuMP."
},

{
    "location": "formulation-details.html#SOCNLPUBFForm-1",
    "page": "Formulation Details",
    "title": "SOCNLPUBFForm",
    "category": "section",
    "text": "The starting point is SDPUBFForm. The SDP constraint can be relaxed to a set of SOC constraints, starting from either the real or complex form of the matrix on which the PSD-ness constraint is applied.Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696\nAndersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863."
},

{
    "location": "formulation-details.html#SOCConicUBFForm-1",
    "page": "Formulation Details",
    "title": "SOCConicUBFForm",
    "category": "section",
    "text": "See SOCNLPUBFForm"
},

{
    "location": "formulation-details.html#LPfullUBFForm-1",
    "page": "Formulation Details",
    "title": "LPfullUBFForm",
    "category": "section",
    "text": "Matrix formulation that generalizes simplified DistFlow equations, as introduced in :Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399Note that this formulation is complex-valued and additional steps are needed to implement this in JuMP."
},

{
    "location": "formulation-details.html#LPdiagUBFForm-1",
    "page": "Formulation Details",
    "title": "LPdiagUBFForm",
    "category": "section",
    "text": "This formulation has originally been developed by Sankur et al.Sankur, M. D., Dobbe, R., Stewart, E., Callaway, D. S., & Arnold, D. B. (2016). A linearized power flow model for optimization in unbalanced distribution systems. https://arxiv.org/abs/1606.04492v2This formulation is here cast as only considering the diagonal elements defined in LPfullUBFForm, which furthermore leads to the imaginary part of the lifted node voltage variable W being redundant and substituted out."
},

{
    "location": "formulation-details.html#LPLinUBFForm-1",
    "page": "Formulation Details",
    "title": "LPLinUBFForm",
    "category": "section",
    "text": "Scalar reformulation of:Sankur, M. D., Dobbe, R., Stewart, E., Callaway, D. S., & Arnold, D. B. (2016). A linearized power flow model for optimization in unbalanced distribution systems. https://arxiv.org/abs/1606.04492v2This formulation was already derived in real variables and parameters."
},

]}
