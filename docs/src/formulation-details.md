# Three-phase formulation details


## `StandardACPForm`
Real-valued formulation from:
- Formulation without shunts: Mahdad, B., Bouktir, T., & Srairi, K. (2006). A three-phase power flow modelization: a tool for optimal location and control of FACTS devices in unbalanced power systems. In IEEE Industrial Electronics IECON (pp. 2238–2243).



## `StandardDCPForm`
Applying all of the standard DC linearization tricks to the `StandardACPForm`

## `SOCWRForm`
Applying the standard BIM voltage cross-product (sine and cosine) substitution tricks to `StandardACPForm` results immediately in a SOC formulation.

## `SDPUBFForm`
The BFM SDP relaxation as described in:
- Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399
Note that this formulation is complex-valued and additional steps are needed to implement this in JuMP.

## `SOCNLPUBFForm`
The starting point is `SDPUBFForm`. The SDP constraint can be relaxed to a set of SOC constraints, starting from either the real or complex form of the matrix on which the PSD-ness constraint is applied.
- Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696
- Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.


## `SOCConicUBFForm`
See `SOCNLPUBFForm`


## `LPfullUBFForm`
Matrix formulation that generalizes `simplified DistFlow equations`, as introduced in :
- Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399

Note that this formulation is complex-valued and additional steps are needed to implement this in JuMP.

## `LPdiagUBFForm`
This formulation has originally been developed by Sankur et al.
- Sankur, M. D., Dobbe, R., Stewart, E., Callaway, D. S., & Arnold, D. B. (2016). A linearized power flow model for optimization in unbalanced distribution systems. https://arxiv.org/abs/1606.04492v2

This formulation is here cast as only considering the diagonal elements defined in `LPfullUBFForm`, which furthermore leads to the imaginary part of the lifted node voltage variable W being redundant and substituted out.

## `LPLinUBFForm`
Scalar reformulation of:
- Sankur, M. D., Dobbe, R., Stewart, E., Callaway, D. S., & Arnold, D. B. (2016). A linearized power flow model for optimization in unbalanced distribution systems. https://arxiv.org/abs/1606.04492v2
This formulation was already derived in real variables and parameters.
