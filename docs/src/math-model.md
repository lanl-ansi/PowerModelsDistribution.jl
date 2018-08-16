# The ThreePhasePowerModels Mathematical Model

As ThreePhasePowerModels implements a variety of power network optimization problems, the implementation is the best reference for precise mathematical formulations.  This section provides a complex number based mathematical specification for a prototypical unbalanced AC Optimal Power Flow problem, to provide an overview of the typical mathematical models in ThreePhasePowerModels.


## Preliminaries

ThreePhasePowerModels implements a  generalized version of the AC Optimal Power Flow problem. These generalizations make it possible for ThreePhasePowerModels to more accurately capture real-world distribution network datasets.  
- The core representation is based on [Matpower](http://www.pserc.cornell.edu/matpower/), but with the generalizations [made in PowerModels](https://lanl-ansi.github.io/PowerModels.jl/stable/math-model.html).
- This is further extended to take into account phase unbalance, largely in line with [^1].  

In the mathematical description below,
- The ordered set of conductors $C$ typically equals $\{ a,b,c\}$; the amount of conductors is denoted $|C|$.
- The set of complex numbers is $\mathbb{C}$ and real numbers is $\mathbb{R}$.
- A bold typeface, .e.g. "$\mathbf{S}$", indicates either a vector ($\in \mathbb{C}^{|C|}$) or matrix ($\in \mathbb{C}^{|C|\times |C|}$); for scalar variables a normal font type, e.g. "$S$" is used.
- Operator $diag$ takes the diagonal (vector) from a square matrix.
- Superscript $H$ indicates complex conjugate transpose ([Hermitian transpose](https://en.wikipedia.org/wiki/Hermitian_adjoint)).
- Note that complex power is defined as $\mathbf{S}_{ij} = \mathbf{V}_{i} \mathbf{I}_{ij}^H$ and is therefore a complex matrix of dimension $|C| \times |C|$.
- The branch series impedance $\mathbf{Z}_{ij}$ or admittance $\mathbf{Y}_{ij}$  is a complex matrix of dimension $|C| \times |C|$.
- Series impedance and admittance are interchangeable through the matrix inverse $\mathbf{Y}_{ij} = ( \mathbf{Z}_{ij})^{-1}$.  Note that $\mathbf{Y}_{ij}$ or $\mathbf{Z}_{ij}$ may not be invertible, e.g. in case of single-phase branches in a three-phase grid. Therefore, case the [pseudo-inverse](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse) is used.
- The branch $\mathbf{Y}^c_{ij}, \mathbf{Y}^c_{ji}$ and bus $\mathbf{Y}^s_{k}$ shunt matrices do not need to be diagonal.


## Sets

The definitions of the sets involved remain unchanged w.r.t. the balanced OPF problem definition, except for the addition of the conductor set:
```math
\begin{align}
%
\mbox{sets:} & \nonumber \\
& N \mbox{ - buses}\nonumber \\
& R \mbox{ - references buses}\nonumber \\
& E, E^R \mbox{ - branches, forward and reverse orientation} \nonumber \\
& G, G_i \mbox{ - generators and generators at bus $i$} \nonumber \\
& L, L_i \mbox{ - loads and loads at bus $i$} \nonumber \\
& S, S_i \mbox{ - shunts and shunts at bus $i$} \nonumber \\
& C \mbox{ - conductors (ordered)} \nonumber \\
%
\end{align}
```

## Parameters
The following data is defined in the same way as [PowerModels](https://lanl-ansi.github.io/PowerModels.jl/stable/math-model.html), except for the addition of a conductor index "$c$":
```math
\begin{align}
\mbox{scalar data:} & \nonumber \\
& S^{gl}_{k,c}, S^{gu}_{k,c} \in \mathbb{C} \;\; \forall k \in G, \forall c \in C \nonumber \\
& c_{2k}, c_{1k}, c_{0k} \in \mathbb{R} \;\; \forall k \in G \nonumber \\
& v^l_{i,c}, v^u_{i,c} \in \mathbb{R} \;\; \forall i \in N, \forall c \in C \nonumber\\
& S^d_{k,c}\in \mathbb{C} \;\; \forall k \in L, \forall c \in C \nonumber \\
& {s^u}_{ij,c}, \in \mathbb{R} \;\; \forall (i,j) \in E, \forall c \in C \nonumber\\
& {i^u}_{ij,c}, \in \mathbb{R} \;\; \forall (i,j) \in E, \forall c \in C \nonumber\\
& \theta^{\Delta l}_{ij,c}, \theta^{\Delta u}_{ij,c} \in \mathbb{R}\;\; \forall (i,j) \in E, \forall c \in C \nonumber\\
& V^{\text{ref}}_{i,c}  \in \mathbb{C} \;\; \forall r \in R, \forall c \in C \nonumber\\
%
\end{align}
```

The following vector symbols are defined analogously for notational simplicity:
```math
\begin{align}\mbox{vector data:} & \nonumber \\
& \mathbf{S}^{gl}_{k}:= [S^{gl}_{k,c}]_{c \in C}, \mathbf{S}^{gu}_{k} := [S^{gu}_{k,c}]_{c \in C}  \nonumber  \mbox{ - generator complex power bounds}\\
& \mathbf{v}^l_{i} := [v^l_{i,c}]_{c \in C}, \mathbf{v}^u_{i} := [v^u_{i,c}]_{c \in C}  \nonumber  \mbox{ - bus voltage bounds}\\
& \mathbf{S}^d_{k} := [S^d_{k,c}]_{c \in C}  \nonumber  \mbox{ - load complex power consumption}\\
& {\mathbf{s}^u}_{ij} := [{s^u}_{ij,c}]_{c \in C}  \nonumber  \mbox{ - branch apparent power flow limit}\\
& {\mathbf{i}^u}_{ij} := [{i^u}_{ij,c}]_{c \in C}  \nonumber \mbox{ - branch current limit}\\
& \mathbf{V}^{\text{ref}}_{i} =  [V^{\text{ref}}_{i,c}]_{c \in C}  \nonumber \mbox{ - reference bus voltage set point}\\
%
\end{align}
```
where the notation $\mathbf{v}^l_{i} := [v^l_{i,c}]_{c \in C}$ reflects that the vector $\mathbf{v}^l_{i}$ is constructed by putting the individual phase values $v^l_{i,c}$ in a vector (in order $a,b,c$).

Impedance and admittance parameters are matrices:
```math
\begin{align}\mbox{matrix data:} & \nonumber \\
& \mathbf{Y}^s_{k}\in \mathbb{C}^{c\times c} \;\; \forall k \in S \nonumber \mbox{ - bus shunt admittance}\\
& \mathbf{Y}_{ij}, \mathbf{Y}^c_{ij}, \mathbf{Y}^c_{ji}\in \mathbb{C}^{c\times c} \;\; \forall (i,j) \in E \nonumber \mbox{ - branch pi-section parameters}
%
\end{align}
```


##  Unbalanced AC Optimal Power Flow  (Bus Injection Model)


### Variables

```math
\begin{align}
& S^g_{k,c}  \in \mathbb{C} \;\; \forall k\in G, \forall c \in C \nonumber\\
& V_{i,c}  \in \mathbb{C} \;\; \forall i\in N, \forall c \in C \nonumber\\
& \mathbf{S}_{ij}  \in \mathbb{C}^{c\times c} \;\; \forall (i,j) \in E \cup E^R \nonumber\\
& \mathbf{S}^g_{k} := [S^g_{k,c}]_{c \in C} \nonumber\\
& \mathbf{V}_{i} := [V_{i,c}]_{c \in C} \nonumber\\
%
\end{align}
```


###  Mathematical Formulation

A complete mathematical model is as follows,

```math

\begin{align}
\mbox{minimize: } & \sum_{k \in G} c_{2k} \left( \sum_{c \in C} \Re(S^g_{k,c}) \right)^2 + c_{1k}  \sum_{c \in C} \Re(S^g_{k,c}) + c_{0k} \label{eq_objective}\\
%
\mbox{subject to: } & \nonumber \\
& \mathbf{V}_{i} = \mathbf{V}^{\text{ref}}_{i}   \;\; \forall r \in R \label{eq_ref_bus}\\
& S^{gl}_{k,c} \leq S^g_{k,c} \leq S^{gu}_{k,c} \;\; \forall k \in G, \forall c \in C  \label{eq_gen_bounds}\\
& v^l_{i,c} \leq |V_{i,c}| \leq v^u_{i,c} \;\; \forall i \in N, \forall c \in C \label{eq_voltage_bounds}\\
& \sum_{\substack{k \in G_i}} \mathbf{S}^g_k - \sum_{\substack{k \in L_i}} \mathbf{S}^d_k - \sum_{\substack{k \in S_i}}  \mathbf{V}_i \mathbf{V}^H_i (\mathbf{Y}^s_k)^H = \sum_{\substack{(i,j)\in E_i \cup E_i^R}} diag(\mathbf{S}_{ij}) \;\; \forall i\in N \label{eq_kcl_shunt}\\
& \mathbf{S}_{ij} =  {\mathbf{V}_i \mathbf{V}_i^H} \left( \mathbf{Y}_{ij} + \mathbf{Y}^c_{ij}\right)^H - {\mathbf{V}_i \mathbf{V}^H_j} \mathbf{Y}^H_{ij}  \;\; \forall (i,j)\in E \label{eq_power_from}\\
& \mathbf{S}_{ji} = \mathbf{V}_j \mathbf{V}_j^H \left( \mathbf{Y}_{ij} + \mathbf{Y}^c_{ji} \right)^H - { \mathbf{V}_j} \mathbf{V}^H_i \mathbf{Y}^H_{ij} \;\; \forall (i,j)\in E \label{eq_power_to}\\
& |diag(\mathbf{S}_{ij})| \leq \mathbf{s}^u_{ij} \;\; \forall (i,j) \in E \cup E^R \label{eq_thermal_limit}\\
& |\mathbf{I}_{ij}| \leq \mathbf{i}^u_{ij} \;\; \forall (i,j) \in E \cup E^R \label{eq_current_limit}\\
& \theta^{\Delta l}_{ij,c} \leq \angle (V_{i,c} V^*_{j,c}) \leq \theta^{\Delta u}_{ij,c} \;\; \forall (i,j) \in E, \forall c \in C \label{eq_angle_difference}
%
\end{align}
```


### Mapping to function names
- Eq. $\eqref{eq_objective}$ - `PowerModels.objective_min_fuel_cost`
- Eq. $\eqref{eq_ref_bus}$ - `constraint_tp_theta_ref` in `constraint_template.jl`
- Eq. $\eqref{eq_gen_bounds}$ - bounds of `PowerModels.variable_generation`
- Eq. $\eqref{eq_voltage_bounds}$ - bounds of `variable_tp_voltage`
- Eq. $\eqref{eq_kcl_shunt}$ - `PowerModels.constraint_kcl_shunt`
- Eq. $\eqref{eq_power_from}$ - `constraint_ohms_tp_yt_from` in `constraint_template.jl`
- Eq. $\eqref{eq_power_to}$ - `constraint_ohms_tp_yt_to` in `constraint_template.jl`
- Eq. $\eqref{eq_thermal_limit}$ - `PowerModels.constraint_thermal_limit_from` and `PowerModels.constraint_thermal_limit_to`
- Eq. $\eqref{eq_current_limit}$ - `PowerModels.constraint_current_limit`
- Eq. $\eqref{eq_angle_difference}$ - `PowerModels.constraint_voltage_angle_difference`


##  Unbalanced AC Optimal Power Flow (Branch Flow Model)

### Variables

```math
\begin{align}
& S^g_{k,c}  \in \mathbb{C}\;\; \forall k\in G, \forall c \in C \nonumber\\
& V_{i,c} \in \mathbb{C} \;\; \forall i\in N, \forall c \in C \nonumber\\
& I^{s}_{ij,c}  \in \mathbb{C}\;\; \forall e \in E, \forall c \in C \nonumber \\
& \mathbf{S}_{ij}  \in \mathbb{C}^{c\times c} \;\; \forall (i,j) \in E \cup E^R \nonumber\\
& \mathbf{S}^g_{k} := [S^g_{k,c}]_{c \in C} \nonumber\\
& \mathbf{V}_{i} := [V_{i,c}]_{c \in C} \nonumber\\
& \mathbf{I}^{s}_{ij} := [{I}^{s}_{ij,c}]_{c \in C} \nonumber\\
%
\end{align}
```


###  Mathematical Formulation

A complete mathematical model is as follows,

```math
\begin{align}
\mbox{minimize: } & \sum_{k \in G} c_{2k} \left( \sum_{c \in C} \Re(S^g_{k,c}) \right)^2 + c_{1k}  \sum_{c \in C} \Re(S^g_{k,c}) + c_{0k} \nonumber \\
%
\mbox{subject to: } & \nonumber \\
& \mathbf{V}_{i} = \mathbf{V}^{\text{ref}}_{i}   \;\; \forall r \in R \nonumber\\
& S^{gl}_{k,c} \leq S^g_{k,c} \leq S^{gu}_{k,c} \;\; \forall k \in G, \forall c \in C   \nonumber\\
& v^l_{i,c} \leq |V_{i,c}| \leq v^u_{i,c} \;\; \forall i \in N, \forall c \in C  \nonumber\\
& \sum_{\substack{k \in G_i}} \mathbf{S}^g_k - \sum_{\substack{k \in L_i}} \mathbf{S}^d_k - \sum_{\substack{k \in S_i}}  \mathbf{V}_i \mathbf{V}^H_i (\mathbf{Y}^s_k)^H = \sum_{\substack{(i,j)\in E_i \cup E_i^R}} diag(\mathbf{S}_{ij}) \;\; \forall i\in N  \nonumber\\
& \mathbf{S}_{ij} + \mathbf{S}_{ji} =  \mathbf{V}_i \mathbf{V}_i^H (\mathbf{Y}^c_{ij})^H + \mathbf{Z}_{ij} \mathbf{I}^{s}_{ij}(\mathbf{I}^{s}_{ij})^H + \mathbf{V}_j \mathbf{V}_j^H (\mathbf{Y}^c_{ji})^H  \;\; \forall (i,j)\in E \label{eq_line_losses}\\
& \mathbf{S}^{s}_{ij} = \mathbf{S}_{ij} + \mathbf{V}_i \mathbf{V}_i^H (\mathbf{Y}^c_{ij})^H  \;\; \forall (i,j) \in E \cup E^R  \label{eq_series_power_flow} \\
& \mathbf{S}^{s}_{ij} = \mathbf{V}_i (\mathbf{I}^{s}_{ij})^H  \;\; \forall (i,j) \in E \cup E^R  \label{eq_complex_power_definition}\\
& \mathbf{V}_i = \mathbf{V}_j - \mathbf{Z}_{ij} \mathbf{I}^{s}_{ij} \forall (i,j)\in E  \label{eq_ohms_bfm}\\
& |diag(\mathbf{S}_{ij})| \leq \mathbf{s}^u_{ij} \;\; \forall (i,j) \in E \cup E^R  \nonumber\\
& |\mathbf{I}_{ij}| \leq \mathbf{i}^u_{ij} \;\; \forall (i,j) \in E \cup E^R \nonumber \\
& \theta^{\Delta l}_{ij,c} \leq \angle (V_{i,c} V^*_{j,c}) \leq \theta^{\Delta u}_{ij,c} \;\; \forall (i,j) \in E, \forall c \in C  \nonumber
%
\end{align}
```

### Mapping to function names
- Eq. $\eqref{eq_line_losses}$ - `constraint_tp_flow_losses` in `constraint_template.jl`
- Eq. $\eqref{eq_series_power_flow}$ - implicit, substituted out in implementation
- Eq. $\eqref{eq_complex_power_definition}$ - `constraint_tp_branch_current` in `constraint_template.jl`
- Eq. $\eqref{eq_ohms_bfm}$ - `constraint_tp_voltage_magnitude_difference` in `constraint_template.jl`


[^1] Gan, L., & Low, S. H. (2014). Convex relaxations and linear approximation for optimal power flow in multiphase radial networks. In PSSC (pp. 1–9). Wroclaw, Poland. https://doi.org/10.1109/PSCC.2014.7038399
