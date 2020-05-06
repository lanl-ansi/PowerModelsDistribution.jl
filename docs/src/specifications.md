# Problem Specifications

In addition to the standard power flow `run_mc_pf`, and optimal power flow `run_mc_opf`, there are several notable problem specifications included in PowerModelsDistribution

## Optimal Power Flow (OPF) with On-Load Tap Changers (OLTC)

This problem is identical to `mc_opf`, except that all transformers are now modelled as on-load tap changers (OLTCs). Each phase has an individual tap ratio, which can be either variable or fixed, as specified in the data model.

### OLTC Objective

```julia
objective_min_fuel_cost(pm)
```

### OLTC Variables

```julia
variable_mc_voltage(pm)
variable_mc_branch_flow(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, cnd=c)
    PMs.variable_dcline_flow(pm, cnd=c)
end
variable_mc_transformer_flow(pm)
variable_mc_oltc_tap(pm)
```

### OLTC Constraints

```julia
constraint_mc_model_voltage(pm)

for i in PMs.ids(pm, :ref_buses)
    constraint_mc_theta_ref(pm, i)
end

for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
    constraint_mc_power_balance(pm, i, cnd=c)
end

for i in PMs.ids(pm, :branch)
    constraint_mc_ohms_yt_from(pm, i)
    constraint_mc_ohms_yt_to(pm, i)

    for c in PMs.conductor_ids(pm)
        PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end
end

for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_dcline(pm, i, cnd=c)
end

for i in PMs.ids(pm, :transformer)
    constraint_mc_oltc(pm, i)
end
```

## Minimal Load Delta (MLD) Problem Specification

Load shed (continuous) problem. See "Relaxations of AC Maximal Load Delivery for Severe Contingency Analysis" by C. Coffrin _et al._ (DOI: [10.1109/TPWRS.2018.2876507](https://ieeexplore.ieee.org/document/8494809)) for single-phase case.

### MLD Variables

```math
\begin{align}
\mbox{variables: } & \nonumber \\
& z^v_i \in \{0,1\}\ \ \forall i \in N \mbox{ - bus voltage on/off variable} \\
& z^g_i \in \{0,1\}\ \ \forall i \in G \mbox{ - generator on/off variable} \\
& z^{b}_i \in \{0,1\}\ \ \forall i \in B\mbox{ - storage on/off variable} \\
& z^d_i \in (0,1)\ \ \forall i \in L \mbox{ - continuous load shedding variable} \\
& z^s_i \in (0,1)\ \ \forall i \in H \mbox{ - continuous shunt shedding variable}
\end{align}
```

### MLD Objective

```math
\begin{align}
\mbox{minimize: }\left (
\sum_{\substack{i\in N,c\in C}}{10 \left (1-z^v_i \right )} + \sum_{\substack{i\in L,c\in C}}{10 \omega_{i,c}\left |\Re{\left (S^d_i\right )}\right |\left ( 1-z^d_i \right ) } + \sum_{\substack{i\in H,c\in C}}{\left | \Re{\left (S^s_i \right )}\right | \left (1-z^s_i \right ) } + \sum_{\substack{i\in G,c\in C}}{\Delta^g_i } + \sum_{\substack{i\in B,c\in C}}{\Delta^b_i} \right )
\end{align}
```

where

```math
\begin{align}
\Delta^g_i &>= \left [\Re{\left (S^g_{i}(0) \right )} - \Re{\left (S^g_i \right )} \right ] \\
\Delta^g_i &>= -\left [\Re{\left (S^g_{i}(0) \right )} - \Re{\left (S^g_i \right )} \right ] \\
\Delta^b_i &>= \left [\Re{\left (S^b_{i}(0) \right )} - \Re{\left (S^b_i \right )} \right ] \\
\Delta^b_i &>= -\left [\Re{\left (S^b_{i}(0) \right )} - \Re{\left (S^b_i \right )} \right ]
\end{align}
```

### MLD Constraints

```math
\begin{align}
\mbox{subject to: } & \nonumber \\
& z^v_i v^l_{i,c} \leq \left | V_{i,c} \right | \leq z_i^v v^u_{i,c}\ \ \forall i \in N,\forall c \in C \\
& z^g_i S^{gl}_{i,c} \leq S^g_{i,c} \leq z^g_i S^{gu}_{i,c}\ \ \forall i \in G,\forall c \in C \\
& \sum_{\substack{k\in G_i,c\in C}} S^g_{k,c} - \sum_{\substack{k\in L_i,c\in C}} z^d_k S^d_{k,c}- \sum_{\substack{k\in H_i,c\in C}} z^s_k Y^s_{k,c}\left | V_{i,c} \right |^2 \nonumber \\
& = \sum_{\substack{(i,j)\in E_i\cup E_i^R,c\in C}} S_{ij,c}\ \forall i \in N
\end{align}
```
