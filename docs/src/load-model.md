# Load Models


$S^d=[S^d_1...S^d_n]^T$ is a column vector $n\times 1$, specifying for each individual load how much power is consumed. By definition, this is

$$S^d=U^d\odot \left(I^d\right)^*,$$

where $U^d$ is the voltage applied across each individual load and $I^d$ is the current drawn by each load. At the same time, the composed load is connected to a bus with voltage $U^\text{bus}$, and draws a current $I^{\text{bus}}$ and power $S^{\text{bus}}$.
How these quantities relate to each other, depends on how the load is connected.

$$(S^d, U^d, I^d) ↔ (S^\text{bus}, U^\text{bus}, I^\text{bus})$$

In the implementations, these variables are referred to as


$$\begin{align}
S^d&=\text{pd}+j.\text{qd} & S^\text{bus}&=\text{pd_bus}+j.\text{qd_bus}\\
I^d&=\text{crd}+j.\text{cid} & I^\text{bus}&=\text{crd_bus}+j.\text{cid_bus}\\
U^d&=\text{vrd}+j.\text{vid} & U^\text{bus}&=\text{vr}+j.\text{vi}\\
\end{align}$$

## Voltage dependency
The general, exponential load model is defined as
$$P^d_i = P^{d,0}_i \left(\frac{V^d_i}{V^{d,0}_i}\right)^{\alpha_i} = a_i \left(V^d_i\right)^{\alpha_i}$$
$$Q^d_i = Q^{d,0}_i \left(\frac{V^d_i}{V^{d,0}_i}\right)^{\beta_i} = b_i \left(V^d_i\right)^{\beta_i}.$$

This might seem overly complicated, but occurs in distribution network data due to experimental model fitting of loads. There are a few cases which get a special name: constant power ($\alpha=\beta=0$), constant current ($\alpha=\beta=1$), and constant impedance ($\alpha=\beta=2$).

## Wye-connected Loads
A wye-connected load connects between a set of phases $\mathcal{P}$ and a neutral conductor $n$. The voltage as seen by each individual load is then
$$U^d = U^\text{bus}_\mathcal{P}-U^\text{bus}_n,$$
whilst the current
$$I^\text{bus}_\mathcal{P} = I^\text{d},\;\;\;I^\text{bus}_n=-1^TI^d$$
We now develop the expression for the power drawn at the bus for the phase conductors
$$
  S^\text{bus}_\mathcal{P} = (U^d+U^\text{bus}_n)\odot(I^d)^*
      = S^d+U^\text{bus}_n S^d\oslash U^d.
$$
From conservation of power or simply the formulas above,
$$
    S^\text{bus}_n = -1^TS^\text{bus}_\mathcal{P}+1^TS^d.
$$
### Grounded neutral
Note that when the neutral is grounded, i.e. $U^\text{bus}_n=0$, these formulas simplify to
$$S^\text{bus}_\mathcal{P}=S^d,\;\;\;S^\text{bus}_n=0,$$
which is why in Kron-reduced unbalanced networks, you can directly insert the power consumed by the loads, in the nodal power balance equations.

## Delta-connected Loads
Firstly, define the three-phase delta transformation matrix
$$M^\Delta_3 = \begin{bmatrix}\;\;\;1 & -1 & \;\;0\\ \;\;\;0 & \;\;\;1 & -1\\ -1 & \;\;\;0 & \;\;\;1\end{bmatrix},$$
which can be extended to more phases in a straight-forward manner. Now,
$$U^d = M^\Delta U^\text{bus},\;\;\; I^\text{bus} = \left(M^\Delta\right)^T I^d.$$
We can related $S^\text{bus}$ to $U^\text{bus}$ and $I^d$
$$
S^\text{bus} = U^\text{bus}\odot \left(I^\text{bus}\right)^*
             = U^\text{bus}\odot \left(M^\Delta\right)^T\left(I^d\right)^*,
$$
and using the fact that $\left(I^d\right)^*=S^d \oslash U^d$, and the expression above for $U^d$,
$$
S^\text{bus} = U^\text{bus}\left(M^\Delta\right)^T S^d \oslash M^\Delta U^\text{bus}
$$
