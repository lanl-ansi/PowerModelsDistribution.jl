
# Shared patterns
- Each component has a unique (amonst components of the same type) identifier `id`.
- Everything is defined in SI units, except when a base is explicitly mentioned in the description.
- The default sometimes only mentions the default element, not the default value (e.g. `true` and not `fill(fill(true, 3), 3)`)

# Data model

## Bus

The data model below allows us to include buses of arbitrary many terminals (i.e., more than the usual four). This would be useful for
- underground lines with multiple neutrals which are not joined at every bus;
- distribution lines that carry several conventional lines in parallel (see for example the quad circuits in NEVTestCase).

name | default | type | description
------|----------|-----|-----------
terminals|[1,2,3,4]|Vector
vm_max | / | Vector | maximum conductor-to-ground voltage magnitude
vm_min | / | Vector | minimum conductor-to-ground voltage magnitude
vm_cd_max | / | Vector{Tuple} | e.g.  [(1,2,210)] means \|U1-U2\|>210
vm_cd_min | / | Vector{Tuple} | e.g.  [(1,2,230)] means \|U1-U2\|<230
grounded | [] | Vector | a list of terminals which are grounded
rg | [] | Vector | resistance of each defined grounding
xg | [] | Vector | reactance of each defined grounding

The tricky part is how to encode bounds for these type of buses. The most general is defining a list of three-tupples. Take for example a typical bus in a three-phase, four-wire network, where `terminals=[a,b,c,n]`. Such a bus might have
- phase-to-neutral bounds `vm_pn_max=250`, `vm_pn_min=210`
- and phase-to-phase bounds `vm_pp_max=440`, `vm_pp_min=360
`.

We can then define this equivalently as
- `vm_cd_max = [(a,n,250), (b,n,250), (c,n,250), (a,b,440), (b,c,440), (c,a,440)]`
- `vm_cd_min = [(a,n,210), (b,n,210), (c,n,210), (a,b,360), (b,c,360), (c,a,360)]`

Since this might be confusing for novice users, we also allow the user to define bounds through the following properties.


name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
phases | [1,2,3] | Vector
neutral | 4 |  | maximum conductor-to-ground voltage magnitude
vm_pn_max | / | Real | maximum phase-to-neutral voltage magnitude for all phases
vm_pn_min | / | Real | minimum phase-to-neutral voltage magnitude for all phases
vm_pp_max | / | Real | maximum phase-to-phase voltage magnitude for all phases
vm_pp_min | / | Real | minimum phase-to-phase voltage magnitude for all phases

## Line

This is a pi-model branch. A linecode implies `rs`, `xs`, `b_fr`, `b_to`, `g_fr` and `g_to`; when these properties are additionally specified, they overwrite the one supplied through the linecode.

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
f_bus | / | |
f_connections | / | | indicates for each conductor, to which terminal of the f_bus it connects
t_bus | / | |
t_connections | / | | indicates for each conductor, to which terminal of the t_bus it connects
linecode | / | | a linecode
rs | / | | series resistance matrix, size of n_conductors x n_conductors
xs | / | | series reactance matrix, size of n_conductors x n_conductors
g_fr | / | | from-side conductance
b_fr | / | | from-side susceptance
g_to | / | | to-side conductance
b_to | / | | to-side susceptance
c_rating | / | | symmetrically applicable current rating
s_rating | / | | symmetrically applicable power rating

## Linecode

- Should the linecode also include a `c_rating` and/`s_rating`?

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
rs | / | | series resistance matrix
xs | / | | series reactance matrix n_conductors
g_fr | / | | from-side conductance
b_fr | / | | from-side susceptance
g_to | / | | to-side conductance
b_to | / | | to-side susceptance

## Shunt

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
bus | / | |
connections | / | |
g | / | | conductance, size should be \|connections\|x\|connections\|
b | / | | susceptance, size should be \|connections\|x\|connections\|

## Capacitor

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
bus | / | |
connections | / | |
qd_ref | / | | conductance, size should be \|connections\|x\|connections\|
vnom | / | | conductance, size should be \|connections\|x\|connections\|

## Load

name | default | type | description
------|----------|----|-----------
id | / | | unique identifier
bus | / | |
connections | / | |
configuration | / | {wye, delta} | if wye-connected, the last connection will indicate the neutral
model | / | | indicates the type of voltage-dependency

### `model=constant_power`

name | default | type | description
------|----------|----|-----------
pd | / | Vector{Real} |
qd | / | Vector{Real} |

### `model=constant_current/impedance`

name | default | type | description
------|----------|----|-----------
pd_ref | / | Vector{Real} |
qd_ref | / | Vector{Real} |
vnom | / | Real |

### `model=exponential`

name | default | type | description
------|----------|----|-----------
pd_ref | / | Vector{Real} |
qd_ref | / | Vector{Real} |
vnom | / | Real |
exp_p | / | Vector{Real} |
exp_q | / | Vector{Real} |

## Generator

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
bus | / | |
connections | / | |
configuration | / | {wye, delta} | if wye-connected, the last connection will indicate the neutral
pg_min | / | | lower bound on active power generation per phase
pg_max | / | | upper bound on active power generation per phase
qg_min | / | | lower bound on reactive power generation per phase
qg_max | / | | upper bound on reactive power generation per phase

## AL2W Transformer
These are transformers are assymetric (A), lossless (L) and two-winding (2W). Assymetric refers to the fact that the secondary is always has a `wye` configuration, whilst the primary can be `delta`.

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
n_phases | size(rs)[1] | Int>0 | number of phases
f_bus | / | |
f_connections | / | | indicates for each conductor, to which terminal of the f_bus it connects
t_bus | / | |
t_connections | / | | indicates for each conductor, to which terminal of the t_bus it connects
configuration | / | {wye, delta} | for the from-side; the to-side is always connected in wye
tm_nom | / | Real | nominal tap ratio for the transformer
tm_max | / | Vector | maximum tap ratio for each phase (base=`tm_nom`)
tm_min | / | Vector | minimum tap ratio for each phase (base=`tm_nom`)
tm_set | fill(1.0, `n_phases`) | Vector | set tap ratio for each phase (base=`tm_nom`)
tm_fix | fill(true, `n_phases`) | Vector | indicates for each phase whether the tap ratio is fixed

TODO: add tm stuff

## Transformer
These are n-winding, n-phase, lossy transformers. Note that most properties are now Vectors (or Vectors of Vectors), indexed over the windings.

name | default | type | description
------|----------|-----|-----------
id | / | | unique identifier
n_phases | size(rs)[1] | Int>0 | number of phases
n_windings | size(rs)[1] | Int>0 | number of windings
bus | / | Vector | list of bus for each winding
connections | | Vector{Vector} | list of connection for each winding
configurations | | Vector{{wye, delta}} | list of configuration for each winding
xsc | 0.0 | Vector | list of short-circuit reactances between each pair of windings; enter as a list of the upper-triangle elements
rs | 0.0 | Vector | list of the winding resistance for each winding
tm_nom | / | Vector{Real} | nominal tap ratio for the transformer
tm_max | / | Vector{Vector} | maximum tap ratio for each winding and phase (base=`tm_nom`)
tm_min | / | Vector{Vector} | minimum tap ratio for for each winding and phase (base=`tm_nom`)
tm_set | 1.0 | Vector{Vector} | set tap ratio for each winding and phase (base=`tm_nom`)
tm_fix | true | Vector{Vector} | indicates for each winding and phase whether the tap ratio is fixed
