# Formulation details

```@meta
CurrentModule = ThreePhasePowerModels
```

## `ACPPowerModel`
Real-valued formulation from (extended with shunts):
```
@INPROCEEDINGS{4153495,
author={B. Mahdad and T. Bouktir and K. Srairi},
booktitle={IECON 2006 - 32nd Annual Conference on IEEE Industrial Electronics},
title={A Three-Phase Power Flow Modelization: A Tool for Optimal Location and Control of FACTS Devices in Unbalanced Power Systems},
year={2006},
pages={2238-2243},
doi={10.1109/IECON.2006.347766},
ISSN={1553-572X},
month={Nov},}
```


## `DCPPowerModel`
This formulation applies all of the standard DC linearization tricks to the `StandardACPForm`:
- Taylor expansion of sine: $\sin(x) \approx x$
- Approximation of cosine: $\cos(x) \approx 1$
- Voltage magnitude approximated as one pu: $|V_{i,a}| \approx |V_{i,b}| \approx |V_{i,c}| \approx 1$pu.

## `SOCWRPowerModel`
This formulation applies the standard BIM voltage cross-product (sine and cosine) substitution tricks to the `StandardACPForm`, which results immediately in a SOC formulation.


## `SDPUBFPowerModel`
```@docs
SDPUBFPowerModel
```

## `SOCNLPUBFPowerModel`
```@docs
SOCNLPUBFPowerModel
```



## `SOCConicUBFPowerModel`
```@docs
SOCConicUBFPowerModel
```



## `LPfullUBFPowerModel`
```@docs
LPfullUBFPowerModel
```

## `LPdiagUBFPowerModel`
```@docs
LPdiagUBFPowerModel
```

## `LPLinUBFPowerModel`
```@docs
LPLinUBFPowerModel
```
