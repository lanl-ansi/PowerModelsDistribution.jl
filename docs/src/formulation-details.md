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
Applying all of the standard DC linearization tricks to the `StandardACPForm`

## `SOCWRPowerModel`
Applying the standard BIM voltage cross-product (sine and cosine) substitution tricks to `StandardACPForm` results immediately in a SOC formulation.


## `SDPUBFPowerModel`
```@docs
Order = [:type]
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
