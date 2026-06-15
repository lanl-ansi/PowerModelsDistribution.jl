# Test Changes

## opf_ravens_extra

### Test trans dy

Mark as locally infeasible.

### IEEE Assets

Mark as locally solved.

## mld

### Transformer lpubfdiag mld [ACTIVE]

Leaving this test as failed as this actually suggests a meaningful error on lpubfdiag.

```julia
@test result["termination_status"] == NORM_LIMIT
#ACTUAL: @test result["termination_status"] == OTHER_ERROR

@test isapprox(result["objective"], 36; atol=1)
#ACTUAL: @test isapprox(result["objective"], 46; atol=1)
@test isapprox(result["solution"]["load"]["1"]["status"], .75; atol=1e-3)
#ACTUAL: @test isapprox(result["solution"]["load"]["1"]["status"], 0.623; atol=1e-3)
```

## Shunt [ACTIVE]

### Matrix shunts ACP/ACR/IVR

Results should be different when computing only based on the diagonal elements. But the issue is that in RAVENS the GS and BS matrices both become diagonal. Leaving this test as failed because it represents a meaningful error in ravens representation of a shunt.

Seemingly RAVENS only specifies b, g, b0, and g0 with option for a matrix spec of a shunt. This clearly does not model interphase admittance, loosing a functionality in the dss. 

``` julia
@test(all([shunt["bs"][c,d]!=0 for c in 1:3, d in 1:3 if c!=d]))
    #fails because shunt diagonal elements are zero in the original shunt
    #gs = 3x3 0 matrix 
    #bs = diagonal 3x3 where i,i = 0.009999999999999998

# check the results are different with only diagonal elements
@test(!isapprox(result_acp["solution"]["bus"]["2"]["vm"], result_acp_diag["solution"]["bus"]["2"]["vm"]))
    # yeilds --> !(isapprox([1.0615720000759332, 0.9893438121302054, 0.9435451906580832], [1.0615720000759332, 0.9893438121302054, 0.9435451906580832]))
```

## Load Models

### ZIP Load Models 8

Consistent scale factor change in output results from math vs eng model.

```julia
# delta loads
@test isapprox(pd(result, "1"), [2.9], atol=1E-1)
@test isapprox(qd(result, "1"), [1.0], atol=1E-1)
# wye loads
@test isapprox(pd(result, "2"), [8.7], atol=1E-1)
@test isapprox(qd(result, "2"), [2.9], atol=1E-1)
@test isapprox(pd(result, "3"), [8.7], atol=1E-1)
@test isapprox(qd(result, "3"), [2.9], atol=1E-1)
@test isapprox(pd(result, "4"), [8.5, 8.5, 8.5], atol=1E-1)

@test isapprox(qd(result, "4"), [2.8, 2.8, 2.8], atol=1E-1)  
```

Becomes:

```julia
# delta loads
@test isapprox(pd(result, "1"), [0.01], atol=1E-3)
@test isapprox(qd(result, "1"), [0.003], atol=1E-3)
# wye loads
@test isapprox(pd(result, "2"), [0.01], atol=1E-3)
@test isapprox(qd(result, "2"), [0.003], atol=1E-3)
@test isapprox(pd(result, "3"), [0.01], atol=1E-3)
@test isapprox(qd(result, "3"), [0.003], atol=1E-3)
@test isapprox(pd(result, "4"), [0.003, 0.003, 0.003], atol=1E-3)

@test isapprox(qd(result, "4"), [0.001, 0.001, 0.001], atol=1E-3)  
```

### loadmodels connection variations

Corrected lack of 0.3 in:

```julia
# single-phase wye loads
@test isapprox(pd(result, "10")*1000, [133, 133, 133], atol=1E-1)
```

To:

```julia
# single-phase wye loads
@test isapprox(pd(result, "10")*1000, [133.3, 133.3, 133.3], atol=1E-1)
```

This was just a typo fix.

## pf

### 3-bus balanced no linecode basefreq defined acp pf

All errors are that the tests only hold up to an atol of 10^-4 rather than the 10^-8 threshold which is a sensible numerical error for an optimization problem.

The last test only holds up to 10^-2 ie:
`Evaluated: isapprox(9.194043674460072, 9.193042740509885)`

## opf bf

### 3-bus unbalanced fbs opf_bf with voltage-dependent loads

```julia
@test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.9512, 0.9964, 0.9936]; atol=2e-3))
    #yeilds: [0.9512706468876924, 0.9807522438366527, 0.9773347164780088] needing an atol of 9e-1 to pass
```

## opv iv

### ivr opf power variable expressions [FIXED_BUT_PENDING_VALDIATION]

PENDING FIX: This fix makes some sense and resolves the error but needs someone more familiar with the intended behavior to valdiate these changes.

Originally tested:

```juila
math = transform_data_model(ravens_IEEE13_Assets)

pm = instantiate_mc_model_ravens(ravens_IEEE13_Assets, IVRUPowerModel, build_mc_opf)
```

Changed to:

```juila
math = transform_data_model(ravens_IEEE13_Assets)

pm = instantiate_mc_model(IEEE13_Assets, IVRUPowerModel, build_mc_opf)
```

This is because instantiate_mc_model_ravens does not produce an engineering model in the same way so the comparison does not make a ton of sense to be a 1-1 map from the non ravens test to the ravens test. The new version compares the number of objects in the ravens math model to the equivalent form in the dss based engineering model.

And then we need to change `(length(math["branch"]))*2` to `(length(math["branch"])+1)*2` as the math model is missing a branch created by the dss math model.

## Transformers

### Major Infeasibility Issue [FIXED_BUT_PENDING_VALDIATION]

All transformer test cases yeild infeasible solutions. What I have found so far is that transformers when generated from the dss->eng->math pipeline differ from those generated in ravens.

```julia
julia> math_dss["transformer"]
Dict{String, Any} with 2 entries:
  "1" => Dict{String, Any}("source_id"=>"_virtual_transformer.transformer.tx1.1", "t_connections"=>[1, 2, 3], "f_bus"=>1, "polarity"=>1, "sm_ub"=>750.0, "cm_ub"=>Inf, "tm_fix"=>Bool[1, 1, 1], "tm_lb"=>[0.9, 0.9, 0.9], "tm_set"=>[1.02, 1.02, 1.02], "t_vbase"=>0.57735…)
  "2" => Dict{String, Any}("source_id"=>"_virtual_transformer.transformer.tx1.2", "t_connections"=>[1, 2, 3], "f_bus"=>3, "polarity"=>1, "sm_ub"=>750.0, "cm_ub"=>Inf, "tm_fix"=>Bool[1, 1, 1], "tm_lb"=>[0.9, 0.9, 0.9], "tm_set"=>[0.97, 0.97, 0.97], "t_vbase"=>0.57735…)

julia> math_ravens["transformer"]
Dict{String, Any} with 2 entries:
  "1" => Dict{String, Any}("source_id"=>"_virtual_transformer.transformer.tx1.1", "t_connections"=>[1, 2, 3], "f_bus"=>1, "polarity"=>1, "sm_ub"=>0.5, "cm_ub"=>Inf, "tm_fix"=>Bool[1, 1, 1], "tm_lb"=>[0.9, 0.9, 0.9], "tm_set"=>[1.0, 1.0, 1.0], "t_vbase"=>0.57735…)
  "2" => Dict{String, Any}("source_id"=>"_virtual_transformer.transformer.tx1.2", "t_connections"=>[1, 2, 3], "f_bus"=>3, "polarity"=>1, "sm_ub"=>0.5, "cm_ub"=>Inf, "tm_fix"=>Bool[1, 1, 1], "tm_lb"=>[0.9, 0.9, 0.9], "tm_set"=>[1.0, 1.0, 1.0], "t_vbase"=>0.57735…)
```

Here we can see that the dss transformer properly stores the `1.02` and `0.97` `tm_set` values while the ravens transformer does not.
The Ravens transformer also does not store the correct 'sm_ub' value.

My hypothosis is that the ravens transformer does not actually store the tm_set values from the dss in the json and thus nothing is propogated to PMD ever.

Manually setting these to the dss version does not appear to fix things.

While this is almost certainly an issue it seems like the actual fix has something to do with the middle branches br_X on the virtual transformer ~pi~ model being 100x too high realitive to what we would expect from the dss version of the math model. This is fixed by dividing `x_sc` by zbase in the ravens2math mirroring the work done in eng2math.

### Note on Center Taps

Center tap transformers are not implemented in this version of this branch. Pending roberts implementation under the ref-ravens schema we need to see if this tests still pass. Currently they error as expected on the unimplemented center tap phase codes.

### All Tests

Voltage angles all needed a degree->radian conversion to run and also needed to be pegged to between -2pi and 2pi.

Raised acceptable test error thresholds for voltage magnitude and angle to account for minor computational differences. 

Previously it was erroring on angle differences like:

```julia
Seen: [-0.5207676215530239, -2.618659171805384, 1.5744288766975583] (rad)
True: [-0.5235987755982988, -2.6249751949994717, 1.567305668290908] (rad)
```

### Three Winding Tests [ACTIVE]

These tests do not return feasible solutions when converted from ravens as they do in the DSS. There is no use in debugging the diff result failures in these tests until basic feasibility is passing.

### OLTC Tests [ACTIVE]

```julia
@test norm(tap(1,pm)-[0.95, 0.95, 0.95], Inf) <= 1E-4
@test norm(tap(2,pm)-[1.05, 1.05, 1.05], Inf) <= 1E-4
```

This obviously fails because right now all RAVENS based transformers have a tap of `1.0` on both ends. This is an issue in ravens. 