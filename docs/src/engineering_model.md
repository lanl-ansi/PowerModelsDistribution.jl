
# Introduction to the PowerModelsDistribution Data Models

In this notebook we introduce the engineering data model added to PowerModelsDistribution in version v0.9.0. We will give serveral examples of how to use this new data model directly, including new transformations that have become easier with it's introduction, how to convert it to the the lower-level mathematical model that was previously the only user interface we offered, and how to get various types of results using this new model.

## Imports

All commands in this document with no package namespace specified are directly exported by PowerModelsDistribution or already available in Julia base. Any commands that are only avaiable via an external package will be specified by including by using `import`, which will require specifying the originating package before the command, _e.g._ `Ipopt.Optimizer` as you will see below.


```julia
using PowerModelsDistribution
```

    ┌ Info: Precompiling PowerModelsDistribution [d7431456-977f-11e9-2de3-97ff7677985e]
    └ @ Base loading.jl:1260


In these examples we will use the following optimization solvers, specified using `optimizer_with_attributes` from JuMP v0.21


```julia
import Ipopt

ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)
```




    MathOptInterface.OptimizerWithAttributes(Ipopt.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute,Any}[MathOptInterface.RawParameter("tol") => 1.0e-6, MathOptInterface.RawParameter("print_level") => 0])



## Parsing Data

Here we give the first example of how to parse data into the `ENGINEERING` data model structure, which is the default data structure type that the user will see without passing additional arguments, as we demonstrate later.

We start with a 3 bus unbalanced load case provided as a dss file in the `test` folder of the PowerModelsDistribution.jl repository


```julia
eng = parse_file("../test/data/opendss/case3_unbalanced.dss")
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "case3_unbalanced.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "case3_unbalanced.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 37 in "case3_unbalanced.dss" is not supported, skipping.[39m





    Dict{String,Any} with 9 entries:
      "voltage_source" => Dict{Any,Any}("source"=>Dict{String,Any}("source_id"=>"vs…
      "name"           => "3bus_example"
      "line"           => Dict{Any,Any}("quad"=>Dict{String,Any}("cm_ub"=>[400.0, 4…
      "settings"       => Dict{String,Any}("sbase_default"=>500.0,"vbases_default"=…
      "files"          => ["case3_unbalanced.dss"]
      "load"           => Dict{Any,Any}("l2"=>Dict{String,Any}("source_id"=>"load.l…
      "bus"            => Dict{Any,Any}("primary"=>Dict{String,Any}("rg"=>Float64[]…
      "linecode"       => Dict{Any,Any}("556mcm"=>Dict{String,Any}("b_fr"=>[25.4648…
      "data_model"     => ENGINEERING



Different information and warning messages will be given depending on the input file. In the case above, these messages all related to various parse notifications that arise during a parse of a dss file, and can be safely ignored

The resulting data structure is a Julia dictionary. The first notable field is `"data_model"` which specifies which data model this data structure corresponds to, in this case `ENGINEERING`. This value is expected to be an `Enum` of type `DataModel`

The next notable field is `"settings"`, which contains some important default/starting values for the distribution network


```julia
eng["settings"]
```




    Dict{String,Any} with 5 entries:
      "sbase_default"        => 500.0
      "vbases_default"       => Dict("sourcebus"=>0.23094)
      "voltage_scale_factor" => 1000.0
      "power_scale_factor"   => 1000.0
      "base_frequency"       => 50.0



- `"sbase_default"` is the starting value for the power base,
- `"vbases_default"` is the starting voltage base for the case, and multiple voltage bases can be specified, which would be useful in cases where there are multiple isolated islands with their own generation,
- `"voltage_scale_factor"` is a scaling factor for all voltage values, which in the case of OpenDSS is in kV by default
- `"power_scale_factor"` is a scaling factor for all power values
- `"base_frequency"` is the base frequency of the network in Hz, which is useful to know for mixed frequency networks

Next we look at the `"bus"` components


```julia
eng["bus"]
```




    Dict{Any,Any} with 3 entries:
      "primary"   => Dict{String,Any}("rg"=>Float64[],"grounded"=>Int64[],"status"=…
      "sourcebus" => Dict{String,Any}("rg"=>Float64[],"grounded"=>Int64[],"status"=…
      "loadbus"   => Dict{String,Any}("rg"=>[0.0],"grounded"=>[4],"status"=>ENABLED…



We can see there are three buses in this system, identified by ids `"primary"`, `"sourcebus"`, and `"loadbus"`. 

__NOTE__: In Julia, order of Dictionary keys is not fixed, nor does it retain the order in which it was parsed like _e.g._ `Vectors`. 

Identifying components by non-integer names is a new feature of the `ENGINEERING` model, and makes network debugging more straightforward. 

__NOTE__: all names are converted to lowercase on parse from the originating dss file.

Each bus component has the following properties in the `ENGINEERING` model


```julia
eng["bus"]["sourcebus"]
```




    Dict{String,Any} with 5 entries:
      "rg"        => Float64[]
      "grounded"  => Int64[]
      "status"    => ENABLED
      "terminals" => [1, 2, 3]
      "xg"        => Float64[]



- `"terminals"` indicates which terminals on the bus have active connections
- `"grounded"` indicates which terminals are grounded
- `"rg"` and `"xg"` indicate the grounding resistance and reactance of the ground
- `"status"` indicates whether a bus is `ENABLED` or `DISABLED`, and is specified for every component in the engineering model

Next, we look at the `"line"` components, which is a generic name for both overhead lines and underground cables, which we do not differentiate between in the nomenclature


```julia
eng["line"]
```




    Dict{Any,Any} with 2 entries:
      "quad"   => Dict{String,Any}("cm_ub"=>[400.0, 400.0, 400.0],"cm_ub_c"=>[600.0…
      "ohline" => Dict{String,Any}("cm_ub"=>[400.0, 400.0, 400.0],"cm_ub_c"=>[600.0…




```julia
eng["line"]["quad"]
```




    Dict{String,Any} with 11 entries:
      "cm_ub"         => [400.0, 400.0, 400.0]
      "cm_ub_c"       => [600.0, 600.0, 600.0]
      "f_connections" => [1, 2, 3]
      "length"        => 1.0
      "status"        => ENABLED
      "source_id"     => "line.quad"
      "t_connections" => [1, 2, 3]
      "f_bus"         => "primary"
      "t_bus"         => "loadbus"
      "cm_ub_b"       => [600.0, 600.0, 600.0]
      "linecode"      => "4/0quad"



Again, we see components identified by their OpenDSS names. A `"line"` is an edge object, which will always have the following properties:

- `"f_bus"`
- `"t_bus"`
- `"f_connections"` - list of terminals to which the line is connected on the from-side
- `"t_connections"` - list of terminals to which the line is connected on the to-side

Here we are also introduced to two important concepts, the `"source_id"`, which is an easy way to identify from where an object originates in the dss file, and a data type element, pointed to by `"linecode"` in this case.

A data type element is an element that does not represent a real engineering object, but only contains data that one of those real objects can refer to, in this case a linecode, which contains information like line resistance/reactance and conductance/susceptance.


```julia
eng["linecode"]["4/0quad"]
```




    Dict{String,Any} with 6 entries:
      "b_fr" => [25.4648 -0.0 -0.0; -0.0 25.4648 -0.0; -0.0 -0.0 25.4648]
      "rs"   => [0.1167 0.0467 0.0467; 0.0467 0.1167 0.0467; 0.0467 0.0467 0.1167]
      "xs"   => [0.0667 0.0267 0.0267; 0.0267 0.0667 0.0267; 0.0267 0.0267 0.0667]
      "b_to" => [25.4648 -0.0 -0.0; -0.0 25.4648 -0.0; -0.0 -0.0 25.4648]
      "g_to" => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
      "g_fr" => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]



Next, we introduce a node element, the `"load"` object, where we also see the first example of a specification of less than three phases at a time


```julia
eng["load"]["l1"]
```




    Dict{String,Any} with 10 entries:
      "source_id"     => "load.l1"
      "qd_nom"        => [3.0]
      "status"        => ENABLED
      "model"         => POWER
      "connections"   => [1, 4]
      "vm_nom"        => 0.23094
      "pd_nom"        => [9.0]
      "dispatchable"  => NO
      "bus"           => "loadbus"
      "configuration" => WYE



We can see that the length of the Vectors for `"pd_nom"` and `"qd_nom"` are only one, although the number of terminals listed in `"connections"` is two. This is because the connection is WYE, and therefore the final connection is a grounded neutral

Here we are also introduced to two new Enums, `WYE`, which gives the connection configuration, and `NO` under dispatchable, which indicates that if this case were used in an MLD problem, _i.e._ with `run_mc_mld` that this load would not be sheddable.

Finally, we show the generation source for this case, which in opendss is a voltage source named `"source"`


```julia
eng["voltage_source"]["source"]
```




    Dict{String,Any} with 8 entries:
      "source_id"   => "vsource.source"
      "rs"          => [4.27691e-8 3.96342e-9 3.96342e-9; 3.96342e-9 4.27691e-8 3.9…
      "va"          => [0.0, -120.0, 120.0]
      "status"      => ENABLED
      "connections" => [1, 2, 3]
      "vm"          => [0.229993, 0.229993, 0.229993]
      "xs"          => [1.54178e-7 -1.04497e-9 -1.04497e-9; -1.04497e-9 1.54178e-7 …
      "bus"         => "sourcebus"



- `"vm"` - specifies the fixed voltage magnitudes per phase at the bus
- `"va"` - specifies the fixed reference angles per phases at the bus
- `"rs"` and `"xs"` specifies internal impedances of the voltage source

### Importing raw dss properties

In case there are additional properties that you want to use from dss, it is possible to import those directly into the `ENGINEERING` (and `MATHEMATICAL`) data structure with the `import_all` keyword argument


```julia
eng_all = parse_file("../test/data/opendss/case3_unbalanced.dss"; import_all=true)

eng_all["line"]
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "case3_unbalanced.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "case3_unbalanced.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 37 in "case3_unbalanced.dss" is not supported, skipping.[39m





    Dict{Any,Any} with 2 entries:
      "quad"   => Dict{String,Any}("cm_ub"=>[400.0, 400.0, 400.0],"cm_ub_c"=>[600.0…
      "ohline" => Dict{String,Any}("cm_ub"=>[400.0, 400.0, 400.0],"cm_ub_c"=>[600.0…



You will note the presence of `"dss"` dictionaries under components, and `"dss_options"` at the root level

## Running Optimal Power Flow

In this section we introduce how to run an optimal power flow (opf) in PowerModelsDistribution on an engineering data model

In order to run an OPF problem you will need

1. a data model
2. a formulation
3. a solver

In these examples we will use the `eng` model we worked with above, the `ACPPowerModel`, which is a AC power flow formulation in polar coordinates, and the `ipopt_solver` we already defined above


```julia
result = run_mc_opf(eng, ACPPowerModel, ipopt_solver)
```

    [35m[warn | PowerModels]: Updated generator 1 cost function with order 3 to a function of order 2: [0.5, 0.0][39m
    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit http://projects.coin-or.org/Ipopt
    ******************************************************************************
    





    Dict{String,Any} with 8 entries:
      "solve_time"         => 3.9884
      "optimizer"          => "Ipopt"
      "termination_status" => LOCALLY_SOLVED
      "dual_status"        => FEASIBLE_POINT
      "primal_status"      => FEASIBLE_POINT
      "objective"          => 0.0214812
      "solution"           => Dict{String,Any}("voltage_source"=>Dict{Any,Any}("sou…
      "objective_lb"       => -Inf



The result of `run_mc_opf` will be very familiar to those who are already familiar with PowerModels and PowerModelsDistribution. The notable difference will be in the `"solution"` dictionary


```julia
result["solution"]
```




    Dict{String,Any} with 6 entries:
      "voltage_source" => Dict{Any,Any}("source"=>Dict{String,Any}("qg_bus"=>[3.177…
      "line"           => Dict{Any,Any}("quad"=>Dict{String,Any}("qf"=>[3.09488, 3.…
      "settings"       => Dict{String,Any}("sbase"=>0.5)
      "load"           => Dict{Any,Any}("l2"=>Dict{String,Any}("qd_bus"=>[0.0, 3.0,…
      "bus"            => Dict{Any,Any}("primary"=>Dict{String,Any}("va"=>[-0.22425…
      "per_unit"       => false



Here you can see that the solution comes back out by default into the same data model as is provided by the user to the `run_` command, as well as being in SI units, as opposed to per unit, which is used during the solve. For example,


```julia
result["solution"]["bus"]["loadbus"]
```




    Dict{String,Any} with 2 entries:
      "va" => [-0.484238, -120.243, 120.274]
      "vm" => [0.222521, 0.226727, 0.225577]



If for some reason you want to return the result in per-unit rather than SI, you can specify this in the `run_` command by


```julia
result_pu = run_mc_opf(eng, ACPPowerModel, ipopt_solver; make_si=false)

result_pu["solution"]["bus"]["loadbus"]
```

    [35m[warn | PowerModels]: Updated generator 1 cost function with order 3 to a function of order 2: [0.5, 0.0][39m





    Dict{String,Any} with 2 entries:
      "va" => [-0.484238, -120.243, 120.274]
      "vm" => [0.963546, 0.981757, 0.976779]



### Branch Flow formulations

Previously, to use a branch flow formulation, such as `SOCNLPUBFPowerModel`, it was required to use a different `run_` command, but now, by using multiple dispatch we have simplified this for the user


```julia
result_bf = run_mc_opf(eng, SOCNLPUBFPowerModel, ipopt_solver)
```

    [35m[warn | PowerModels]: Updated generator 1 cost function with order 3 to a function of order 2: [0.5, 0.0][39m





    Dict{String,Any} with 8 entries:
      "solve_time"         => 0.190441
      "optimizer"          => "Ipopt"
      "termination_status" => LOCALLY_SOLVED
      "dual_status"        => FEASIBLE_POINT
      "primal_status"      => FEASIBLE_POINT
      "objective"          => 0.0211303
      "solution"           => Dict{String,Any}("voltage_source"=>Dict{Any,Any}("sou…
      "objective_lb"       => -Inf



## Engineering Model Transformations

One of the power things about the engineering model is that data transformations are much more simple. Here we illustrate two examples that are currently included in PowerModelsDistribution, but writing your own data transformation functions will be trivial, as we will show

First, there are several objects that have loss models by default when parsing from dss files, such as voltage sources, transformers, and switches. To remove these loss models, therefore making these components lossless, we can use the included `make_lossess!` function. Here we use a basic 2-winding wye-wye connected transformer case from `test` to illustrate this


```julia
eng_ut = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")

eng_ut["transformer"]["tx1"]
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "ut_trans_2w_yy.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 38 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m





    Dict{String,Any} with 16 entries:
      "polarity"      => [1, 1]
      "sm_nom"        => [500.0, 500.0]
      "tm_lb"         => [[0.9, 0.9, 0.9], [0.9, 0.9, 0.9]]
      "connections"   => [[1, 2, 3, 4], [1, 2, 3, 4]]
      "tm_set"        => [[1.02, 1.02, 1.02], [0.97, 0.97, 0.97]]
      "tm_step"       => [[0.03125, 0.03125, 0.03125], [0.03125, 0.03125, 0.03125]]
      "bus"           => ["1", "2"]
      "configuration" => ConnConfig[WYE, WYE]
      "noloadloss"    => 0.05
      "xsc"           => [0.05]
      "source_id"     => "transformer.tx1"
      "rw"            => [0.01, 0.02]
      "tm_fix"        => Array{Bool,1}[[1, 1, 1], [1, 1, 1]]
      "vm_nom"        => [11.0, 4.0]
      "tm_ub"         => [[1.1, 1.1, 1.1], [1.1, 1.1, 1.1]]
      "imag"          => 0.11



We can see that `"noloadloss"`, `"rw"`, and `"imag"` are non-zero, but if we apply the `make_lossless!` function we can see these parameters are set to zero, effectively eliminating the losses


```julia
make_lossless!(eng_ut)

eng_ut["transformer"]["tx1"]
```




    Dict{String,Any} with 16 entries:
      "polarity"      => [1, 1]
      "sm_nom"        => [500.0, 500.0]
      "tm_lb"         => [[0.9, 0.9, 0.9], [0.9, 0.9, 0.9]]
      "connections"   => [[1, 2, 3, 4], [1, 2, 3, 4]]
      "tm_set"        => [[1.02, 1.02, 1.02], [0.97, 0.97, 0.97]]
      "tm_step"       => [[0.03125, 0.03125, 0.03125], [0.03125, 0.03125, 0.03125]]
      "bus"           => ["1", "2"]
      "configuration" => ConnConfig[WYE, WYE]
      "noloadloss"    => 0.0
      "xsc"           => [0.0]
      "source_id"     => "transformer.tx1"
      "rw"            => [0.0, 0.0]
      "tm_fix"        => Array{Bool,1}[[1, 1, 1], [1, 1, 1]]
      "vm_nom"        => [11.0, 4.0]
      "tm_ub"         => [[1.1, 1.1, 1.1], [1.1, 1.1, 1.1]]
      "imag"          => 0.0



Alternatively, we can apply this function at parse


```julia
eng_ut = parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[make_lossless!])

eng_ut["transformer"]["tx1"]
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "ut_trans_2w_yy.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 38 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m





    Dict{String,Any} with 16 entries:
      "polarity"      => [1, 1]
      "sm_nom"        => [500.0, 500.0]
      "tm_lb"         => [[0.9, 0.9, 0.9], [0.9, 0.9, 0.9]]
      "connections"   => [[1, 2, 3, 4], [1, 2, 3, 4]]
      "tm_set"        => [[1.02, 1.02, 1.02], [0.97, 0.97, 0.97]]
      "tm_step"       => [[0.03125, 0.03125, 0.03125], [0.03125, 0.03125, 0.03125]]
      "bus"           => ["1", "2"]
      "configuration" => ConnConfig[WYE, WYE]
      "noloadloss"    => 0.0
      "xsc"           => [0.0]
      "source_id"     => "transformer.tx1"
      "rw"            => [0.0, 0.0]
      "tm_fix"        => Array{Bool,1}[[1, 1, 1], [1, 1, 1]]
      "vm_nom"        => [11.0, 4.0]
      "tm_ub"         => [[1.1, 1.1, 1.1], [1.1, 1.1, 1.1]]
      "imag"          => 0.0



Another transformation function included in PowerModelsDistribution is the `apply_voltage_bounds!` function, which will apply some voltage bounds in SI units, given some percent value, _e.g._ if we want the lower bound on voltage to be `0.9` and upper bound `1.1` after per-unit conversion


```julia
apply_voltage_bounds!(eng_ut; vm_lb=0.9, vm_ub=1.1)

eng_ut["bus"]["2"]
```




    Dict{String,Any} with 7 entries:
      "rg"        => [0.0]
      "grounded"  => [4]
      "status"    => ENABLED
      "terminals" => [1, 2, 3, 4]
      "vm_ub"     => [2.54034, 2.54034, 2.54034, 2.54034]
      "vm_lb"     => [2.07846, 2.07846, 2.07846, 2.07846]
      "xg"        => [0.0]



Alternatively, this can be specified at parse by


```julia
eng_ut = parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[make_lossless!, (apply_voltage_bounds!, "vm_lb"=>0.9, "vm_ub"=>1.1)])

eng_ut["bus"]["2"]
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "ut_trans_2w_yy.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 38 in "ut_trans_2w_yy.dss" is not supported, skipping.[39m





    Dict{String,Any} with 7 entries:
      "rg"        => [0.0]
      "grounded"  => [4]
      "status"    => ENABLED
      "terminals" => [1, 2, 3, 4]
      "vm_ub"     => [2.54034, 2.54034, 2.54034, 2.54034]
      "vm_lb"     => [2.07846, 2.07846, 2.07846, 2.07846]
      "xg"        => [0.0]



## Mathematical Model

In this section we introduce the mathematical model, which was the previous user-facing model in PowerModelsDistribution, explain how conversions between the model happen in practice, and give an example of how to do this conversion manually

In practice, unless the user is interested, the conversion between the `ENGINEERING` and `MATHEMATICAL` models should be seemless and invisible to the user. By providing an `ENGINEERING` model to a `run_` command the `run_mc_model` command will know to convert the model to `MATHEMATICAL`, which will be used to the generate the JuMP model that will actually be optimized. Similarly, the solution generated by this optimization will be automatically converted back to the format of the `ENGINEERING` model.

Let's first take a look at how to convert to the `MATHEMATICAL` model


```julia
math = transform_data_model(eng)
```

    [35m[warn | PowerModels]: Updated generator 1 cost function with order 3 to a function of order 2: [0.5, 0.0][39m





    Dict{String,Any} with 18 entries:
      "bus"         => Dict{String,Any}("4"=>Dict{String,Any}("grounded"=>Bool[0, 0…
      "name"        => "3bus_example"
      "dcline"      => Dict{String,Any}()
      "map"         => Dict{String,Any}[Dict("unmap_function"=>"_map_math2eng_root!…
      "settings"    => Dict{String,Any}("sbase_default"=>500.0,"vbases_default"=>Di…
      "gen"         => Dict{String,Any}("1"=>Dict{String,Any}("pg"=>[0.0, 0.0, 0.0]…
      "branch"      => Dict{String,Any}("1"=>Dict{String,Any}("br_r"=>[1.09406 0.43…
      "storage"     => Dict{String,Any}()
      "switch"      => Dict{String,Any}()
      "basekv"      => 0.23094
      "baseMVA"     => 0.5
      "conductors"  => 3
      "per_unit"    => true
      "data_model"  => MATHEMATICAL
      "shunt"       => Dict{String,Any}()
      "transformer" => Dict{String,Any}()
      "bus_lookup"  => Dict{Any,Int64}("primary"=>1,"sourcebus"=>2,"loadbus"=>3)
      "load"        => Dict{String,Any}("1"=>Dict{String,Any}("model"=>POWER,"conne…



There are a couple of things to notice right away. First, the data model transform automatically converts the model to per-unit. Second, there are a lot of empty component sets, whereas in the `ENGINEERING` model, only component types that had components in them were listed. In the `MATHEMATICAL` model certain component dictionaries are always expected to exist, and the `eng2math` conversion functions will automatically populate these.

Next, there are a few unusal fields, such as `"settings"`, which previously didn't exist in the `MATHEMATICAL` model. This is used for the per-unit conversion specifically in PowerModelsDistribution. Also, is the `"map"` field, which is a `Vector` of Dictionaries that enable the conversion back to `ENGINEERING` from `MATHEMATICAL`. Without this it would be impossible to convert back, and in fact only the solution can be converted, because some properties are combined destructively during the conversion to the `MATHEMATICAL` model, and therefore cannot be reverse engineered. However, since the conversion to `MATHEMATICAL` is not in-place, you will always have a copy of `eng` alongside `math`.

Here is an example of one of the `"map"` entries


```julia
math["map"][end]
```




    Dict{String,Any} with 3 entries:
      "to"             => ["gen.1", "bus.4", "branch.3"]
      "from"           => "source"
      "unmap_function" => "_map_math2eng_voltage_source!"



Alternatively, the `MATHEMATICAL` model can be returned directly from the `parse_file` command with the `data_model` keyword argument


```julia
math = parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model=MATHEMATICAL)
```

    [32m[info | PowerModels]: Circuit has been reset with the "clear" on line 1 in "case3_unbalanced.dss"[39m
    [35m[warn | PowerModels]: Command "calcvoltagebases" on line 35 in "case3_unbalanced.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Command "solve" on line 37 in "case3_unbalanced.dss" is not supported, skipping.[39m
    [35m[warn | PowerModels]: Updated generator 1 cost function with order 3 to a function of order 2: [0.5, 0.0][39m





    Dict{String,Any} with 18 entries:
      "bus"         => Dict{String,Any}("4"=>Dict{String,Any}("grounded"=>Bool[0, 0…
      "name"        => "3bus_example"
      "dcline"      => Dict{String,Any}()
      "map"         => Dict{String,Any}[Dict("unmap_function"=>"_map_math2eng_root!…
      "settings"    => Dict{String,Any}("sbase_default"=>500.0,"vbases_default"=>Di…
      "gen"         => Dict{String,Any}("1"=>Dict{String,Any}("pg"=>[0.0, 0.0, 0.0]…
      "branch"      => Dict{String,Any}("1"=>Dict{String,Any}("br_r"=>[1.09406 0.43…
      "storage"     => Dict{String,Any}()
      "switch"      => Dict{String,Any}()
      "basekv"      => 0.23094
      "baseMVA"     => 0.5
      "conductors"  => 3
      "per_unit"    => true
      "data_model"  => MATHEMATICAL
      "shunt"       => Dict{String,Any}()
      "transformer" => Dict{String,Any}()
      "bus_lookup"  => Dict{Any,Int64}("primary"=>1,"sourcebus"=>2,"loadbus"=>3)
      "load"        => Dict{String,Any}("1"=>Dict{String,Any}("model"=>POWER,"conne…



### Running `MATHEMATICAL` models

There is very little difference from the user point-of-view in running `MATHEMATICAL` models other than the results will not be automatically converted back to the the format of the `ENGINEERING` model


```julia
result_math = run_mc_opf(math, ACPPowerModel, ipopt_solver)

result_math["solution"]
```




    Dict{String,Any} with 7 entries:
      "baseMVA"    => 0.5
      "branch"     => Dict{String,Any}("1"=>Dict{String,Any}("qf"=>[0.00618975, 0.0…
      "gen"        => Dict{String,Any}("1"=>Dict{String,Any}("qg_bus"=>[0.00635515,…
      "load"       => Dict{String,Any}("1"=>Dict{String,Any}("qd_bus"=>[0.0, 0.006,…
      "conductors" => 3
      "bus"        => Dict{String,Any}("4"=>Dict{String,Any}("va"=>[0.0, -2.0944, 2…
      "per_unit"   => true



It is also possible to manually convert the solution back to the `ENGINEERING` format, provided you have the __map__


```julia
result_eng = transform_solution(result_math["solution"], math)
```




    Dict{String,Any} with 6 entries:
      "voltage_source" => Dict{Any,Any}("source"=>Dict{String,Any}("qg_bus"=>[3.177…
      "line"           => Dict{Any,Any}("quad"=>Dict{String,Any}("qf"=>[3.09488, 3.…
      "settings"       => Dict{String,Any}("sbase"=>0.5)
      "load"           => Dict{Any,Any}("l2"=>Dict{String,Any}("qd_bus"=>[0.0, 3.0,…
      "bus"            => Dict{Any,Any}("primary"=>Dict{String,Any}("va"=>[-0.22425…
      "per_unit"       => false



## Conclusion

This concludes the introduction to the `ENGINEERING` data model and conversion to the `MATHEMATICAL` model. We hope that you will find this new data model abstraction easy to use and simple to understand
