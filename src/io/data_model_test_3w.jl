BASE_DIR = "/Users/sclaeys/code/PowerModelsDistribution.jl/src/io"
include("$BASE_DIR/data_model_util.jl")
include("$BASE_DIR/data_model_components.jl")
include("$BASE_DIR/../core/data_model_mapping.jl")
include("$BASE_DIR/../core/data_model_pu.jl")

## BUILD DATA MODEL

data_model = create_data_model()

add!(data_model, "linecode", create_linecode(id="3_conds", rs=LinearAlgebra.diagm(0=>fill(1.0, 3)), xs=LinearAlgebra.diagm(0=>fill(1.0, 3))))

add!(data_model, "voltage_source", create_voltage_source(id="source", bus="sourcebus", vm=[0.23, 0.23, 0.23], va=[0.0, -2*pi/3, 2*pi/3],
    #pg_max=fill(10E10,3), pg_min=fill(-10E10,3), qg_max=fill(10E10,3), qg_min=fill(-10E10,3),
))

add!(data_model, "line", create_line(id=:test, f_bus="sourcebus", t_bus="tr_prim", linecode="3_conds", length=1.3, f_connections=collect(1:3), t_connections=collect(1:3)))

add!(data_model, "bus", create_bus(id="sourcebus",  terminals=collect(1:3)))
add!(data_model, "bus", create_bus(id="tr_prim",    terminals=collect(1:4)))
add!(data_model, "bus", create_bus(id="tr_sec",     terminals=collect(1:4)))
#add!(data_model, "bus", create_bus("4", terminals=collect(1:4)))

# add!(data_model, "transformer_nw", create_transformer_nw("1", 3, ["2", "3", "4"], [[1,2,3], [1,2,3,4], [1,2,3]],
#     [0.230, 0.230, 0.230], [0.230, 0.230, 0.230],
#         configuration=["delta", "wye", "delta"],
#         xsc=[0.0, 0.0, 0.0],
#         rs=[0.0, 0.0, 0.0],
#         loadloss=0.00,
#         imag=0.00,
# ))

add!(data_model, "transformer_nw", create_transformer_nw(id="1", bus=["tr_prim", "tr_sec"], connections=[[1,2,3,4], [1,2,3,4]],
    vnom=[0.230, 100], snom=[0.230, 0.230],
    configuration=["wye", "wye"],
    xsc=[0.0],
    rs=[0.0, 0.0],
    noloadloss=0.00,
    imag=0.00,
))

add!(data_model, "load", create_load(id="1", bus="tr_prim", connections=collect(1:4), pd=[1.0, 1.0, 1.0]))

# add!(data_model, "generator", create_generator("1", "source",
#     connections=[1, 2, 3, 4],
#     pg_min=fill(-100, 3),
#     pg_max=fill( 100, 3),
#     qg_min=fill(-100, 3),
#     qg_max=fill( 100, 3),
# ))

#add!(data_model, "capacitor", create_capacitor(1, "tr_sec", 0.230*sqrt(3), qd_ref=[1.0, 1.0, 1.0]*1E-3, connections=[1,2,3], configuration="delta"))

check_data_model(data_model)
# PREP THE MODEL FOR OPTIMIZATION

data_model_map!(data_model)
#bsh = data_model["shunt"]["_virtual_1"]["b_sh"]
#
data_model_make_pu!(data_model, vbases=Dict("sourcebus"=>0.230))
#
data_model_index!(data_model)
data_model_make_compatible_v8!(data_model)

# OPTIMIZE THE MODEL

import PowerModelsDistribution
PMD = PowerModelsDistribution
import PowerModels
PMs = PowerModels
import InfrastructureModels
IM = InfrastructureModels

import JuMP, Ipopt

ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer)
pm = PMs.instantiate_model(data_model, PMs.ACPPowerModel, PMD.build_mc_opf, multiconductor=true, ref_extensions=[PMD.ref_add_arcs_trans!])
sol = PMs.optimize_model!(pm, optimizer=ipopt_solver)

solution_unmake_pu!(sol["solution"], data_model)
solution_identify!(sol["solution"], data_model)
solution_unmap!(sol["solution"], data_model)

sol["solution"]
