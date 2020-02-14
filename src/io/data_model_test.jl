BASE_DIR = "/Users/sclaeys/code/PowerModelsDistribution.jl/src/io"
include("$BASE_DIR/data_model_util.jl")
include("$BASE_DIR/data_model_components.jl")
include("$BASE_DIR/../core/data_model_mapping.jl")
include("$BASE_DIR/../core/data_model_pu.jl")

function make_test_data_model()

    data_model = create_data_model()

    add!(data_model, "linecode", create_linecode("6_conds", rs=ones(6, 6), xs=ones(6, 6)))
    add!(data_model, "linecode", create_linecode("4_conds", rs=ones(4, 4), xs=ones(4, 4)))
    add!(data_model, "linecode", create_linecode("3_conds", rs=ones(3, 3), xs=ones(3, 3)))
    add!(data_model, "linecode", create_linecode("2_conds", rs=ones(2, 2), xs=ones(2, 2)))

    # 3 phase + 3 neutral conductors
    add!(data_model, "line", create_line("1", "1", "2", "6_conds", 1; f_connections=[1,2,3,4,4,4], t_connections=collect(1:6)))
    add!(data_model, "line", create_line("2", "2", "3", "6_conds", 1; f_connections=[1,2,3,4,5,6], t_connections=[1,2,3,4,4,4]))
    # 3 phase + 1 neutral conductors
    add!(data_model, "line", create_line("3", "3", "4", "4_conds", 1.2))
    # 3 phase conductors
    add!(data_model, "line", create_line("4", "4", "5", "3_conds", 1.3; f_connections=collect(1:3), t_connections=collect(1:3)))
    # 2 phase + 1 neutral conductors
    add!(data_model, "line", create_line("5", "4", "6", "3_conds", 1.3, f_connections=[1,3,4], t_connections=[1,3,4]))
    # 1 phase + 1 neutral conductors
    add!(data_model, "line", create_line("6", "4", "7", "2_conds", 1.7, f_connections=[2,4], t_connections=[2,4]))
    # 2 phase conductors
    add!(data_model, "line", create_line("7", "4", "8", "2_conds", 1.3, f_connections=[1,2], t_connections=[1,2]))
    for i in 8:1000
        add!(data_model, "line", create_line("$i", "4", "8", "2_conds", 1.3, f_connections=[1,2], t_connections=[1,2]))
    end

    add!(data_model, "bus", create_bus("1", terminals=collect(1:4)))
    add!(data_model, "bus", create_bus("2", terminals=collect(1:6)))
    add!(data_model, "bus", create_bus("3", terminals=collect(1:4)))
    add!(data_model, "bus", create_bus("4"))
    add!(data_model, "bus", create_bus("5", terminals=collect(1:3)))
    add!(data_model, "bus", create_bus("6", terminals=[1,3,4]))
    add!(data_model, "bus", create_bus("7", terminals=[2,4]))
    add!(data_model, "bus", create_bus("8", terminals=[1,2]))
    add!(data_model, "bus", create_bus("9", terminals=[1,2,3,4]))
    add!(data_model, "bus", create_bus("10", terminals=[1,2,3]))

    #
    add!(data_model, "load", create_load("1", "7", connections=[2,4], pd=[1.0], qd=[1.0]))
    add!(data_model, "load", create_load("2", "8", connections=[1,2], pd_ref=[1.0], qd_ref=[1.0], model="constant_current", vnom=[230*sqrt(3)]))
    add!(data_model, "load", create_load("3", "6", connections=[1,4], pd_ref=[1.0], qd_ref=[1.0], model="constant_impedance", vnom=[230]))
    add!(data_model, "load", create_load("4", "6", connections=[3,4], pd_ref=[1.0], qd_ref=[1.0], model="exponential", vnom=[230], alpha=[1.2], beta=[1.5]))
    add!(data_model, "load", create_load("5", "4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3)))
    add!(data_model, "load", create_load("6", "4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_current", vnom=fill(230, 3)))
    add!(data_model, "load", create_load("7", "4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_impedance", vnom=fill(230, 3)))
    add!(data_model, "load", create_load("8", "4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="exponential", vnom=fill(230, 3), alpha=[2.1,2.4,2.5], beta=[2.1,2.4,2.5]))
    add!(data_model, "load", create_load("9", "5", configuration="delta", pd=fill(1.0, 3), qd=fill(1.0, 3)))

    add!(data_model, "generator", create_generator("1", "1", configuration="wye"))

    add!(data_model, "transformer_nw", create_transformer_nw("1", 3, ["5", "9", "10"], [[1,2,3], [1,2,3,4], [1,2,3]],
        [0.230, 0.230, 0.230], [0.230, 0.230, 0.230],
            configuration=["delta", "wye", "delta"],
            xsc=[0.0, 0.0, 0.0],
            rs=[0.0, 0.0, 1.0],
            loadloss=0.05,
            imag=0.05,
    ))

    add!(data_model, "capacitor", create_capacitor("cap_3ph", "3", 0.230*sqrt(3), qd_ref=[1, 2, 3]))
    add!(data_model, "capacitor", create_capacitor("cap_3ph_delta", "4", 0.230*sqrt(3), qd_ref=[1, 2, 3], configuration="delta", connections=[1,2,3]))
    add!(data_model, "capacitor", create_capacitor("cap_2ph_yg", "6", 0.230*sqrt(3), qd_ref=[1, 2], connections=[1,2],  configuration="wye-grounded"))
    add!(data_model, "capacitor", create_capacitor("cap_2ph_yfl", "6", 0.230*sqrt(3), qd_ref=[1, 2], connections=[1,2],  configuration="wye-floating"))
    add!(data_model, "capacitor", create_capacitor("cap_2ph_y", "5", 0.230*sqrt(3), qd_ref=[1, 2], connections=[1,3,4]))

    return data_model
end


function make_3wire_data_model()

    data_model = create_data_model()

    add!(data_model, "linecode", create_linecode("3_conds", rs=LinearAlgebra.diagm(0=>fill(1.0, 3)), xs=LinearAlgebra.diagm(0=>fill(1.0, 3))))

    # 3 phase conductors
    add!(data_model, "line", create_line(:test, "source", "tr_prim", "3_conds", 1.3; f_connections=collect(1:3), t_connections=collect(1:3)))

    add!(data_model, "bus", create_bus("source", terminals=collect(1:4)))
    add!(data_model, "bus", create_bus("tr_prim", terminals=collect(1:4)))
    add!(data_model, "bus", create_bus("tr_sec", terminals=collect(1:4)))
    #add!(data_model, "bus", create_bus("4", terminals=collect(1:4)))

    # add!(data_model, "transformer_nw", create_transformer_nw("1", 3, ["2", "3", "4"], [[1,2,3], [1,2,3,4], [1,2,3]],
    #     [0.230, 0.230, 0.230], [0.230, 0.230, 0.230],
    #         configuration=["delta", "wye", "delta"],
    #         xsc=[0.0, 0.0, 0.0],
    #         rs=[0.0, 0.0, 0.0],
    #         loadloss=0.00,
    #         imag=0.00,
    # ))

    add!(data_model, "transformer_nw", create_transformer_nw("1", 3, ["tr_prim", "tr_sec"], [[1,2,3,4], [1,2,3,4]],
        [0.230, 0.230], [0.230, 0.230],
            configuration=["wye", "wye"],
            xsc=[0.0],
            rs=[0.0, 0.0],
            loadloss=0.00,
            imag=0.00,
    ))



    #
    add!(data_model, "load", create_load("1", "tr_sec", connections=collect(1:4), pd=[1.0, 2.0, 3.0]/100, qd=[1.0, 2.0, 3.0]/100))

    add!(data_model, "generator", create_generator("1", "source",
        connections=[1, 2, 3, 4],
        pg_min=fill(-100, 3),
        pg_max=fill( 100, 3),
        qg_min=fill(-100, 3),
        qg_max=fill( 100, 3),
    ))

    return data_model
end


@time dm_hl = make_3wire_data_model()

@time check_data_model(dm_hl)
dm_hl
#
dm = map_down_data_model(dm_hl)
make_pu!(dm, vbases=Dict("source"=>230.0))
data_model_index!(dm)
##
import PowerModelsDistribution
PMD = PowerModelsDistribution
import PowerModels
PMs = PowerModels
import InfrastructureModels
IM = InfrastructureModels

import JuMP, Ipopt

ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer)

#dm_comp = PMD.parse_file("/Users/sclaeys/code/PowerModelsDistribution.jl/test/data/opendss/case3_balanced.dss")
dm = make_compatible_v8!(dm)
dm
##
pm = PMs.instantiate_model(dm, PMs.ACPPowerModel, PMD.build_mc_opf, multiconductor=true, ref_extensions=[PMD.ref_add_arcs_trans!])
sol = PMs.optimize_model!(pm, optimizer=ipopt_solver)

sol_remove_pu!(sol["solution"], dm)
solution_ind2id!(sol["solution"], dm)
