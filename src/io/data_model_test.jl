BASE_DIR = "/Users/sclaeys/code/PowerModelsDistribution.jl/src/io"
include("$BASE_DIR/data_model_util.jl")
include("$BASE_DIR/data_model_components.jl")
include("$BASE_DIR/../core/data_model_mapping.jl")
include("$BASE_DIR/../core/data_model_pu.jl")

function make_test_data_model()

    data_model = create_data_model()

    add!(data_model, "linecode", create_linecode(id="6_conds", rs=ones(6, 6), xs=ones(6, 6)))
    add!(data_model, "linecode", create_linecode(id="4_conds", rs=ones(4, 4), xs=ones(4, 4)))
    add!(data_model, "linecode", create_linecode(id="3_conds", rs=ones(3, 3), xs=ones(3, 3)))
    add!(data_model, "linecode", create_linecode(id="2_conds", rs=ones(2, 2), xs=ones(2, 2)))

    # 3 phase + 3 neutral conductors
    add!(data_model, "line", create_line(id="1",    f_bus="1", t_bus="2", linecode="6_conds", length=1, f_connections=[1,2,3,4,4,4], t_connections=collect(1:6)))
    add!(data_model, "line", create_line(id="2",    f_bus="2", t_bus="3", linecode="6_conds", length=1, f_connections=[1,2,3,4,5,6], t_connections=[1,2,3,4,4,4]))
    # 3 phase + 1 neutral conductors
    add!(data_model, "line", create_line(id="3",    f_bus="3", t_bus="4", linecode="4_conds", length=1.2))
    # 3 phase conductors
    add!(data_model, "line", create_line(id="4",    f_bus="4", t_bus="5", linecode="3_conds", length=1.3, f_connections=collect(1:3), t_connections=collect(1:3)))
    # 2 phase + 1 neutral conductors
    add!(data_model, "line", create_line(id="5",    f_bus="4", t_bus="6", linecode="3_conds", length=1.3, f_connections=[1,3,4], t_connections=[1,3,4]))
    # 1 phase + 1 neutral conductors
    add!(data_model, "line", create_line(id="6",    f_bus="4", t_bus="7", linecode="2_conds", length=1.7, f_connections=[2,4], t_connections=[2,4]))
    # 2 phase conductors
    add!(data_model, "line", create_line(id="7",    f_bus="4", t_bus="8", linecode="2_conds", length=1.3, f_connections=[1,2], t_connections=[1,2]))
    for i in 8:1000
    add!(data_model, "line", create_line(id="$i",   f_bus="4", t_bus="8", linecode="2_conds", length=1.3, f_connections=[1,2], t_connections=[1,2]))
    end

    add!(data_model, "bus", create_bus(id="1",  terminals=collect(1:4)))
    add!(data_model, "bus", create_bus(id="2",  terminals=collect(1:6)))
    add!(data_model, "bus", create_bus(id="3",  terminals=collect(1:4)))
    add!(data_model, "bus", create_bus(id="4"))
    add!(data_model, "bus", create_bus(id="5",  terminals=collect(1:3)))
    add!(data_model, "bus", create_bus(id="6",  terminals=[1,3,4]))
    add!(data_model, "bus", create_bus(id="7",  terminals=[2,4]))
    add!(data_model, "bus", create_bus(id="8",  terminals=[1,2]))
    add!(data_model, "bus", create_bus(id="9",  terminals=[1,2,3,4]))
    add!(data_model, "bus", create_bus(id="10", terminals=[1,2,3]))

    #
    add!(data_model, "load", create_load(id="1", bus="7", connections=[2,4], pd=[1.0], qd=[1.0]))
    add!(data_model, "load", create_load(id="2", bus="8", connections=[1,2], pd_ref=[1.0], qd_ref=[1.0], model="constant_current", vnom=[230*sqrt(3)]))
    add!(data_model, "load", create_load(id="3", bus="6", connections=[1,4], pd_ref=[1.0], qd_ref=[1.0], model="constant_impedance", vnom=[230]))
    add!(data_model, "load", create_load(id="4", bus="6", connections=[3,4], pd_ref=[1.0], qd_ref=[1.0], model="exponential", vnom=[230], alpha=[1.2], beta=[1.5]))
    add!(data_model, "load", create_load(id="5", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3)))
    add!(data_model, "load", create_load(id="6", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_current", vnom=fill(230, 3)))
    add!(data_model, "load", create_load(id="7", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_impedance", vnom=fill(230, 3)))
    add!(data_model, "load", create_load(id="8", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="exponential", vnom=fill(230, 3), alpha=[2.1,2.4,2.5], beta=[2.1,2.4,2.5]))
    add!(data_model, "load", create_load(id="9", bus="5", configuration="delta", pd=fill(1.0, 3), qd=fill(1.0, 3)))

    add!(data_model, "generator", create_generator(id="1", bus="1", configuration="wye"))

    add!(data_model, "transformer_nw", create_transformer_nw(id="1", bus=["5", "9", "10"], connections=[[1,2,3], [1,2,3,4], [1,2,3]],
        vnom=[0.230, 0.230, 0.230], snom=[0.230, 0.230, 0.230],
        configuration=["delta", "wye", "delta"],
        xsc=[0.0, 0.0, 0.0],
        rs=[0.0, 0.0, 1.0],
        loadloss=0.05,
        imag=0.05,
    ))

    add!(data_model, "capacitor", create_capacitor(id="cap_3ph",        bus="3", vnom=0.230*sqrt(3), qd_ref=[1, 2, 3]))
    add!(data_model, "capacitor", create_capacitor(id="cap_3ph_delta",  bus="4", vnom=0.230*sqrt(3), qd_ref=[1, 2, 3], configuration="delta", connections=[1,2,3]))
    add!(data_model, "capacitor", create_capacitor(id="cap_2ph_yg",     bus="6", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,2],  configuration="wye-grounded"))
    add!(data_model, "capacitor", create_capacitor(id="cap_2ph_yfl",    bus="6", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,2],  configuration="wye-floating"))
    add!(data_model, "capacitor", create_capacitor(id="cap_2ph_y",      bus="5", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,3,4]))

    return data_model
end


function make_3wire_data_model()

    data_model = create_data_model()

    add!(data_model, "linecode", create_linecode(id="3_conds", rs=LinearAlgebra.diagm(0=>fill(1.0, 3)), xs=LinearAlgebra.diagm(0=>fill(1.0, 3))))

    add!(data_model, "voltage_source", create_voltage_source(id="source", bus="sourcebus", vm=[0.23, 0.23, 0.23], va=[0.0, -2*pi/3, 2*pi/3],
        #pg_max=fill(10E10,3), pg_min=fill(-10E10,3), qg_max=fill(10E10,3), qg_min=fill(-10E10,3),
    ))

    # 3 phase conductors
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
        vnom=[0.230, 0.230], snom=[0.230, 0.230],
        configuration=["wye", "wye"],
        xsc=[0.0],
        rs=[0.0, 0.0],
        noloadloss=0.00,
        imag=0.00,
    ))



    #
    add!(data_model, "load", create_load(id="1", bus="tr_prim", connections=collect(1:4), pd=[1.0, 1.0, 1.0]))

    # add!(data_model, "generator", create_generator("1", "source",
    #     connections=[1, 2, 3, 4],
    #     pg_min=fill(-100, 3),
    #     pg_max=fill( 100, 3),
    #     qg_min=fill(-100, 3),
    #     qg_max=fill( 100, 3),
    # ))

    #add!(data_model, "capacitor", create_capacitor(1, "tr_sec", 0.230*sqrt(3), qd_ref=[1.0, 1.0, 1.0]*1E-3, connections=[1,2,3], configuration="delta"))

    return data_model
end


data_model = make_3wire_data_model()
check_data_model(data_model)

data_model_map!(data_model)
#bsh = data_model["shunt"]["_virtual_1"]["b_sh"]
#
data_model_make_pu!(data_model, vbases=Dict("sourcebus"=>0.230))
#
data_model_index!(data_model)
data_model_make_compatible_v8!(data_model)
#
import PowerModelsDistribution
PMD = PowerModelsDistribution
import PowerModels
PMs = PowerModels
import InfrastructureModels
IM = InfrastructureModels

import JuMP, Ipopt

ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer)
pm = PMs.instantiate_model(data_model, PMs.IVRPowerModel, PMD.build_mc_opf_iv, multiconductor=true, ref_extensions=[PMD.ref_add_arcs_trans!])
sol = PMs.optimize_model!(pm, optimizer=ipopt_solver)

solution_unmake_pu!(sol["solution"], data_model)
solution_identify!(sol["solution"], data_model)
solution_unmap!(sol["solution"], data_model)
#vm = sol["solution"]["bus"]["tr_sec"]["vm"]
#va = sol["solution"]["bus"]["tr_sec"]["va"]
#v = vm.*exp.(im*va/180*pi)
sol["solution"]
