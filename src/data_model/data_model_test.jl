using PowerModelsDistribution
import LinearAlgebra

function make_test_data_model()

    data_model = create_data_model()

    add_linecode!(data_model, "6_conds", rs=ones(6, 6), xs=ones(6, 6))
    add_linecode!(data_model, "4_conds", rs=ones(4, 4), xs=ones(4, 4))
    add_linecode!(data_model, "3_conds", rs=ones(3, 3), xs=ones(3, 3))
    add_linecode!(data_model, "2_conds", rs=ones(2, 2), xs=ones(2, 2))

    # 3 phase + 3 neutral conductors
    add_line!(data_model, "1",    f_bus="1", t_bus="2", linecode="6_conds", length=1, f_connections=[1,2,3,4,4,4], t_connections=collect(1:6))
    add_line!(data_model, "2",    f_bus="2", t_bus="3", linecode="6_conds", length=1, f_connections=[1,2,3,4,5,6], t_connections=[1,2,3,4,4,4])
    # 3 phase + 1 neutral conductors
    add_line!(data_model, "3",    f_bus="3", t_bus="4", linecode="4_conds", length=1.2)
    # 3 phase conductors
    add_line!(data_model, "4",    f_bus="4", t_bus="5", linecode="3_conds", length=1.3, f_connections=collect(1:3), t_connections=collect(1:3))
    # 2 phase + 1 neutral conductors
    add_line!(data_model, "5",    f_bus="4", t_bus="6", linecode="3_conds", length=1.3, f_connections=[1,3,4], t_connections=[1,3,4])
    # 1 phase + 1 neutral conductors
    add_line!(data_model, "6",    f_bus="4", t_bus="7", linecode="2_conds", length=1.7, f_connections=[2,4], t_connections=[2,4])
    # 2 phase conductors
    add_line!(data_model, "7",    f_bus="4", t_bus="8", linecode="2_conds", length=1.3, f_connections=[1,2], t_connections=[1,2])
    for i in 8:1000
    add_line!(data_model, id="$i",   f_bus="4", t_bus="8", linecode="2_conds", length=1.3, f_connections=[1,2], t_connections=[1,2])
    end

    add_bus!(data_model, "1",  terminals=collect(1:4))
    add_bus!(data_model, "2",  terminals=collect(1:6))
    add_bus!(data_model, "3",  terminals=collect(1:4))
    add_bus!(data_model, "4")
    add_bus!(data_model, "5",  terminals=collect(1:4))
    add_bus!(data_model, "6",  terminals=[1,3,4])
    add_bus!(data_model, "7",  terminals=[2,4])
    add_bus!(data_model, "8",  terminals=[1,2])
    add_bus!(data_model, "9",  terminals=[1,2,3,4])
    add_bus!(data_model, "10", terminals=[1,2,3])

    #
    add_load!(data_model, "1", bus="7", connections=[2,4], pd=[1.0], qd=[1.0])
    add_load!(data_model, "2", bus="8", connections=[1,2], pd_ref=[1.0], qd_ref=[1.0], model="constant_current", vnom=[230*sqrt(3)])
    add_load!(data_model, "3", bus="6", connections=[1,4], pd_ref=[1.0], qd_ref=[1.0], model="constant_impedance", vnom=[230])
    add_load!(data_model, "4", bus="6", connections=[3,4], pd_ref=[1.0], qd_ref=[1.0], model="exponential", vnom=[230], alpha=[1.2], beta=[1.5])
    add_load!(data_model, "5", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3))
    add_load!(data_model, "6", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_current", vnom=fill(230, 3))
    add_load!(data_model, "7", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_impedance", vnom=fill(230, 3))
    add_load!(data_model, "8", bus="4", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="exponential", vnom=fill(230, 3), alpha=[2.1,2.4,2.5], beta=[2.1,2.4,2.5])
    add_load!(data_model, "9", bus="5", configuration="delta", pd=fill(1.0, 3), qd=fill(1.0, 3))

    add_generator!(data_model, "1", bus="1", configuration="wye")

    add_transformer_nw!(data_model, "1", bus=["5", "9", "10"], connections=[[1,2,3], [1,2,3,4], [1,2,3]],
        vnom=[0.230, 0.230, 0.230], snom=[0.230, 0.230, 0.230],
        configuration=["delta", "wye", "delta"],
        xsc=[0.0, 0.0, 0.0],
        rs=[0.0, 0.0, 1.0],
        noloadloss=0.05,
        imag=0.05,
    )

    add_capacitor!(data_model, "cap_3ph",        bus="3", vnom=0.230*sqrt(3), qd_ref=[1, 2, 3])
    add_capacitor!(data_model, "cap_3ph_delta",  bus="4", vnom=0.230*sqrt(3), qd_ref=[1, 2, 3], configuration="delta", connections=[1,2,3])
    add_capacitor!(data_model, "cap_2ph_yg",     bus="6", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,3],  configuration="wye-grounded")
    add_capacitor!(data_model, "cap_2ph_yfl",    bus="6", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,3],  configuration="wye-floating")
    add_capacitor!(data_model, "cap_2ph_y",      bus="5", vnom=0.230*sqrt(3), qd_ref=[1, 2], connections=[1,3,4])

    return data_model
end


function make_3wire_data_model()

    data_model = create_data_model()

    add_linecode!(data_model, "3_conds", rs=LinearAlgebra.diagm(0=>fill(1.0, 3)), xs=LinearAlgebra.diagm(0=>fill(1.0, 3)))

    add_voltage_source!(data_model, "source", bus="sourcebus", connections=collect(1:3),
        vm=[0.23, 0.23, 0.23], va=[0.0, -2*pi/3, 2*pi/3],
        pg_max=fill(10E10,3), pg_min=fill(-10E10,3), qg_max=fill(10E10,3), qg_min=fill(-10E10,3),
        rs=ones(3,3)/10
    )

    # 3 phase conductors
    add_line!(data_model, id=:test, f_bus="sourcebus", t_bus="tr_prim", linecode="3_conds", length=1.3, f_connections=collect(1:3), t_connections=collect(1:3))

    add_bus!(data_model, "sourcebus",  terminals=collect(1:3))
    add_bus!(data_model, "tr_prim",    terminals=collect(1:4))
    add_bus!(data_model, "tr_sec",     terminals=collect(1:4))
    #add!(data_model, "bus", create_bus("4", terminals=collect(1:4)))

    # add!(data_model, "transformer_nw", create_transformer_nw("1", 3, ["2", "3", "4"], [[1,2,3], [1,2,3,4], [1,2,3]],
    #     [0.230, 0.230, 0.230], [0.230, 0.230, 0.230],
    #         configuration=["delta", "wye", "delta"],
    #         xsc=[0.0, 0.0, 0.0],
    #         rs=[0.0, 0.0, 0.0],
    #         loadloss=0.00,
    #         imag=0.00,
    # ))

    add_transformer_nw!(data_model, "1", bus=["tr_prim", "tr_sec"], connections=[[1,2,3,4], [1,2,3,4]],
        vnom=[0.230, 0.230], snom=[0.230, 0.230],
        configuration=["wye", "wye"],
        xsc=[0.0],
        rs=[0.0, 0.0],
        noloadloss=0.00,
        imag=0.00,
    )



    #
    add_load!(data_model, "1", bus="tr_prim", connections=collect(1:4), pd=[1.0, 1.0, 1.0])

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
data_model
#
data_model_map!(data_model)
#bsh = data_model["shunt"]["_virtual_1"]["b_sh"]
##
data_model_make_pu!(data_model, vbases=Dict("sourcebus"=>0.230))
#
data_model_index!(data_model)
data_model_make_compatible_v8!(data_model)
##
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

solution_make_si!(sol["solution"], data_model)
solution_identify!(sol["solution"], data_model)
solution_unmap!(sol["solution"], data_model)
#vm = sol["solution"]["bus"]["tr_sec"]["vm"]
#va = sol["solution"]["bus"]["tr_sec"]["va"]
#v = vm.*exp.(im*va/180*pi)
sol["solution"]
