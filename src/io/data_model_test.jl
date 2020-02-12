BASE_DIR = "/Users/sclaeys/code/PowerModelsDistribution.jl/src/io"
include("$BASE_DIR/data_model_util.jl")
include("$BASE_DIR/data_model_components.jl")
include("$BASE_DIR/data_model_mapping.jl")
#include("$BASE_DIR/data_model_pu.jl")

function make_test_data_model()

    data_model = create_data_model()

    data_model["linecode"] = component_dict_from_list!([
        create_linecode("1", 6, rs=ones(6, 6), xs=ones(6, 6)),
        create_linecode("2", 4, rs=ones(4, 4), xs=ones(4, 4)),
        create_linecode("3", 3, rs=ones(3, 3), xs=ones(3, 3)),
        create_linecode("4", 2, rs=ones(2, 2), xs=ones(2, 2)),
    ])

    data_model["line"] = component_dict_from_list!([
        # 3 phase + 3 neutral conductors
        create_line("1", "1", "2", 6; t_terminals=[1,2,3,4,4,4]),
        # 3 phase + 1 neutral conductors
        create_line("2", "2", "3", 4),
        # 3 phase conductors
        create_line_with_linecode("3", "3", "4", data_model["linecode"]["4"], 1.3),
        # 2 phase + 1 neutral conductors
        create_line("4", "3", "5", 3, f_terminals=[1,3,4], t_terminals=[1,3,4]),
        # 1 phase + 1 neutral conductors
        create_line_with_linecode("5", "3", "6", data_model["linecode"]["4"], 1.7, f_terminals=[2,4], t_terminals=[2,4]),
        # 2 phase conductors
        create_line("6", "3", "7", 2, f_terminals=[1,2], t_terminals=[1,2]),
    ])

    data_model["voltage_zone"] = component_dict_from_list!([
        create_voltage_zone("1", 0.400)
    ])

    data_model["bus"] = component_dict_from_list!([
        create_bus_in_zone("1", "1", data_model, grounded=true, rg=0.1, xg=0.1, vm_fix=[fill(230, 3)..., 0], va_fix=[0, -pi*2/3, pi*2/3, 0]),
        create_bus_in_zone("2", "1", data_model),
        create_bus_in_zone("3", "1", data_model),
        create_bus_in_zone("4", "1", data_model, terminals=[1,2,3], neutral=missing, vm_ll_min=fill(0.9, 3), vm_ll_max=fill(1.1, 3)),
        create_bus_in_zone("5", "1", data_model, terminals=[1,3,4], phases=[1,3], vm_ln_min=fill(0.9, 2), vm_ln_max=fill(1.1, 2)),
        create_bus_in_zone("6", "1", data_model, terminals=[2,4], phases=[2], vm_ln_min=fill(0.9, 1), vm_ln_max=fill(1.1, 1)),
        create_bus_in_zone("7", "1", data_model, terminals=[1,2], phases=[1,2], neutral=missing, vm_ll_min=fill(0.9, 1), vm_ll_max=fill(1.1, 1)),
        create_bus_in_zone("8", "1", data_model, vm_ln_min=fill(230*0.9, 1), vm_ln_max=fill(1.1, 1)),
        create_bus_in_zone("9", "1", data_model, terminals=[1,2,3], neutral=missing, vm_ll_min=fill(0.9, 1), vm_ll_max=fill(1.1, 1)),
    ])

    data_model["load"] = component_dict_from_list!([
        create_load("1", "6", terminals=[2,4], pd=[1.0], qd=[1.0]),
        create_load("2", "7", terminals=[1,2], pd=[1.0], qd=[1.0], model="constant_current", vm_nom=[230*sqrt(3)]),
        create_load("3", "5", terminals=[1,4], pd=[1.0], qd=[1.0], model="constant_impedance", vm_nom=[230]),
        create_load("4", "5", terminals=[3,4], pd=[1.0], qd=[1.0], model="exponential", vm_nom=[230], alpha=[1.2], beta=[1.5]),
        create_load("5", "3", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3)),
        create_load("6", "3", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_current", vm_nom=fill(230, 3)),
        create_load("7", "3", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="constant_impedance", vm_nom=fill(230, 3)),
        create_load("8", "3", configuration="wye", pd=fill(1.0, 3), qd=fill(1.0, 3), model="exponential", vm_nom=fill(230, 3), alpha=[2.1,2.4,2.5], beta=[2.1,2.4,2.5]),
        create_load("9", "4", configuration="delta", pd=fill(1.0, 3), qd=fill(1.0, 3)),
    ])

    data_model["generator"] = component_dict_from_list!([
        create_generator("1", "1", configuration="wye"),
    ])

    data_model["transformer_nw3ph_lossy"] = component_dict_from_list!([
        create_transformer_nw3ph_lossy("1", 3, ["4", "8", "9"], [[1,2,3], [1,2,3,4], [1,2,3]],
        [0.230, 0.230, 0.230], [0.230, 0.230, 0.230],
            configuration=["delta", "wye", "delta"],
            xsc=[0.0, 0.0, 0.0],
            rs=[0.0, 0.0, 1.0],
            loadloss=0.05,
            imag=0.05,
        ),
    ])

    data_model["capacitor"] = component_dict_from_list!([
        create_capacitor("cap_3ph", "3", 0.230*sqrt(3), qd_ref=[1, 2, 3]),
        create_capacitor("cap_3ph_delta", "4", 0.230*sqrt(3), qd_ref=[1, 2, 3], configuration="delta", terminals=[1,2,3]),
        create_capacitor("cap_2ph_yg", "6", 0.230*sqrt(3), qd_ref=[1, 2], terminals=[1,2],  configuration="wye-grounded"),
        create_capacitor("cap_2ph_yfl", "6", 0.230*sqrt(3), qd_ref=[1, 2], terminals=[1,2],  configuration="wye-floating"),
        create_capacitor("cap_2ph_y", "5", 0.230*sqrt(3), qd_ref=[1, 2], terminals=[1,3,4]),
    ])

    return data_model
end

#dm = make_test_data_model()
#index_data_model(dm, components=["capacitor"])
