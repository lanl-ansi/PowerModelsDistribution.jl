using JSON
using Random
using UUIDs
using XLSX
using CSV, DataFrames
import Dates

rng = MersenneTwister(4689);

#primary workflow function:

function fix_cyme!(
    ravens_data::RavensModel;
    load_file::Union{String, Nothing} = nothing,
    solar_file::Union{String, Nothing} = nothing,
    meter_json::Union{String, Nothing} = nothing,
    split_data::Bool = false,
    )    
    """
    input: RAVENS object file that has been run though the native RAVENS cyme converter fixing:
        - PECs that were modeled as RotatingMachines
        - Impedance/Admittance matrices lacking rows/cols

    Parameters:
        - load_file: CSV of load profiles per meter
        - solar_file: CSV estimated solar generation per meter
        - meter_file: JSON meter lookup table

    Output: RAVENS object that resolves all further downstream issues needed to analyze RAVENS files including: 
        - Adds defined post-transformer loads to the RAVENS file 
        - Adds defined meter objects for fault analysis to the RAVENS file
    """
    if load_file !== nothing && solar_file !== nothing && meter_json !== nothing
        ravens_data = _generate_profiles!(ravens_data; load_file, solar_file, meter_json, split_data)
    end

    return ravens_data
end


#component workflow functions:

function _generate_profiles!(
    ravens_data::RavensModel;
    load_file::String,
    solar_file::String,
    meter_json::String,
    split_data::Bool,
    )


    # parse ravens object
    energy_consumers = ravens_data["EnergyConsumer"]
    PECs = ravens_data["PowerElectronicsConnection"]

    # Scale/multiplier for power values for load and PVs
    to_kw = 1000.0
    

    # Create Dictionary for BasicIntervalSchedule
    BasicIntervalSchedule = Dict()

    # Create Dictionary for Curve, if it does not exists in ravens
    if !haskey(ravens_data, "Curve")
        Curve = Dict()
    else
        Curve = ravens_data["Curve"]
    end

    # Read CSV file into a df & loop through all timesteps
    df_load = DataFrame(CSV.File(load_file))
    df_solar = DataFrame(CSV.File(solar_file))
    meters = JSON.parsefile(meter_json)


    transformer_load = Dict{String,Any}()
    transformer_gen = Dict{String,Any}()
    found = []

    #generate a dict from Transformer Names to Load Profiles at timestamps (meter name --> trans name)
    for meter in filter(!=("timestamp"), names(df_load))
        vals = df_load[!, meter]
        
        meter_info = get(meters, meter, nothing)
        if isnothing(meter_info)
            println("Meter $(meter) not found in model")
            continue
        end
        
        push!(found, meter)
        trans_name = meter_info["transformer"]
        meter_phase = meter_info["phase"][1]
        
        if haskey(transformer_load, trans_name)
            trans_data = transformer_load[trans_name]
            trans_data["vals"] += vals
            
            # Merge phases if not already present
            existing_phase = trans_data["phase"][1]
            if !occursin(meter_phase, existing_phase)
                trans_data["phase"][1] = existing_phase * meter_phase
            end
        else
            transformer_load[trans_name] = Dict{String,Any}(
                "vals" => vals,
                "phase" => copy(meter_info["phase"])
            )
        end
    end
    
    #generate a dict from Transformer Names to Solar Gen Profiles at timestamps (meter name --> trans name)
    for meter in filter(!=("timestamp"), names(df_solar))
        vals = df_solar[!, meter]
        
        meter_info = get(meters, meter, nothing)
        if isnothing(meter_info)
            println("Meter $(meter) not found in model")
            continue
        end
        
        push!(found, meter)
        trans_name = meter_info["transformer"]
        meter_phase = meter_info["phase"][1]
        
        if haskey(transformer_gen, trans_name)
            trans_data = transformer_gen[trans_name]
            trans_data["vals"] += vals
            
            # Merge phases if not already present
            existing_phase = trans_data["phase"][1]
            if !occursin(meter_phase, existing_phase)
                trans_data["phase"][1] = existing_phase * meter_phase
            end
        else
            transformer_gen[trans_name] = Dict{String,Any}(
                "vals" => vals,
                "phase" => copy(meter_info["phase"])
            )
        end
    end


    # #---------------------------------------------

    df_len = length(df_load[!, "timestamp"])


    # Loop through the entire timeseries
    for i in 1:1:df_len
        # P_sub_t = df[i, "Airport_Ckt_MW"] * 1000.0              # converted to kW
        # Q_sub_t = df[i, "Airport_Ckt_MVAR"] * (-1.0) * 1000.0   # converted to kVAR & negative due to excel file val error.
        # Irrad_t = df[i, "GEN_PITKIN_SOLAR,W3 (Average)"]/scale_solar

        # Load Forecasts
        for (name, ravens_obj) in energy_consumers
            _name = replace(name, "_L" => "")
            # Create EnergyConsumerSchedule dictionary if none exists
            if !haskey(BasicIntervalSchedule, "$(name)_Forecast")
                BasicIntervalSchedule["$(name)_Forecast"] = Dict()
                BasicIntervalSchedule["$(name)_Forecast"]["Ravens.cimObjectType"] = "EnergyConsumerSchedule"
                BasicIntervalSchedule["$(name)_Forecast"]["IdentifiedObject.mRID"] = "#_$(uppercase(string(UUIDs.uuid1(rng))))"
                BasicIntervalSchedule["$(name)_Forecast"]["IdentifiedObject.name"] = "$(name)_Forecast"
                BasicIntervalSchedule["$(name)_Forecast"]["BasicIntervalSchedule.value1Unit"] = "UnitSymbol.W"
                BasicIntervalSchedule["$(name)_Forecast"]["BasicIntervalSchedule.value1Multiplier"] = "UnitMultiplier.k"
                BasicIntervalSchedule["$(name)_Forecast"]["BasicIntervalSchedule.value2Unit"] = "UnitSymbol.VAr"
                BasicIntervalSchedule["$(name)_Forecast"]["BasicIntervalSchedule.value2Multiplier"] = "UnitMultiplier.k"
                BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.timeStep"] = 3600
                BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.startDay"] = "Day.Sunday"

                # Init vector
                BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.RegularTimePoints"] = Vector{Dict}(undef, df_len)

                # Add "EnergyConsumer.LoadProfile" to EnergyConsumer as reference
                ravens_obj["EnergyConsumer.LoadProfile"] = "EnergyConsumerSchedule::'$(name)_Forecast'"

            end


            # Compute P and Q Forecast at time t
            pf = .95
            if haskey(transformer_load, _name)
                P_calc = transformer_load[_name]["vals"][i] 
                Q_calc = transformer_load[_name]["vals"][i] ./ pf .* sin(acos(pf))
            else
                println("Load missing $(_name)")
                P_calc = 0.0
                Q_calc = 0.0
            end

            # 2) Compute load allocation
            BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.RegularTimePoints"][i] = Dict()
            BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.RegularTimePoints"][i]["RegularTimePoint.sequenceNumber"] = i

            
            BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.RegularTimePoints"][i]["RegularTimePoint.value1"] = P_calc
            BasicIntervalSchedule["$(name)_Forecast"]["EnergyConsumerSchedule.RegularTimePoints"][i]["RegularTimePoint.value2"] = Q_calc

            ###############  Correct the load values and use the max in the timeseries #########################
            if (i == 1) # Assing the first timeseries value to start
                ravens_obj["EnergyConsumer.p"] = P_calc*to_kw
                ravens_obj["EnergyConsumer.q"] = Q_calc*to_kw

                if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
                    num_phases = length(ravens_obj["EnergyConsumer.EnergyConsumerPhase"])
                    for p in 1:1:num_phases
                        ravens_obj["EnergyConsumer.EnergyConsumerPhase"][p]["EnergyConsumerPhase.p"] = ravens_obj["EnergyConsumer.p"]/num_phases
                        ravens_obj["EnergyConsumer.EnergyConsumerPhase"][p]["EnergyConsumerPhase.q"] = ravens_obj["EnergyConsumer.q"]/num_phases
                    end
                end
            else
                if (P_calc*to_kw > ravens_obj["EnergyConsumer.p"])
                    ravens_obj["EnergyConsumer.p"] = P_calc*to_kw

                    if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
                        num_phases = length(ravens_obj["EnergyConsumer.EnergyConsumerPhase"])
                        for p in 1:1:num_phases
                            ravens_obj["EnergyConsumer.EnergyConsumerPhase"][p]["EnergyConsumerPhase.p"] = ravens_obj["EnergyConsumer.p"]/num_phases
                        end
                    end

                end

                if (Q_calc*to_kw > ravens_obj["EnergyConsumer.q"])
                    ravens_obj["EnergyConsumer.q"] = Q_calc*to_kw

                    if haskey(ravens_obj, "EnergyConsumer.EnergyConsumerPhase")
                        num_phases = length(ravens_obj["EnergyConsumer.EnergyConsumerPhase"])
                        for p in 1:1:num_phases
                            ravens_obj["EnergyConsumer.EnergyConsumerPhase"][p]["EnergyConsumerPhase.q"] = ravens_obj["EnergyConsumer.q"]/num_phases
                        end
                    end
                end
            end
            ####################################

        end
        PECS_TO_IGNORE = []
        # Generation (Solar PV) Forecasts
        for (name, ravens_obj) in PECs
            if (name ∉ PECS_TO_IGNORE)
                _name = replace(name, "_G" => "")
                # Get rated S power from pec
                pecs_srated = ravens_obj["PowerElectronicsConnection.ratedS"]
                # Create values that go inside Curve
                if !haskey(Curve, "$(name)_Profile")
                    Curve["$(name)_Profile"] = Dict()
                    Curve["$(name)_Profile"]["Ravens.cimObjectType"] = "DispatchCurve"
                    Curve["$(name)_Profile"]["IdentifiedObject.mRID"] = "#_$(uppercase(string(UUIDs.uuid1(rng))))"
                    Curve["$(name)_Profile"]["IdentifiedObject.name"] = "$(name)_Profile"
                    Curve["$(name)_Profile"]["Curve.xUnit"] = "UnitSymbol.h"
                    Curve["$(name)_Profile"]["Curve.y1Unit"] = "UnitSymbol.W"
                    Curve["$(name)_Profile"]["Curve.y1Multiplier"] = "UnitMultiplier.k"

                    # Init vector
                    Curve["$(name)_Profile"]["Curve.CurveDatas"] = Vector{Dict}(undef, df_len)

                    # Add "EnergyConsumer.LoadProfile" to EnergyConsumer as reference
                    ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["PhotoVoltaicUnit.GenerationProfile"] = "Curve::'$(name)_Profile'"
                end
            
                # Compute dispatch
                Curve["$(name)_Profile"]["Curve.CurveDatas"][i] = Dict()
                Curve["$(name)_Profile"]["Curve.CurveDatas"][i]["CurveData.xvalue"] = i
                if haskey(transformer_gen, _name)
                    P_calc = transformer_gen[_name]["vals"][i] #.* 1000 
                else
                    println("Solar missing $(_name)")
                    P_calc = 0.0
                end
                Curve["$(name)_Profile"]["Curve.CurveDatas"][i]["CurveData.y1value"] = P_calc

                ###############  Correct the dispatch values and use the max in the timeseries #########################
                if (i == 1) # Assign the first timeseries value to start
                    ravens_obj["PowerElectronicsConnection.p"] = -P_calc*to_kw
                    if haskey(ravens_obj, "PowerElectronicsConnection.PowerElectronicsUnit")
                        ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["PowerElectronicsUnit.maxP"] = ravens_obj["PowerElectronicsConnection.ratedS"]
                    end
                else
                    if (P_calc*to_kw > (-1*ravens_obj["PowerElectronicsConnection.p"]))
                        ravens_obj["PowerElectronicsConnection.p"] = -P_calc*to_kw
                        if haskey(ravens_obj, "PowerElectronicsConnection.PowerElectronicsUnit")
                            ravens_obj["PowerElectronicsConnection.PowerElectronicsUnit"]["PowerElectronicsUnit.maxP"] = ravens_obj["PowerElectronicsConnection.ratedS"]
                        end
                    end
                end
                ####################################
            end

        end
    end
    if split_data
        return ravens_data, RavensModel(BasicIntervalSchedule), RavensModel(Curve)
    end
    ravens_data.data["BasicIntervalSchedule"] = BasicIntervalSchedule
    ravens_data.data["Curve"] = Curve
    return RavensModel(ravens_data.data)
end

