"""
    parse_file(file)

Parses a matlab .m `file` into a ThreePhase PowerModels data structure.
"""
function parse_file(file::String)
    if endswith(file, ".m")
        tppm_data = ThreePhasePowerModels.parse_matlab(file)
    else
        error(LOGGER, "only .m files are supported")
    end

    # TODO integrate network checks
    #check_network_data(tppm_data)

    return tppm_data
end