"MetaFormatter for ConsoleLogger for PMD"
function _pmd_metafmt(level::Logging.LogLevel, _module, group, id, file, line)
    @nospecialize
    color = Logging.default_logcolor(level)
    prefix = "$(_module) | " * (level == Logging.Warn ? "Warning" : string(level)) * " ] :"
    suffix = ""
    Logging.Info <= level < Logging.Warn && return color, prefix, suffix
    _module !== nothing && (suffix *= "$(_module)")
    if file !== nothing
        _module !== nothing && (suffix *= " ")
        suffix *= Base.contractuser(file)
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end
    !isempty(suffix) && (suffix = "@ " * suffix)

    return color, prefix, suffix
end


"Sets loglevel for PMD to :Error, silencing Info and Warn"
function silence!()
    set_logging_level!(:Error)
end


"Resets the log level to Info"
function reset_logging_level!()
    Logging.global_logger(_LOGGER)

    return
end


"Restores the global logger to its default state (before PMD was loaded)"
function restore_global_logger!()
    Logging.global_logger(_DEFAULT_LOGGER)

    return
end


"Sets the logging level for PMD: :Info, :Warn, :Error"
function set_logging_level!(level::Symbol)
    Logging.global_logger(_make_filtered_logger(getfield(Logging, level)))

    return
end


"Helper function to create the filtered logger for PMD"
function _make_filtered_logger(level)
    LoggingExtras.EarlyFilteredLogger(_LOGGER) do log
        if log._module == PowerModelsDistribution && log.level < level
            return false
        else
            return true
        end
    end
end
