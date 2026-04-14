const _LOGGER = Ref{Logging.ConsoleLogger}()

"""
    silence!()

Silence logging within PowerModelsDistribution.

This is equivalent to calling `set_logging_level!(:Error)`.
"""
silence!() = set_logging_level!(:Error)

"""
    reset_logging_level!()

Resets the log level to `:Info`.

This is equivalent to calling `set_logging_level!(:Info)`.
"""
reset_logging_level!() = set_logging_level!(:Info)

# A no-op. We keep this method for backward compatibility.
restore_global_logger!() = nothing

function _meta_formatter(level::Logging.LogLevel, _module, args...)
    return Logging.default_logcolor(level), "$(_module) | $level]:", ""
end

function set_logging_level!(level::Logging.LogLevel)
    _LOGGER[] =
        Logging.ConsoleLogger(stdout, level; meta_formatter = _meta_formatter)
    return
end

"""
    set_logging_level!(level::Symbol)

Set the logging level within PowerModelsDistribution.

`level` must be one of `:Error`, `:Warn`, `:Info`, or `:Debug`.
"""
set_logging_level!(level::Symbol) = set_logging_level!(getfield(Logging, level))

function _log_if_level(f, level, logger = _LOGGER[])
    if level >= Logging.min_enabled_level(logger)
        Logging.with_logger(f, logger)
    end
    return
end

# Currently unused.
# macro _error(msg)
#     return quote
#         $_log_if_level(() -> @error($msg), $(Logging.Error))
#         error($msg)
#     end |> esc
# end

macro _warn(msg)
    return :($_log_if_level(() -> @warn($msg), $(Logging.Warn)) )|> esc
end

macro _debug(msg)
    return :($_log_if_level(() -> @debug($msg), $(Logging.Debug))) |> esc
end

macro _info(msg)
    return :($_log_if_level(() -> @info($msg), $(Logging.Info)) )|> esc
end
