# PowerModelsDistribution exports everything except internal symbols, which are defined as
# those whose name starts with an underscore. If you don't want all of these
# symbols in your environment, then use `import PowerModelsDistribution` instead of
# `using PowerModelsDistribution`.

# Do not add PowerModelsDistribution-defined symbols to this exclude list. Instead, rename
# them with an underscore.

const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]
for sym in names(@__MODULE__, all=true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
         Base.isidentifier(sym_string[2:end])))
       continue
    end
    @eval export $sym
end

# the follow items are also exported for user-friendlyness when calling
# `using PowerModelsDistribution`

# so that users do not need to import JuMP to use a solver with PowerModelsDistribution
import JuMP: optimizer_with_attributes
export optimizer_with_attributes

import MathOptInterface: TerminationStatusCode
export TerminationStatusCode

import MathOptInterface: ResultStatusCode
export ResultStatusCode

for status_code_enum in [TerminationStatusCode, ResultStatusCode]
    for status_code in instances(status_code_enum)
        @eval import MathOptInterface: $(Symbol(status_code))
        @eval export $(Symbol(status_code))
    end
end

# InfrastructureModels exports
export nw_id_default, optimize_model!, ismultinetwork, update_data!
