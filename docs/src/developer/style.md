# Style Conventions

In general, the following conventions should be adhered to when making changes or additions to the code base. These conventions should include any conventions applied across the InfrastructureModels ecosystem specific to power engineering (_i.e_ conventions from InfrastructureModels, PowerModels, PowerModelsRestoration, etc.) with some additions specific to PowerModelsDistribution.

## Functions

Function additions should meeting the following criteria:

- All functions should be clearly named, without abbreviations, and with underscores between words, _e.g._ `parse_file` or `constraint_bus_voltage_magnitude`; in Python this is known as [`lower_case_with_underscores`](https://legacy.python.org/dev/peps/pep-0008/#descriptive-naming-styles). The exception to the abbreviate rule is cases where abbreviations would be expected in the modeling of power systems.
- All functions that are not prepended by an underscore `_` will be exported by default (_i.e._ when a user uses `using PowerModelsDistribution`). Public functions should have a detailed docstring instructing on usage
- All functions that modify data in place should end with an exclamation point `!` and the function input that is being modified should be the first argument (or first arguments in the case where multiple inputs are being modified in place). The exceptions to this rule are constraint and variable creation functions (_i.e._ those functions related to JuMP model creation), which do not include the exclamation point
- All function arguments, including keyword arguments, should have their types specified.
- Private functions, _i.e._ those intended to be for internal use only, should follow the same descriptive naming conventions as functions exported by default, and should always include docstrings to describe their purpose.
- Functions should be separated by two blank lines

```julia
"this function demonstrates how an internal, in-place data altering function should be defined"
function _concise_descriptive_name!(data::Dict{String,<:Any}, a::Real, b::Vector{<:Real}, c::Matrix{<:Complex}; d::Bool=false, e::Vector{Function}=Function[])
end
```

## Types & Enums

When specifying types, _i.e._ when specifying the type of a function argument, or creating enums, these guidelines are recommended:

- Prefer to use `Vector{T}` instead of `Array{T,1}`
- Prefer to use `Matrix{T}` instead of `Array{T,2}`
- Enums should __only__ be used in the `ENGINEERING` data model, never the `MATHEMATICAL` data model
- Enums must be added to the JSON parser when introduced

## Constants

Whenever possible, `const` should be used to eliminate unnecessary re-evaluations of code, and every `const` should have a docstring, whether internal or public.

## JuMP Variables and Constraints

For functions that create JuMP variables and constraints in particular, we follow the following naming convention as originally adopted by PowerModels:

`<jump macro id>(_<phase variant>)_<comp short name>_<quantity name>(_real|_imaginary|_magnitude|_angle|_factor)(_fr|_to)(_sqr)(_on_off)`

in the interest of intuitive names for users, the following special cases are also acceptable,

- `_power_real` -(can be replaced with)-> `_active`
- `_power_imaginary` -(can be replaced with)-> `_reactive`

In the case of PowerModelsDistribution, there are additional tags indicating that a function is a multiconductor variant, three-phase specific, etc.:

- `mc` multi-conductor, with an explicit neutral (which is the last conductor by convention?)
- `mp` multi-phase, for constraints that have no (explicit) neutral and multiple phases
- `3p` three-phase, when a constraint is hard-coded for three phases

Currently, all phase-aware functions use `mc`, but this is subject to change in the future as we refactor. If the function is not multiphase specific, these are not needed in the function name.

## Formulation Styles

- All new formulations should have __clear__ error messages when they do not support existing components. For example, if a formulation addition which is intended to work with OPF does not support delta-wye transformers, the `constraint_mc_transformer_power_dy`
- Formulation `abstract type` and `mutable struct` must be specified in [CapitalizedWords](https://legacy.python.org/dev/peps/pep-0008/#descriptive-naming-styles), which is a subtype of [camelCase](https://en.wikipedia.org/wiki/Camel_case) with the first word also capitalized.

## Problem Specification Styles

- If a new problem specification is only needed due to the requirements of a new formulation, and is not a new type of problem, _e.g._ another OPF formulation, a `build_` function with the same name as the existing formulation should be created that accepts a specific `PowerModel` (multiple dispatch)
- If a new problem specification is a new type of problem that will _e.g._ accept multiple formulations, new `build_` and `run_` functions should be created that do not collide with existing problem specification functions

## Metaprogramming

In general, it is better to avoid metaprogramming patterns, like creating functions algorithmically, in order to aid in the debugging of code. Metaprogramming can create significant challenges in interpreting stacktraces upon errors.

## Markdown

Markdown files should be properly formatted, particularly when including tables. Developers are encouraged to use [markdownlint](https://github.com/markdownlint/markdownlint) and a markdown formatter (such as in VSCode).

# File Structure

It is important that new functions, variables, constraints, etc. all go into appropriate places in the code base so that future maintenance and debugging is easier. Pay attention to the current file structure and attempt to conform as best as possible to it. In general

- `/src/core` contains the core logic of the package, including variable creation and constraint templates, _i.e._ things that are agnostic to the formulation
- `/src/data_model` contains all of the logic to transform between the `ENGINEERING` and `MATHEMATICAL` data models and model creation helper tools
- `src/form` contains formulation specific variable and constraint functions, organized under separate files for different formulations
- `src/io` contains all of the tools to parse and save files, in particular all of the logic necessary to parse dss files and output json files
- `src/prob` contains all problem specifications
- `docs/src` contains all source markdown files for the documentation
- `examples` contains Jupyter notebooks with walkthroughs of PowerModelsDistribution for new users

# Dependencies (Project.toml)

All new dependencies should be carefully considered before being added. It is important to keep the number of external dependencies low to avoid reliance on features that may not be maintained in the future. If possible, Julia Standard Library should be used, particularly in the case where reproducing the desired feature is trivial. There will be cases where it is not simple to duplicate a feature and subsequently maintain it within the package, so adding a dependency would be appropriate in such cases.

All new dependencies are are ultimately approved should also include an entry under `[compat]` indicating the acceptable versions (Julia automerge requirement). This includes test-only dependencies that appear under `[extras]`

The `Manifest.toml` __should not__ be included in the repo.
