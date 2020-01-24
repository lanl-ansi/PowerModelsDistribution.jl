import PowerModels.ref

ref(pm::_PMs.AbstractPowerModel, key::Symbol, idx, param::String, cnd::Int; nw::Int=pm.cnw) = pm.ref[:nw][nw][key][idx][param][cnd]
ref(pm::_PMs.AbstractPowerModel, nw::Int, key::Symbol, idx, param::String, cnd::Int) = pm.ref[:nw][nw][key][idx][param][cnd]
