bus_name2id(tppm_data, name) = [bus["index"] for (_,bus) in tppm_data["bus"] if haskey(bus, "name") && bus["name"]==name][1]
load_name2id(tppm_data, name) = [load["index"] for (_,load) in tppm_data["load"] if haskey(load, "name") && load["name"]==name][1]
va(sol, tppm_data, name) = TPPMs.wraptopi(sol["solution"]["bus"][string(bus_name2id(tppm_data, name))]["va"][:])
vm(sol, tppm_data, name) = sol["solution"]["bus"][string(bus_name2id(tppm_data, name))]["vm"][:]
pdvar(pm, tppm_data, name) = [PMs.var(pm, pm.cnw, c, :pd, load_name2id(tppm_data, name)) for c in 1:3]
pd(pm, tppm_data, name) = [isa(x, Number) ? x : JuMP.value(x) for x in pdvar(pm, tppm_data, name)]
qdvar(pm, tppm_data, name) = [PMs.var(pm, pm.cnw, c, :qd, load_name2id(tppm_data, name)) for c in 1:3]
qd(pm, tppm_data, name) = [isa(x, Number) ? x : JuMP.value(x) for x in qdvar(pm, tppm_data, name)]
sd(pm, tppm_data, name) = pd(sol, tppm_data, name)+im*qd(sol, tppm_data, name)

@testset "loadmodels pf" begin
    @testset "connection variations" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_lm_1230.dss")
        pm = PMs.build_generic_model(tppm, PMs.ACPPowerModel, TPPMs.post_tp_pf_lm, multiconductor=true)
        sol = PMs.solve_generic_model(pm, ipopt_solver)
        # voltage magnitude at load bus
        @test isapprox(vm(sol, tppm, "loadbus"), [1, 1, 1], atol=1E-5)
        # single-phase delta loads
        @test isapprox(pd(pm, tppm, "d1ph"), [0.4000, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1ph"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1ph2"), [0, 0.4000, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1ph2"), [0, 0.3000, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1ph23"), [0, 0.2866, 0.1134], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1ph23"), [0, 0.0345, 0.2655], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1ph0"), [0, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1ph0"), [0, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1ph00"), [0, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1ph00"), [0, 0, 0], atol=1E-4)
        # Leads to an error now instead; order of conductors needed for this
        # @test isapprox(pd(pm, tppm, "d1ph123"), [0.4000, 0, 0], atol=1E-4)
        # @test isapprox(qd(pm, tppm, "d1ph123"), [0.3000, 0, 0], atol=1E-4)
        # single-phase wye loads
        # Leads to an error now instead; order of conductors needed for this
        # @test isapprox(pd(pm, tppm, "y1ph"), [0.4000, 0, 0], atol=1E-4)
        # @test isapprox(qd(pm, tppm, "y1ph"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1ph2"), [0, 0.4000, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1ph2"), [0, 0.3000, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1ph23"), [0, 0.2866, 0.1134], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1ph23"), [0, 0.0345, 0.2655], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1ph0"), [0, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1ph0"), [0, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1ph00"), [0, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1ph00"), [0, 0, 0], atol=1E-4)
        # Leads to an error now instead; order of conductors needed for this
        # @test isapprox(pd(pm, tppm, "y1ph123"), [0.4000, 0, 0], atol=1E-4)
        # @test isapprox(qd(pm, tppm, "y1ph123"), [0.3000, 0, 0], atol=1E-4)
        # three-phase loads
        @test isapprox(pd(pm, tppm, "d3ph"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3ph"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3ph123"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3ph123"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3ph1230"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3ph1230"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3ph213"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3ph213"), [0.100, 0.100, 0.100], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3ph3120"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3ph3120"), [0.100, 0.100, 0.100], atol=1E-4)
    end
    @testset "models 1/2/5 in acp pf" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_lm_models.dss")
        pm = PMs.build_generic_model(tppm, PMs.ACPPowerModel, TPPMs.post_tp_pf_lm, multiconductor=true)
        sol = PMs.solve_generic_model(pm, ipopt_solver)
        # voltage magnitude at load bus
        @test isapprox(vm(sol, tppm, "loadbus"), [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(pm, tppm, "d1phm1"), [0.4000, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm1"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1phm2"), [0.2783, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm2"), [0.2087, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1phm5"), [0.3336, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm5"), [0.2502, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm1"), [0.4000, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm1"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm2"), [0.2783, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm2"), [0.2087, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm5"), [0.3336, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm5"), [0.2502, 0, 0], atol=1E-4)
        # delta three-phase loads
        @test isapprox(pd(pm, tppm, "d3phm1"), [0.1160, 0.1465, 0.1375], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm1"), [0.0896, 0.0977, 0.1127], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3phm2"), [0.1005, 0.1348, 0.1212], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm2"), [0.0771, 0.0854, 0.1048], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3phm5"), [0.1080, 0.1405, 0.1291], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm5"), [0.0831, 0.0914, 0.1087], atol=1E-4)
        # wye three-phase loads
        @test isapprox(pd(pm, tppm, "y3phm1"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm1"), [0.1000, 0.1000, 0.1000], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y3phm2"), [0.0920, 0.1324, 0.1349], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm2"), [0.0690, 0.0993, 0.1012], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y3phm5"), [0.1108, 0.1329, 0.1341], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm5"), [0.0831, 0.0997, 0.1006], atol=1E-4)
    end
    @testset "models 1/2/5 in acr pf" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_lm_models.dss")
        pm = PMs.build_generic_model(tppm, PMs.ACRPowerModel, TPPMs.post_tp_pf_lm, multiconductor=true)
        sol = PMs.solve_generic_model(pm, ipopt_solver)
        # voltage magnitude at load bus
        @test isapprox(vm(sol, tppm, "loadbus"), [0.83072, 0.99653, 1.0059], atol=1.5E-4)
        # delta and wye single-phase load models
        @test isapprox(pd(pm, tppm, "d1phm1"), [0.4000, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm1"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1phm2"), [0.2783, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm2"), [0.2087, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d1phm5"), [0.3336, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d1phm5"), [0.2502, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm1"), [0.4000, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm1"), [0.3000, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm2"), [0.2783, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm2"), [0.2087, 0, 0], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y1phm5"), [0.3336, 0, 0], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y1phm5"), [0.2502, 0, 0], atol=1E-4)
        # delta three-phase loads
        @test isapprox(pd(pm, tppm, "d3phm1"), [0.1160, 0.1465, 0.1375], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm1"), [0.0896, 0.0977, 0.1127], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3phm2"), [0.1005, 0.1348, 0.1212], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm2"), [0.0771, 0.0854, 0.1048], atol=1E-4)
        @test isapprox(pd(pm, tppm, "d3phm5"), [0.1080, 0.1405, 0.1291], atol=1E-4)
        @test isapprox(qd(pm, tppm, "d3phm5"), [0.0831, 0.0914, 0.1087], atol=1E-4)
        # wye three-phase loads
        @test isapprox(pd(pm, tppm, "y3phm1"), [0.1333, 0.1333, 0.1333], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm1"), [0.1000, 0.1000, 0.1000], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y3phm2"), [0.0920, 0.1324, 0.1349], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm2"), [0.0690, 0.0993, 0.1012], atol=1E-4)
        @test isapprox(pd(pm, tppm, "y3phm5"), [0.1108, 0.1329, 0.1341], atol=1E-4)
        @test isapprox(qd(pm, tppm, "y3phm5"), [0.0831, 0.0997, 0.1006], atol=1E-4)
    end
end
