function _resolve_sparlectra_casefile(casefile::String, path::Union{Nothing,String})::String
  ext = lowercase(splitext(casefile)[2])
  ext in (".m", ".jl") || throw(ArgumentError("run_sparlectra: file extension $(ext) is not supported; use .m or .jl."))
  filename = path === nothing ? joinpath(pwd(), "data", "mpower", strip(casefile)) : joinpath(path, strip(casefile))
  isfile(filename) || error("File $(filename) not found")
  return filename
end

function _import_sparlectra_net(casefile::String, path::Union{Nothing,String}, cfg::SparlectraConfig; performance_profile = nothing)::Net
  filename = _resolve_sparlectra_casefile(casefile, path)
  pf_cfg = cfg.powerflow
  mat_cfg = cfg.matpower
  mpc = _perf_profile_time!(performance_profile, :matpower_case_parse) do
    MatpowerIO.read_case(filename; legacy_compat = true)
  end
  net = _perf_profile_time!(performance_profile, :network_construction) do
    createNetFromMatPowerCase(
      mpc = mpc,
      flatstart = pf_cfg.start_mode.flatstart,
      cooldown = pf_cfg.qlimits.cooldown_iters,
      q_hyst_pu = pf_cfg.qlimits.hysteresis_pu,
      enable_pq_gen_controllers = mat_cfg.enable_pq_gen_controllers,
      bus_shunt_model = mat_cfg.bus_shunt_model,
      matpower_shift_sign = mat_cfg.shift_sign,
      matpower_shift_unit = mat_cfg.shift_unit,
      matpower_ratio = mat_cfg.ratio,
      matpower_pv_voltage_source = mat_cfg.pv_voltage_source,
      matpower_pv_voltage_mismatch_tol_pu = mat_cfg.pv_voltage_mismatch_tol_pu,
      preallocate_network = mat_cfg.preallocate_network,
      preallocate_min_buses = mat_cfg.preallocate_min_buses,
      profile = performance_profile,
    )
  end
  MatpowerIO.apply_mp_isolated_buses!(net, mpc; verbose = 0)
  MatpowerIO.apply_mp_bus_vmva_init!(net, mpc; flatstart = pf_cfg.start_mode.flatstart, verbose = 0)
  if pf_cfg.start_mode.flatstart && (pf_cfg.start_mode.voltage_mode != :classic || pf_cfg.start_mode.angle_mode != :classic)
    _apply_matpower_start_modes!(net, mpc, pf_cfg.start_mode, mat_cfg; performance_profile = performance_profile)
    if pf_cfg.start_mode.voltage_mode in (:all_bus_vm, :profile_blend) || pf_cfg.start_mode.angle_mode in (:bus_va_blend, :matpower_va)
      net.flatstart = false
    end
  end
  return net
end
