% Term-by-term comparison: dynamic 0-D model (instant-kinetics limit)
% vs the Simplified DAC energy model dac.energyModel.
raw = fileread('dac_consist.json');
J = jsondecode(raw);
P = J.P;  dyn = J.dyn;
C = dac.constants();

[E_total, E_thermal, E_elec, E_terms] = dac.energyModel(P, C);
% E_terms columns: sensible | CO2 des | H2O des | blower | vacuum  [MJ/kg]

fprintf('%-22s %14s %14s %12s\n', 'term', 'simplified', 'dynamic', 'rel.dev');
cmp('E_sensible', E_terms(1), dyn.E_sens);
cmp('E_CO2_des',  E_terms(2), dyn.E_CO2);
cmp('E_H2O_des',  E_terms(3), dyn.E_H2O);
cmp('E_thermal',  E_thermal,  dyn.E_th);
cmp('E_blower',   E_terms(4), dyn.E_fan);
cmp('E_vacuum',   E_terms(5), dyn.E_vac_heat);
cmp('E_total', E_total, dyn.E_th + dyn.E_fan + dyn.E_vac_heat);

function cmp(name, a, b)
  fprintf('%-22s %14.6f %14.6f %11.2e\n', name, a, b, abs(a-b)/max(abs(a),1e-12));
end
