[Mesh]
  type = GeneratedMesh
  dim = 2
  nx =
  ny =
  nz = 0
  xmin =
  xmax =
  ymin =
  ymax =
  elem_type = QUAD4
[]

[Variables]
  active = 'phi psi c'

  [./phi]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = phi_txt
    [../]
  [../]
  [./psi]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = psi_txt
    [../]
  [../]
  [./c]
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

[Kernels]
  # Order parameter phi
  [./dphidt]
    type = TimeDerivative
    variable = phi
  [../]
  [./ACBulk_phi]
    type = AllenCahn
    variable = phi
    coupled_variables = 'c'
    mob_name = L_phi
    f_name = g_phi
  [../]
  [./ACInterface_phi]
    type = ACInterface
    variable = phi
    mob_name = L_phi
    kappa_name = kappa_phi
  [../]
  [./noise_conserved_langevin]
    type = ConservedLangevinNoise
    amplitude =
    variable = phi
    noise = uniform_noise
  []

  # Order parameter psi
  [./dpsidt]
    type = TimeDerivative
    variable = psi
  [../]
  [./ACBulk_psi]
    type = AllenCahn
    variable = psi
    coupled_variables = 'c'
    mob_name = L_psi
    f_name = g_psi
  [../]
  [./ACInterface_psi]
    type = ACInterface
    variable = psi
    mob_name = L_psi
    kappa_name = kappa_psi
  [../]

  # Order parameter c
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./phi_c]
    type = CoefCoupledTimeDerivative
    v = 'phi'
    variable = c
    coef =
  [../]
  [./psi_c]
    type = CoefCoupledTimeDerivative
    v = 'psi'
    variable = c
    coef =
  [../]
  [./c_diffusion]
    type = ACInterface
    kappa_name = kappa_c
    mob_name = L_c
    variable = c
  [../]
[]

[BCs]
  [./Periodic]
    [./per_phi_xy]
      variable = phi
      auto_direction = 'x y'
    [../]
    [./per_psi_xy]
      variable = psi
      auto_direction = 'x y'
    [../]
    [./per_c_xy]
      variable = c
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./consts]
    # L_psi or L_phi can be changed to play on the influence of the dissolution/precipitation kinetics
    type = GenericConstantMaterial
    prop_names  = 'L_psi kappa_psi L_phi kappa_phi L_c'
    prop_values =
  [../]

  [./free_energy_phi]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_phi
    coupled_variables = 'phi psi c'
    material_property_names = 'g_phi_a(phi,psi,c) g_phi_b(phi,psi,c)'
    expression = 'if(psi<0.5, g_phi_a, g_phi_b)'
    enable_jit = true
    derivative_order = 1
    #outputs = exodus
  [../]
  [./free_energy_phi_a]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_phi_a
    coupled_variables = 'phi psi c'
    constant_names = 'W x_c c_eq'
    constant_expressions =
    expression = 'W*(phi^2)*((1-phi)^2) + x_c*(c-c_eq)*(1-(3*phi^2-2*phi^3))'
    enable_jit = true
    derivative_order = 1
    #outputs = exodus
  [../]
  [./free_energy_phi_b]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_phi_b
    coupled_variables = 'phi psi c'
    constant_names = 'W'
    constant_expressions =
    expression = 'W*(phi^2)*((1-phi)^2) - W*0.1*(1-(3*phi^2-2*phi^3))'
    enable_jit = true
    derivative_order = 1
    #outputs = exodus
  [../]

  [./free_energy_psi]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_psi
    coupled_variables = 'phi psi c'
    constant_names = 'W x_c c_eq'
    constant_expressions =
    expression = 'W*(psi^2)*((1-psi)^2) - x_c*(c-c_eq)*(3*psi^2-2*psi^3)'
    enable_jit = true
    derivative_order = 2
    #outputs = exodus
  [../]

  [./var]
    type = ParsedMaterial
    property_name = kappa_c
    coupled_variables = 'psi phi'
    constant_names = 'k_c_0 k_c_exp'
    constant_expressions =
    expression = 'k_c_0*(1-psi)*exp(-k_c_exp*phi)'
  [../]

  [./pore]
    type = ParsedMaterial
    property_name = pore_mat
    coupled_variables = 'psi phi'
    expression = '(1-psi)*(1-phi)'
  [../]
[]

[Functions]
  [phi_txt]
    type = PiecewiseMultilinear
    data_file = txt/phi.txt
  []
  [psi_txt]
    type = PiecewiseMultilinear
    data_file = txt/psi.txt
  []
  [c_txt]
    type = PiecewiseMultilinear
    data_file = txt/c.txt
  []
[]

[Preconditioning]
  # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
  # kernels compute all Jacobian matrix entries.
  # This allows us to use the Newton solver below.
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  # Automatic differentiation provides a _full_ Jacobian in this example
  # so we can safely use NEWTON for a fast solve
  solve_type = 'NEWTON'

  l_max_its = 20
  l_tol =
  l_abs_tol =

  nl_max_its = 10
  nl_rel_tol =
  nl_abs_tol =

  start_time = 0.0
  num_steps =

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt =
  [../]
[]

[UserObjects]
  [./uniform_noise]
    type = ConservedMaskedUniformNoise
    mask = pore_mat
  [../]
[]

[Postprocessors]
  [phi_pp]
    type = ElementAverageValue
    variable = phi
  []
  [psi_pp]
    type = ElementAverageValue
    variable = psi
  []
  [c_pp]
    type = ElementAverageValue
    variable = c
  []
  [sum_mat_pp]
    type = LinearCombinationPostprocessor
    pp_names = 'phi_pp psi_pp c_pp'
    pp_coefs =
  []
[]

[UserObjects]
  [john_connor]
    type = Terminator
    expression = 'phi_pp > 0.9'
    fail_mode = HARD
    execute_on = TIMESTEP_END
  []
  [sarah_connor]
    type = Terminator
    expression = 'phi_pp < 1e-5'
    fail_mode = HARD
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  [./other]
    type = VTK
    execute_on = 'TIMESTEP_END'
  [../]
  [console]
    type = Console
    execute_on = 'nonlinear'
    all_variable_norms = true
  []
  [./csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
    show = 'phi_pp psi_pp c_pp sum_mat_pp'
  [../]
[]
