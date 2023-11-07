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
  active = 'psi phi c'

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
  active = 'dpsidt ACBulk_psi ACInterface_psi psi_phi dphidt ACBulk_phi ACInterface_phi dcdt psi_c phi_c c_diffusion'

  #
  # Order parameter phi
  #
  [./dphidt]
    type = TimeDerivative
    variable = phi
  [../]
  [./ACBulk_phi]
    type = AllenCahn
    variable = phi
    args = 'c'
    mob_name = L_phi
    f_name = g_phi
  [../]
  [./ACInterface_phi]
    type = ACInterface
    variable = phi
    mob_name = L_phi
    kappa_name = kappa_phi
  [../]
  [./psi_phi]
    type = CoefCoupledTimeDerivative
    v = 'psi'
    variable = phi
    coef = -1
  [../]

  #
  # Order parameter psi
  #
  [./dpsidt]
    type = TimeDerivative
    variable = psi
  [../]
  [./ACBulk_psi]
    type = AllenCahn
    variable = psi
    args = 'c'
    mob_name = L_psi
    f_name = g_psi
  [../]
  [./ACInterface_psi]
    type = ACInterface
    variable = psi
    mob_name = L_psi
    kappa_name = kappa_psi
  [../]

  #
  # Order parameter c
  #
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./phi_c]
    type = CoefCoupledTimeDerivative
    v = 'phi'
    variable = c
    coef = 1
  [../]
  [./psi_c]
    type = CoefCoupledTimeDerivative
    v = 'psi'
    variable = c
    coef = 2.35
  [../]
  [./c_diffusion]
    type = ACInterface
    kappa_name = kappa_c
    mob_name = L_c
    variable = c
  [../]
[]

[Materials]
  active = 'consts free_energy_psi free_energy_phi'

  [./consts]
    # L_psi or L_phi can be changed to play on the influence of the dissolution/precipitation kinetics
    type = GenericConstantMaterial
    prop_names  = 'L_psi kappa_psi L_phi kappa_phi L_c kappa_c'
    prop_values = '1 0.03 1 0.03 1 0.05'
  [../]

  [./free_energy_phi]
    type = DerivativeParsedMaterial
    block = 0
    f_name = g_phi
    args = 'phi c'
    constant_names = 'W x_c'
    constant_expressions =
    function = 'W*16*(phi^2)*((1-phi)^2) - x_c*(c-1)*(phi^2)*(3-2*phi)'
    enable_jit = true
    derivative_order = 2
  [../]

  [./free_energy_psi]
    type = DerivativeParsedMaterial
    block = 0
    f_name = g_psi
    args = 'psi c'
    constant_names = 'A_psi B_psi C_psi x_c'
    constant_expressions =
    function = 'A_psi*psi^4+B_psi*psi^3+C_psi*psi^2 - x_c*(c-1)*(psi^2)*(3-2*psi)'
    enable_jit = true
    derivative_order = 2
  [../]

  [./var]
    type = ParsedMaterial
    property_name = kappa_psi
    coupled_variables = 'psi'
    function = 'if(psi>0.5, 0.03, 0.01)'
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
  l_tol = 1.0e-4
  l_abs_tol = 1.0e-6

  nl_max_its = 10
  nl_rel_tol = 1.0e-4
  nl_abs_tol = 1.0e-6

  start_time = 0.0
  end_time   = 10

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.1
  [../]
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  [./other]
    type = VTK
  [../]
[]
