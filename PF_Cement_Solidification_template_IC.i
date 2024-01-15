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
  active = 'psi'

  [./psi]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function =
    [../]
  [../]
[]

[Kernels]
  active = 'dpsidt CHBulk_psi CHInterface_psi'

  [./dpsidt]
    type = TimeDerivative
    variable = psi
  [../]
  [./CHBulk_psi]
    type = CahnHilliard
    variable = psi
    mob_name = L_psi
    f_name = g_psi
  [../]
  [./CHInterface_psi]
    type = CHInterface
    variable = psi
    mob_name = L_psi
    kappa_name = kappa_psi
  [../]
[]

[Materials]
  active = 'consts free_energy_psi'

  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'L_psi kappa_psi'
    prop_values =
  [../]
  [./free_energy_psi]
    type = DerivativeParsedMaterial
    block = 0
    property_name = g_psi
    coupled_variables = 'psi'
    constant_names = 'W'
    constant_expressions =
    expression = 'W*16*(psi^2)*((1-psi)^2)'
    enable_jit = true
    derivative_order = 1
  [../]
[]

[Functions]
  active = 'psi_txt'

  [psi_txt]
    type = PiecewiseMultilinear
    data_file = txt/psi.txt
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
  l_tol = 1.0e-3
  l_abs_tol = 1.0e-5

  nl_max_its = 10
  nl_rel_tol = 1.0e-3
  nl_abs_tol = 1.0e-5

  start_time = 0.0
  num_steps = 100

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 10000
  [../]
[]

[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true
  [./other]
    type = VTK
  [../]
[]
