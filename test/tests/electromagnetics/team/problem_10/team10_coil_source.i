# TEAM Problem 10: the coil current
# =================================
#
# Sub-app of team10.i, which is the input to run. It is launched automatically,
# once, at the start of the run.
#
# The benchmark's coil is a winding of 162 turns carrying a uniform current
# density, and that field can be written down analytically. Projecting it
# straight onto the mesh is not good enough for the main solve: the curl-curl
# operator there has every gradient in its null space, so its right hand side has
# to be orthogonal to every discrete gradient. The projected analytic current is
# not, because the meshed coil boundary is faceted and a little current leaks
# across it.
#
# This app removes that leak. It solves, over the whole mesh,
#
#   (grad(phi), grad(psi)) = (J_applied, grad(psi))_coil    for all psi,
#
# so that the current the main app actually uses,
#
#   J = J_applied restricted to the coil, minus grad(phi),
#
# satisfies (J, grad(psi)) = 0 for every psi by construction: a Helmholtz
# projection, the smallest correction that makes the source discretely divergence
# free.
#
# The current is computed here at the peak ampere-turns. The main app scales both
# pieces by the ramp 1 - exp(-t / tau); the projection is linear, so the scaled
# current is still divergence free at every time.
#
# The two pieces are handed to the main app separately, as applied_current and
# correction_potential, because the first is integrated over the coil alone and the
# second over the whole domain.
#
# phi is pinned wherever the main app holds the tangential vector potential to
# zero: the far field and the x = 0 and y = 0 symmetry planes. That matches the
# gradients left in the main problem's null space, and lets the current cross
# those planes, as it has to: only the eighth x, y, z >= 0 is meshed. The z = 0
# plane carries the natural condition in both apps, so phi is left free there; the
# coil current has no z component and does not cross it.

ampere_turns = 913.68 # AT, 162 turns at 5.64 A; overridden by the main app
fe_order = FIRST # likewise; it has to match the main app's H(curl) space
mesh_file = team10.e # likewise

coil_core = 0.05 # m, half-width of the square the corner arcs are centred on
coil_area = '${fparse 0.025 * 0.100}' # m2, 25 mm wide by 100 mm tall
current_density = '${fparse ampere_turns / coil_area}' # A/m2

[Mesh]
  type = MFEMFileMesh
  file = ${mesh_file}
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = ${fe_order}
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = ${fe_order}
  []
[]

[Functions]
  # The specified winding current. d is the distance from the square core the corner
  # arcs are centred on, so grad(d) points radially outwards and has unit length, and
  # z_hat x grad(d) runs anticlockwise around the coil, seen from above, at uniform
  # magnitude. That drives flux up the centre plate. The floor on the denominator
  # only guards the core, where the numerators vanish too; outside the coil the
  # expression is not used.
  [applied_current_density]
    type = ParsedVectorFunction
    expression_x = '-J0 * (max(y-a,0)+min(y+a,0)) / max(sqrt((max(x-a,0)+min(x+a,0))^2 + (max(y-a,0)+min(y+a,0))^2), 1e-30)'
    expression_y = 'J0 * (max(x-a,0)+min(x+a,0)) / max(sqrt((max(x-a,0)+min(x+a,0))^2 + (max(y-a,0)+min(y+a,0))^2), 1e-30)'
    expression_z = '0'
    symbol_names = 'J0 a'
    symbol_values = '${current_density} ${coil_core}'
  []
[]

[Variables]
  [correction_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [applied_current]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  # The applied current is the right hand side of the solve, so it has to be in
  # place before the solve runs rather than at the end of the step.
  [project_applied_current]
    type = MFEMVectorProjectionAux
    variable = applied_current
    vector_coefficient = applied_current_density
    execute_on = TIMESTEP_BEGIN
  []
[]

[BCs]
  # symmetry_z is deliberately absent: phi is free on the plane where the main
  # app's condition is natural.
  [pin_tangential_a_boundary]
    type = MFEMScalarDirichletBC
    variable = correction_potential
    boundary = 'far_field symmetry_x symmetry_y'
    coefficient = 0.0
  []
[]

[Kernels]
  [correction_diffusion]
    type = MFEMDiffusionKernel
    variable = correction_potential
  []
  # (-J_applied, grad(psi)) over the coil, which is where the leak comes from.
  [applied_current_divergence]
    type = MFEMMixedVectorWeakDivergenceKernel
    variable = correction_potential
    trial_variable = applied_current
    block = 'Coil'
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [gmres]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 300
  []
[]

[Executioner]
  type = MFEMSteady
[]

# No Outputs block: both fields are transferred to the main app and written out
# from there.
