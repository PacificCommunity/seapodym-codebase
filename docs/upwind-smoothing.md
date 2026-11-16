# Optional C1 smoothing of upwind velocity splitting

The coupled model accepts an optional parameter:

```xml
<upwind_smoothing_epsilon value="0"/>
```

The default is **zero**, which retains the original transport stencil and its
adjoint branches. No change to an existing parameter file is required.
For a positive epsilon, the direct and manually coded adjoint use the same split.
The epsilon is a fixed numerical setting, not an estimated biological parameter.

## Configuration

Configuration is read only from the XML parameter file. If epsilon is not set,
the smoothing is not done. An explicit XML value of zero disables smoothing.
The startup log records epsilon and whether it came from the parameter file
or the default. The value must be finite and nonnegative; invalid values stop
initialization. Include the XML setting explicitly when recording a run.

The epsilon has the units of the velocity passed to the transport stencil.
There is no universal recommended positive value. Assess epsilon against
the chosen grid, timestep, forcings and scientific application.

## Definition

For velocity v and epsilon > 0, replace its absolute value by

```text
phi(v, epsilon) = |v|                         if |v| >= epsilon
                  (v*v + epsilon*epsilon)
                  -----------------------    if |v| < epsilon.
                         2*epsilon

v_positive = (v + phi)/2
v_negative = (v - phi)/2
```

The smoothing band is |v| < epsilon, where phi has its minimum epsilon/2 at
v = 0. The added numerical diffusion is at most epsilon*dx/4, at v = 0.

Inside the band, the derivatives with respect to v are
`(1 + v/epsilon)/2` and `(1 - v/epsilon)/2`, respectively. The value and first
derivative agree at the joins. This construction is C1; its second derivative
jumps at the joins and it does not make the complete objective globally C2.

The original arithmetic and branch choices are kept when smoothing is disabled
or the velocity is outside the band. The small positive/negative splitting
inside the band introduces additional numerical diffusion. It preserves the
conservative stencil structure but can alter the simulated solution.

## Build and use

Use the usual model build procedure after updating the sources. Existing
parameter files remain usable. To make a run reproducible, include the XML
setting explicitly in its parameter file.

Changing epsilon can change the simulated solution through additional numerical
diffusion. Compare identical parameter vectors before attributing a difference
to reoptimization.
