# Optional stock likelihood weight

The existing stock constraint accepts a fixed, nonnegative weight for each
configured species. Add the optional `weight` attribute to the existing species
element under `/mean_stock_obs`. Keep its `value` and geographical attributes.

Illustrative XML fragment (replace the placeholder and retain the existing
geographical attributes of your parameter file):

```xml
<mean_stock_obs>
    <alb value="YOUR_EXISTING_REFERENCE" weight="1" />
</mean_stock_obs>
```

The default is **1** when the attribute is absent or empty, preserving the unweighted
calculation. An explicit weight must be a finite real number greater than or
equal to zero. Invalid values stop initialization. The attribute is read when
the species stock likelihood is enabled; explicit weights appear in the startup
log. It is a configuration setting, not an estimated biological parameter.

## Effect on the objective and derivatives

The existing one-sided penalty is retained:

```text
J_stock = 0.5 * weight * max(B - B_reference, 0)^2
```

Increasing the weight strengthens this contribution relative to the other
objective components. A zero weight removes its finite contribution. The stock
reference, geographical domain and activation flag retain their existing roles.
Changing the weight does not change the reference or introduce a new penalty.

The multiplication remains inside an ADMB differentiable expression. Its
gradient contribution is therefore scaled by the same weight; no separate
manual adjoint formula is required. The existing transition at the stock
reference is retained, including its second-derivative limitation.

Weights change the objective being optimized. Choose and report them for the
application; this feature does not prescribe any particular weight above one.
