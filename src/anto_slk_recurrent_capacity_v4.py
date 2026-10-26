"""Two-type branching of expressed and silent-but-competent descendants.

Theory-only, constant environment, independent reproductive slots and lineages.
D expresses the focal candidate; C is reproductively active but does NOT express
it. C is not a dormant seed. Offspring losing the necessary toolkit exit this
process; their organismal descendants may still survive outside it.

Unlike v3, expression is not perfectly inherited, and capacity can be lost
AFTER first expression. No finite census, arbitrary abundance threshold, priming
state, SLK-sign eligibility mask, or additional positive field result is used.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F
import math


def probability(value):
    if isinstance(value, bool):
        raise ValueError('Boolean is not a probability')
    try:
        x = F(str(value))
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ValueError('Expected a finite rational probability') from exc
    if not 0 <= x <= 1:
        raise ValueError('Probability outside [0,1]')
    return x


@dataclass(frozen=True)
class CapacityProcess:
    """Two offspring slots; reproduction, retention, then expression assignment.

    Each slot of parent type i makes a viable child with probability q_i. That
    child keeps the required toolkit with probability 1-loss. Conditional on
    retention, D-parent children become silent C with probability silencing;
    C-parent children express D with probability activation. All other slot
    outcomes are counted as no COMPETENT descendant, not necessarily no child.
    Rates are assigned, not inferred from a SLK payoff or pigment observation.

    The declared domain is irreducible and nonsingular: 0<q_i<1, loss<1,
    0<silencing,activation<=1. Reducible limits need separate type-reachability
    analysis and are intentionally not accepted by this survival API.
    """
    q_expressed: F
    q_silent: F
    loss: F
    silencing: F
    activation: F

    def __post_init__(self):
        for name in ('q_expressed','q_silent','loss','silencing','activation'):
            object.__setattr__(self, name, probability(getattr(self, name)))
        if not 0 < self.q_expressed < 1 or not 0 < self.q_silent < 1:
            raise ValueError('Require 0<q_D,q_C<1 for the nonsingular model')
        if self.loss == 1 or self.silencing == 0 or self.activation == 0:
            raise ValueError('This API requires an irreducible two-type process')

    def mean_matrix(self):
        d = 2*self.q_expressed*(1-self.loss)
        c = 2*self.q_silent*(1-self.loss)
        return ((d*(1-self.silencing), d*self.silencing),
                (c*self.activation, c*(1-self.activation)))

    def pgf(self, z):
        """Offspring PGFs F_i(z) = (1 - row_sum(M)/2 + M_i z/2)^2."""
        if len(z) != 2 or any(not 0 <= x <= 1 for x in z):
            raise ValueError('Supply two probability arguments')
        return tuple((1-(row[0]+row[1])/2+(row[0]*z[0]+row[1]*z[1])/2)**2
                     for row in self.mean_matrix())

    def spectral_radius(self):
        (a,b),(c,d) = (tuple(map(float,row)) for row in self.mean_matrix())
        return (a+d+math.hypot(a-d,2*math.sqrt(b*c)))/2

    def regime(self):
        """Exact rational classification, including exactly critical cases."""
        (a,b),(c,d) = self.mean_matrix()
        if a >= 1 or d >= 1:
            return 'supercritical'  # both cross-terms are strictly positive
        determinant = (1-a)*(1-d)-b*c
        return ('supercritical' if determinant < 0 else
                'critical' if determinant == 0 else 'subcritical')

    def first_expression_from_silent(self):
        """P(at least one D child next generation | one competent C parent)."""
        return 1-(1-self.q_silent*(1-self.loss)*self.activation)**2

    def loss_boundary(self):
        """Signed d*=1-1/rho(M at zero loss); survival requires loss<d*."""
        return 1-(1-float(self.loss))/self.spectral_radius()


def finite_extinction(p, generations):
    """Exact probability of no competent descendant at generation n.

    These approach ultimate extinction from below; a positive finite survival
    probability must never be interpreted as positive permanent survival.
    Bound n to avoid exponentially large exact-rational denominators.
    """
    if type(generations) is not int or not 0 <= generations <= 10:
        raise ValueError('Exact-rational finite oracle requires integer 0<=n<=10')
    e = (F(0),F(0))
    for _ in range(generations):
        e = p.pgf(e)
    return e


def survival(p, tolerance=1e-12, max_iterations=100000):
    """Infinite competent-lineage survival via two-sided fixed-point iteration.

    Extinction e is the smallest PGF fixed point. Start its lower sequence at
    zero. For Perron vector v normalized to max(v)=1 and eigenvalue r>1,
    u=1-eps*v is a supersolution when eps<=4*(r-1)/r^2. Both sequences converge
    to e in this nonsingular irreducible model. Computation uses 50-digit
    Decimal arithmetic to avoid near-critical binary64 fixed-point plateaus.
    Returned brackets are numerical, not outward-rounded interval certificates.
    Exact critical/subcritical classification takes precedence over iteration.
    """
    if isinstance(tolerance,bool) or not math.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError('Tolerance must be finite and between zero and one')
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')
    regime = p.regime()
    common = {'regime':regime, 'spectral_radius':p.spectral_radius(),
              'decimal_digits':50, 'numerical_brackets_not_interval_certificates':True}
    if regime != 'supercritical':
        return dict(common, converged=True, survival_lower=(0.0,0.0),
                    survival_upper=(0.0,0.0), width=0.0, residual=0.0, iterations=0)
    with localcontext() as ctx:
        ctx.prec=50
        def dec(x):
            return Decimal(x.numerator)/Decimal(x.denominator)
        (a,b),(c,d) = (tuple(dec(x) for x in row) for row in p.mean_matrix())
        r = (a+d+((a-d)**2+4*b*c).sqrt())/2
        if r <= 1:  # exact positive gap smaller than current Decimal resolution
            return dict(common, converged=False, survival_lower=(0.0,0.0),
                        survival_upper=(1.0,1.0), width=1.0, residual=None, iterations=0)
        v = (b, r-a)
        scale = max(v)
        v = tuple(x/scale for x in v)
        eps = min(Decimal('0.5'), 2*(r-1)/(r*r))
        lower = (Decimal(0),Decimal(0))
        upper = tuple(1-eps*x for x in v)
        tol=Decimal(str(tolerance))/2  # reserve room for binary64 serialization

        def g(z):
            return ((1-(a+b)/2+(a*z[0]+b*z[1])/2)**2,
                    (1-(c+d)/2+(c*z[0]+d*z[1])/2)**2)

        for iteration in range(1,max_iterations+1):
            lower,upper = g(lower),g(upper)
            width = max(upper[i]-lower[i] for i in (0,1))
            if width < -Decimal('1e-45'):
                raise ArithmeticError('Numerical fixed-point brackets crossed')
            if width <= tol:
                break
        slo = tuple(float(max(0,min(1,1-x))) for x in upper)
        shi = tuple(max(slo[i],float(min(1,1-lower[i]))) for i in (0,1))
        middle = tuple((lower[i]+upper[i])/2 for i in (0,1))
        check = g(middle)
        return dict(common, converged=width<=tol, survival_lower=slo,
                    survival_upper=shi, width=max(shi[i]-slo[i] for i in (0,1)),
                    residual=float(max(abs(check[i]-middle[i]) for i in (0,1))), iterations=iteration)


def activation_boundary(p):
    """Exact interior threshold when retained C output exceeds 1 and D does not.

    Write v_D=2q_D(1-loss), v_C=2q_C(1-loss). For v_D<1<v_C,
    a*=[1-v_D(1-s)]*(v_C-1)/[v_C*(1-v_D)]. If 0<a*<1,
    irreducible capacity survives with positive probability iff activation<a*.
    This is not a universal direction of the reactivation effect.
    """
    vd,vc = 2*p.q_expressed*(1-p.loss), 2*p.q_silent*(1-p.loss)
    if not vd < 1 < vc:
        raise ValueError('This boundary requires v_D<1<v_C')
    result = (1-vd*(1-p.silencing))*(vc-1)/(vc*(1-vd))
    if not 0 < result < 1:
        raise ValueError('No interior activation crossing in this parameter slice')
    return result


def initial_capacity_survival(p, theta, expressed_fraction=0):
    """Initial resource frequency changes chance of survival, not criticality."""
    theta,expressed_fraction = probability(theta),probability(expressed_fraction)
    result=survival(p)
    def weighted(values):
        return float(theta)*(float(expressed_fraction)*values[0]+(1-float(expressed_fraction))*values[1])
    return {'lower':weighted(result['survival_lower']), 'upper':weighted(result['survival_upper']),
            'regime':result['regime'], 'converged':result['converged']}


def panel():
    rows=[]
    for i in range(1,20):
        p=CapacityProcess('13/20','9/10','3/10','1/4',F(i,20))
        result=survival(p)
        if not result['converged']:
            raise ArithmeticError('Refuse a scenario table with unconverged survival')
        rows.append({'activation':float(p.activation),
                     'first_expression_probability':float(p.first_expression_from_silent()),
                     'mean_matrix':[[float(x) for x in row] for row in p.mean_matrix()],
                     **result})
    p=CapacityProcess('13/20','9/10','3/10','1/4','3/10')
    return {'schema':'anto_slk_recurrent_capacity_v4', 'date':'2026-10-20',
            'process':'constant-environment irreducible two-type Galton-Watson',
            'natural_regain_demonstrated':False, 'organism_extinction_inferred':False,
            'fixation_inferred':False, 'slk_payoff_calibrated':False,
            'silent_state_reproduces_and_is_not_seed_dormancy':True,
            'empirical_owner':'zuizui0223/aza3',
            'parameters':{'q_expressed':'13/20','q_silent':'9/10','loss':'3/10','silencing':'1/4'},
            'activation_critical_exact':str(activation_boundary(p)),
            'activation_critical':float(activation_boundary(p)),
            'loss_critical_at_activation_03':p.loss_boundary(),
            'activation_sweep':rows}
