"""Periodic ecology with the UNCHANGED v4 expressed/silent capacity process.

D expresses a focal structure; C reproduces without expressing it but retains
its functional toolkit. Losing that toolkit exits this genealogy, not the
organism's population. Deterministic environmental phases are shared by all
lineages. This is neither an empirically fitted flower model nor a random-
environment process. Matrix products/branching criteria are established theory.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F
import math

from anto_slk_recurrent_capacity_v4 import CapacityProcess

IDENTITY=((F(1),F(0)),(F(0),F(1)))


def matmul(a,b):
    """Exact parent-row/child-column two-type mean multiplication."""
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def _radius(matrix):
    (a,b),(c,d)=(tuple(map(float,row)) for row in matrix)
    return (a+d+math.hypot(a-d,2*math.sqrt(b*c)))/2


@dataclass(frozen=True)
class PeriodicCapacity:
    """Known finite phase list, repeated indefinitely from phase zero.

    Every phase must have strictly positive mean entries. The embedding at
    full periods is then a nonsingular, irreducible multitype branching
    process. Restricting this API avoids silently applying its survival
    criterion to reducible or deterministic special cases.
    """
    phases: tuple

    def __post_init__(self):
        phases=tuple(self.phases)
        if not 1 <= len(phases) <= 64:
            raise ValueError('Supply 1..64 periodic phases')
        if any(not isinstance(x,CapacityProcess) for x in phases):
            raise TypeError('Each phase must be the v4 CapacityProcess')
        if any(not 0 < x.activation < 1 or not 0 < x.silencing < 1 for x in phases):
            raise ValueError('Strictly positive two-type mean entries required')
        object.__setattr__(self,'phases',phases)

    def mean_cycle(self):
        out=IDENTITY
        for phase in self.phases:
            out=matmul(out,phase.mean_matrix())
        return out

    def average_mean(self):
        """Arithmetic comparison only; NOT the periodic growth operator."""
        return tuple(tuple(sum(p.mean_matrix()[i][j] for p in self.phases)/len(self.phases)
                           for j in range(2)) for i in range(2))

    def regime(self):
        (a,b),(c,d)=self.mean_cycle()
        if a >= 1 or d >= 1:
            return 'supercritical'
        delta=(1-a)*(1-d)-b*c
        return 'supercritical' if delta < 0 else 'critical' if delta == 0 else 'subcritical'

    def growth_per_generation(self):
        """Cycle Perron root normalized by period; not a fixation probability."""
        matrix=self.mean_cycle()
        scale=max(x for row in matrix for x in row)
        normalized=tuple(tuple(x/scale for x in row) for row in matrix)
        log_radius=(math.log(scale.numerator)-math.log(scale.denominator)
                    +math.log(_radius(normalized)))
        return math.exp(log_radius/len(self.phases))

    def common_loss_boundary(self):
        """Signed d*=1-rho(P_0)^(-1/L), only for common per-child loss d."""
        if len({p.loss for p in self.phases}) != 1:
            raise ValueError('The boundary assumes the SAME loss in all phases')
        return 1-float(1-self.phases[0].loss)/self.growth_per_generation()

    def exact_pgf(self,z):
        """Chronological cycle PGF F0(F1(...F_last(z))); short rational oracle.

        Truncate no offspring states. Restrict exact Fraction length because
        denominators double their number of digits each generation.
        """
        if len(self.phases)>10:
            raise ValueError('Exact rational PGF limited to <=10 phases')
        if len(z)!=2 or any(not 0<=v<=1 for v in z):
            raise ValueError('Supply two probabilities')
        for phase in reversed(self.phases):
            z=phase.pgf(z)
        return z

    def survival(self,tolerance=1e-12,max_iterations=100000):
        """Ultimate survival of COMPETENT genealogy from D or C at phase zero.

        Reuse v4's offspring PGFs but compose in reverse computational order.
        Extinction is the smallest fixed point. Lower iteration starts at 0;
        a componentwise supersolution below 1 is found along a Perron vector
        by backtracking. Both bounds are 70-digit Decimal approximations, NOT
        outward-rounded interval certificates. Return unresolved if the
        bracket has not closed; residual alone is not a convergence proof.
        """
        if isinstance(tolerance,bool) or not math.isfinite(tolerance) or not 0<tolerance<1:
            raise ValueError('Tolerance must be finite and in (0,1)')
        if type(max_iterations) is not int or max_iterations<1:
            raise ValueError('max_iterations must be a positive integer')
        regime=self.regime()
        common={'regime':regime,'period':len(self.phases),
                'growth_per_generation':self.growth_per_generation(),
                'decimal_digits':70,'numerical_brackets_not_interval_certificates':True}
        if regime!='supercritical':
            return dict(common,converged=True,survival_lower=(0.,0.),
                        survival_upper=(0.,0.),width=0.,residual=0.,iterations=0)
        with localcontext() as ctx:
            ctx.prec=70
            def dec(x): return Decimal(x.numerator)/Decimal(x.denominator)
            (a,b),(c,d)=(tuple(dec(v) for v in row) for row in self.mean_cycle())
            radius=(a+d+((a-d)**2+4*b*c).sqrt())/2
            if radius<=1:
                return dict(common,converged=False,survival_lower=(0.,0.),
                            survival_upper=(1.,1.),width=1.,residual=None,iterations=0)
            phase_rows=[tuple(tuple(dec(x)/2 for x in row) for row in p.mean_matrix())
                        for p in self.phases]
            def compose(z):
                for rows in reversed(phase_rows):
                    z=tuple((1-row[0]-row[1]+row[0]*z[0]+row[1]*z[1])**2
                            for row in rows)
                return z
            v=(b,radius-a)
            scale=max(v)
            v=tuple(x/scale for x in v)
            eps=Decimal('.5')
            for _ in range(200):
                upper=tuple(1-eps*x for x in v)
                trial=compose(upper)
                if all(trial[i]<=upper[i] for i in (0,1)) and all(x<1 for x in upper):
                    break
                eps/=2
            else:
                return dict(common,converged=False,survival_lower=(0.,0.),
                            survival_upper=(1.,1.),width=1.,residual=None,iterations=0)
            lower=(Decimal(0),Decimal(0))
            tol=Decimal(str(tolerance))/2
            for iteration in range(1,max_iterations+1):
                lower,upper=compose(lower),compose(upper)
                width=max(upper[i]-lower[i] for i in (0,1))
                if min(upper[i]-lower[i] for i in (0,1)) < -Decimal('1e-60'):
                    raise ArithmeticError('Extinction numerical brackets crossed')
                if width<=tol: break
            slo=tuple(float(max(0,min(1,1-x))) for x in upper)
            shi=tuple(max(slo[i],float(min(1,1-lower[i]))) for i in (0,1))
            middle=tuple((lower[i]+upper[i])/2 for i in (0,1))
            image=compose(middle)
            return dict(common,converged=width<=tol,survival_lower=slo,
                        survival_upper=shi,width=max(shi[i]-slo[i] for i in (0,1)),
                        residual=float(max(abs(image[i]-middle[i]) for i in (0,1))),
                        iterations=iteration)


def symmetric_period(run_length,switching='1/10',loss='1/5'):
    """A^k B^k: same environment proportions, opposite expression advantages."""
    if type(run_length) is not int or not 1<=run_length<=32:
        raise ValueError('Run length must be an integer in 1..32')
    a=CapacityProcess('9/10','2/5',loss,switching,switching)
    b=CapacityProcess('2/5','9/10',loss,switching,switching)
    return PeriodicCapacity((a,)*run_length+(b,)*run_length)


def panel():
    rows=[]
    for switch in ('1/10','1/2','9/10'):
        for k in (1,2,3,4,8):
            p=symmetric_period(k,switch)
            result=p.survival()
            if not result['converged']:
                raise ArithmeticError('Refuse to publish an unconverged survival value')
            rows.append({'phenotype_switching':switch,'environment_run_length':k,
                         'common_loss_boundary':p.common_loss_boundary(),
                         'first_D_child_from_C_at_phase0':float(p.phases[0].first_expression_from_silent()),
                         **result})
    return {'schema':'anto_slk_periodic_capacity_v5','date':'2026-10-20',
            'environment_process':'deterministic periodic, shared across descendants',
            'natural_regain_demonstrated':False,'organism_extinction_inferred':False,
            'fixation_or_optimal_strategy_inferred':False,'slk_payoff_calibrated':False,
            'empirical_owner':'zuizui0223/aza3',
            'q_A':{'D':'9/10','C':'2/5'},'q_B':{'D':'2/5','C':'9/10'},
            'common_loss':'1/5','environment_proportions':[.5,.5],
            'arithmetic_mean_growth':1.04,'pure_D_or_pure_C_geometric_growth':.96,
            'starting_unit':'one silent competent C at first A generation',
            'scenarios':rows}
