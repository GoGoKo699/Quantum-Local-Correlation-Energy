#!/usr/bin/env python3
"""Validate substantive numerical identities used in the resource note.

This validates shell normalizations, fixed-point cosine arithmetic against
rational enclosures, the comparator amplitude, and the aggregate budget.
It does not compile or simulate the complete quantum circuit.
"""
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from resource_ledger import constants, instance, taylor_plan


def atan_interval(q, bits):
    value=Fraction(0)
    j=0
    while True:
        term=Fraction((-1)**j,(2*j+1)*q**(2*j+1))
        value+=term
        j+=1
        next_term=Fraction((-1)**j,(2*j+1)*q**(2*j+1))
        if abs(next_term)<Fraction(1,1<<bits):
            return min(value,value+next_term),max(value,value+next_term)


def pi_interval(bits):
    a,b=atan_interval(5,bits+8)
    c,d=atan_interval(239,bits+8)
    return 16*a-4*d,16*b-4*c


def fixed_cos(m,L,d,f,pi_mid):
    a=abs(m)
    sign=1 if 4*a<=L else -1
    r=a if sign==1 else L//2-a
    scale=1<<f
    x=int(2*pi_mid*r*scale/L)
    t=x*x//scale
    coefficients=[((-1)**j*scale)//math.factorial(2*j) for j in range(d+1)]
    y=coefficients[-1]
    for j in range(d-1,-1,-1):
        y=y*t//scale+coefficients[j]
    return Fraction(sign*y,scale),r,sign


def cosine_validation(L,b):
    n=L.bit_length()-1
    plan=taylor_plan(n,b)
    d,f,s=plan['cosine_degree_in_x_squared'],plan['fraction_bits'],plan['cosine_error_bits']
    lo,hi=pi_interval(f+30)
    pi_mid=(lo+hi)/2
    worst=Fraction(0)
    modes=range(-L//2,L//2) if L<=32 else sorted(set([-L//2,-L//2+1,-L//4,-1,0,1,L//4,L//2-1]))
    for m in modes:
        computed,r,sign=fixed_cos(m,L,d,f,pi_mid)
        x=2*pi_mid*r/L
        reference=sum((Fraction((-1)**j,math.factorial(2*j))*x**(2*j) for j in range(d+7)),Fraction(0))
        # Alternating Taylor enclosure on 0<=x<=pi/2, plus input-pi interval.
        remainder=x**(2*(d+7))/math.factorial(2*(d+7))
        error=abs(computed-sign*reference)+remainder+(hi-lo)
        worst=max(worst,error)
    assert worst<Fraction(1,1<<s)
    return dict(L=L,amplitude_bits=b,tested_modes=len(list(modes)),
                cosine_absolute_error_certified_on_tested_modes=float(worst),
                required_uniform_cosine_error=2.**(-s),
                scope="Rational enclosures on tested modes; general uniform proof is in the note")


def shell_validation(L):
    total_ref=total_green=0.
    worst=0.
    largest_correction=0.
    for m in itertools.product(range(-L//2,L//2),repeat=3):
        if m==(0,0,0):
            continue
        s=0
        while not all(-2**s<=u<2**s for u in m):
            s+=1
        lam=4*sum(math.sin(math.pi*u/L)**2 for u in m)
        reference=1/(8*L*4**s)
        correction=4**(s+1)/(L*L*lam)
        target=1/(2*L**3*lam)
        worst=max(worst,abs(reference*correction/target-1))
        largest_correction=max(largest_correction,correction)
        total_ref+=reference
        total_green+=reference*correction
    assert worst<1e-13 and largest_correction<=1+1e-13
    assert abs(total_ref-7/8*(1-1/L))<1e-12
    return dict(L=L,reference_success_norm_squared=total_ref,
                maximum_shell_correction_squared=largest_correction,
                maximum_relative_amplitude_squared_identity_error=worst,
                Green_half_success_norm_squared=total_green)


def comparator_validation():
    b=8
    scale=1<<b
    worst=0.
    for p in [Fraction(0),Fraction(1,4),Fraction(3,8),Fraction(1,2),Fraction(3,4),Fraction(1)]:
        radicand=4*p*(1-p)
        accepted=sum(Fraction(j*j,scale*scale)<radicand for j in range(scale))
        amplitude=accepted/scale
        error=abs(amplitude-math.sqrt(float(radicand)))
        assert error<=1/scale+1e-15
        worst=max(worst,error)
    return dict(bits=b,maximum_amplitude_error=worst,bound=1/scale,
                convention="strict squared predicate yields ceil, with zero endpoint handled exactly")


def main():
    eps=constants()['default_epsilon']
    budget_checks=[]
    for epsilon in [eps,1e-8,1e-10]:
        r=instance(8,epsilon,.01)
        assert r['total_scalar_error_budget']<epsilon
        assert r['amplitude_estimation']['median_failure_upper']<.01
        table=r['geometric_weight_table']
        assert table['observed_L1_weight_error']<table['analytic_L1_weight_error_upper']
        assert table['integer_endpoints']==sorted(table['integer_endpoints'])
        budget_checks.append(dict(epsilon=epsilon,total_error_bound=r['total_scalar_error_budget'],
            failure_bound=r['amplitude_estimation']['median_failure_upper'],
            degree=r['precision']['neumann_degree']))
    out=dict(status="passed_mathematical_primitive_and_budget_checks_not_full_circuit_simulation",
        shell_checks=[shell_validation(L) for L in [4,8,16]],
        fixed_point_cosine_checks=[cosine_validation(L,b) for L,b in [(8,8),(16,26),(256,26)]],
        comparator_check=comparator_validation(), budget_checks=budget_checks)
    path=(Path(__file__).parent / 'results' / 'resource_validation.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
