#!/usr/bin/env python3
"""Conservative analytic logical-resource ledger, not a circuit compiler.

All floating-point output is an evaluation of the formulas in the companion
note. No quantum circuit, rotation synthesis, physical runtime, or speedup is
claimed. The only input field used by the priced example is the explicit
bit-reversal program shared with classical_benchmark.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class Gates:
    # Elementary reversible Boolean gates before external controls.
    boolean: int = 0
    hadamard: int = 0
    # Single or controlled exact dyadic phase gates, not synthesized T counts.
    phase: int = 0

    def __add__(self, other: "Gates") -> "Gates":
        return Gates(self.boolean + other.boolean,
                     self.hadamard + other.hadamard,
                     self.phase + other.phase)

    def __mul__(self, k: int) -> "Gates":
        return Gates(k*self.boolean, k*self.hadamard, k*self.phase)

    __rmul__ = __mul__


def add(w: int) -> Gates:
    # Cuccaro ripple adder, padded to cover endpoint/carry conventions.
    return Gates(boolean=9*w+8)


def compare(w: int) -> Gates:
    # Compute carry/borrow, retain result, reverse arithmetic.
    return 2*add(w) + Gates(boolean=4*w+8)


def multiply(w: int) -> Gates:
    # Signed operands extended to 2w bits, 2w controlled shifted additions.
    # Each controlled 2w-bit adder costs <= 28w+20 elementary Boolean gates.
    # Padding includes copies and sign extension; full 2w-bit product retained.
    return Gates(boolean=64*w*w+64*w)


def cdf_select(entries: int, bits: int, word: int) -> Gates:
    # Uniform register; cumulative-interval comparisons; SELECT; uncompute.
    return 2*entries*compare(max(bits, word)) + Gates(
        boolean=4*entries*max(bits, word), hadamard=2*bits)


def constants(g=.01, E0=2., to=.02, tv=.02, gamma=.5) -> dict:
    Z = 8*math.pi
    Gbar = 3*math.pi+math.pi**3/24
    alphaJ = 4*math.pi*g*gamma**2
    alphaK = 4*g*gamma*math.sqrt(2*math.pi*Z)
    physicalK = 4*g*gamma*math.sqrt(2*math.pi*Gbar)
    hopping = 6*(abs(to)+abs(tv))/E0
    beta_terms = [hopping, 3*g*Z/E0, 2*alphaJ/E0,
                  alphaK/(math.sqrt(2)*E0)]
    q_terms = [hopping, 3*g*Gbar/E0, 2*alphaJ/E0,
               physicalK/(math.sqrt(2)*E0)]
    B02 = alphaJ**2/2
    Bc2 = B02/(2*E0)
    Gamma = 6*g*Gbar+4*alphaJ+physicalK
    Delta = E0-6*(abs(to)+abs(tv))-Gamma
    old_alphaA = 2*E0+12*(abs(to)+abs(tv))+6*g*7*math.pi+4*alphaJ + 4*g*gamma*math.sqrt(14*math.pi**2)
    Lambda = 2*E0+12*(abs(to)+abs(tv))+Gamma
    return dict(g=g, E0=E0, t_o=to, t_v=tv, gamma_max=gamma, Z=Z,
                Gbar=Gbar, alpha_J=alphaJ, alpha_K=alphaK,
                physical_K_bound=physicalK, beta_terms=beta_terms,
                q_terms=q_terms, beta=sum(beta_terms), q=sum(q_terms),
                B0_squared=B02, Bc_squared=Bc2,
                original_Delta=Delta, original_alpha_A=old_alphaA,
                preconditioned_Delta=E0*(1-sum(q_terms)), Lambda=Lambda,
                default_epsilon=B02/(1600*Lambda*10))


def select_precision(c: dict, epsilon: float) -> dict:
    q, beta, Bc2 = c["q"], c["beta"], c["Bc_squared"]
    if not (0 <= q < 1 and 0 < beta < 1):
        raise ValueError("This bounded-geometric-LCU ledger requires q<1 AND beta<1.")
    p = 0
    while Bc2*q**(p+1)/(1-q) > epsilon/4:
        p += 1
    Lp = sum(beta**k for k in range(p+1))
    # 64 and 32 are deliberately padded normalized primitive composition bounds.
    b = 4
    while True:
        eta = 2.**(-b)
        dM, dB = 64*beta*eta, 32*eta
        r = q+dM
        error = math.inf if r >= 1 else Bc2*(dM/(1-r)**2+2*dB/(1-r))
        if error <= epsilon/8:
            break
        b += 1
    h = max(1, math.ceil(math.log2(16*Bc2*Lp*(p+1)/epsilon)))
    cdf_error = 2*Bc2*Lp*(p+1)*2.**(-h)
    probability_error = epsilon/(4*Bc2*Lp)
    M = 2
    while math.pi/M+math.pi**2/M**2 > probability_error:
        M *= 2
    return dict(neumann_degree=p, neumann_lcu_norm=Lp,
                tail_error=Bc2*q**(p+1)/(1-q), amplitude_bits=b,
                primitive_tolerance=eta, implemented_M_error=dM,
                implemented_source_error=dB, perturbed_norm_bound=r,
                primitive_scalar_error=error, geometric_cdf_bits=h,
                geometric_scalar_error=cdf_error,
                AE_probability_error=probability_error, AE_M=M,
                AE_error_bound=math.pi/M+math.pi**2/M**2,
                AE_energy_error_bound=2*Bc2*Lp*(math.pi/M+math.pi**2/M**2))


def confidence_repetitions(delta: float, qft_error: float = 1e-6) -> tuple[int, float]:
    # Approximate inverse QFT changes any outcome event by at most 2||dU||.
    fail = 1-8/math.pi**2+2*qft_error
    for repetitions in range(1, 10001, 2):
        tail = sum(math.comb(repetitions, j)*fail**j*(1-fail)**(repetitions-j)
                   for j in range((repetitions+1)//2, repetitions+1))
        if tail <= delta:
            return repetitions, tail
    raise ValueError("Failure target outside calculator range")


def taylor_plan(n: int, b: int) -> dict:
    # Uniform cosine error <= 2^-s. Squared-predicate amplitude error then
    # <= L sqrt(2^-s) <= 2^(-b-5), before comparator discretization.
    s = 2*b+2*n+10
    d = 0
    while (2*d+2)*math.log(1.6)-math.lgamma(2*d+3) > -(s+2)*math.log(2):
        d += 1
    # Horner argument t<=3, each truncation/coefficient error <= 2^-f.
    f = s+math.ceil((d+2)*math.log2(3))+8
    # Uniform arithmetic width also covers exact squared predicates and
    # the exact 3n-bit p field. This is intentionally much wider than needed.
    word = f+2*b+6*n+16
    return dict(cosine_error_bits=s, cosine_degree_in_x_squared=d,
                fraction_bits=f, uniform_arithmetic_word_bits=word)


def operation_budget(n: int, prec: dict, arithmetic: dict) -> dict:
    b, h, p = prec["amplitude_bits"], prec["geometric_cdf_bits"], prec["neumann_degree"]
    w, d = arithmetic["uniform_arithmetic_word_bits"], arithmetic["cosine_degree_in_x_squared"]
    cell = 3*n-1
    pbits = 3*n
    a, m, cmp = add(w), multiply(w), compare(w)
    # One-way cosine computation retains Horner iterates. Multiplication
    # product copying and uncomputation are included in each Horner step.
    cos_forward = (2*d+4)*m + (d+8)*a + Gates(boolean=(4*d+20)*w) + 4*cmp
    # Three cosines; exact lambda assembly; predicate; complete uncompute.
    frequency_predicate = 6*cos_forward + 6*m + 24*a + 8*cmp + Gates(hadamard=2*b)
    # x -> (5x+3) mod V and bit reversal; p numerator has exact pbits width.
    p_program = 3*add(pbits)+Gates(boolean=4*pbits)
    gamma_predicate = 2*p_program+4*m+4*a+4*cmp+Gates(hadamard=2*b)
    qft = Gates(hadamard=3*n, phase=3*n*(n-1)//2)
    # e^(-ik_z)-1 polar phase uses n+2 dyadic/sign phases; no inverse trig.
    W = gamma_predicate+frequency_predicate+qft+Gates(boolean=20*w, phase=n+2)
    # Highest-set-bit shell, 64 outer-cube selector, conditional coordinate
    # Hadamards, shell recovery and masks. Bound covers <=16(n+1) comparisons
    # and <=64(n+1) word operations; no list of L^3 modes.
    shell_route = 16*(n+1)*cmp+Gates(boolean=64*(n+1)*w, hadamard=12*n+16)
    position_phases = Gates(phase=3*n*n)
    primitive_q = frequency_predicate+shell_route+position_phases
    # Exact convex mixture p q_s+(1-p)q_(s+ez), with uniform-register
    # projection. A controlled site increment and its inverse are charged.
    density_q = primitive_q+2*p_program+4*compare(pbits)+4*a+Gates(hadamard=2*pbits)
    # Generous routing budget: 128 word comparisons, 256 word controlled
    # swaps (3 Boolean gates/bit), 24 neighbor shifts; include cleanup.
    rank_route = 128*cmp+Gates(boolean=768*w)+48*a
    # Four operator types, finite dyadic coefficients, all ranks/masks.
    top_lcu = cdf_select(4, b+8, w)
    # Real averaging / Hermitian off-diagonal dilation charged explicitly:
    # direct <=4 q; exchange <=4 W; K+K† <=4 q+4 W.
    Mblock = 8*density_q+8*W+rank_route+top_lcu
    # Real Gram source: two W plus conjugate alternative, collision/order mask.
    Bblock = 4*W+rank_route+Gates(hadamard=cell+1, boolean=cell)
    finite_lcu = cdf_select(p+1, h, w)
    power_controls = 2*p*compare(max(w, (p+1).bit_length()))
    scalar = 2*Bblock+p*Mblock+finite_lcu+power_controls+Gates(hadamard=2)
    # All nested conditions can be combined in clean predicate ancillas.
    # At most six extra controls are padded on EVERY elementary primitive.
    # An X with at most eight controls costs <=13 Toffoli using clean ancillas;
    # controlled H additionally <=2 T; controlled phase <=3 P rotations.
    ccx = 13*(scalar.boolean+scalar.hadamard+scalar.phase)
    exact_T = 2*scalar.hadamard
    phase_count = 3*scalar.phase
    # Fixed compiled access circuits: allocate synthesis error to their
    # block, not once per AE invocation. Inverses are literal circuit reversals.
    phase_per_primitive = max(W.phase, density_q.phase, 1)
    phase_tolerance = prec["primitive_tolerance"]/(32*phase_per_primitive)
    # Fresh signal banks for every projected product. Clean arithmetic scratch
    # is shared; retained Horner intermediates are local to one primitive.
    aW = 2*b+8*n+16
    aq = b+pbits+8*n+32
    route_flags = 128+24*cell
    aM = 8*aq+8*aW+route_flags+(b+8)
    aB = 4*aW+route_flags
    scratch = (3*(d+12)+20)*w+16*w+64
    width = 4*cell+2*aB+p*aM+h+(p+1).bit_length()+scratch+16
    return dict(
        primitive_operation_budgets={"cosine_forward":asdict(cos_forward),
            "frequency_predicate":asdict(frequency_predicate),
            "gamma_predicate":asdict(gamma_predicate), "W":asdict(W),
            "density_half_Green":asdict(density_q), "M_block":asdict(Mblock),
            "B_source":asdict(Bblock), "scalar_unitary":asdict(scalar)},
        calls_per_scalar={"B_source":2, "M_block":p,
                          "W_upper":8+8*p, "density_half_Green_upper":8*p},
        analytic_controlled_scalar_budget={"Toffoli_upper":ccx,
            "additional_exact_T_upper":exact_T,
            "unsynthesized_phase_rotations_upper":phase_count,
            "phase_operator_error_per_rotation":phase_tolerance},
        analytic_logical_qubit_upper=width,
        fresh_signal_qubits_per_M=aM, fresh_signal_qubits_per_B=aB,
        reusable_arithmetic_qubits=scratch)


def decimal_pi() -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 100
        def atan_inverse(x: int) -> Decimal:
            z=Decimal(1)/x
            term=z
            value=term
            j=1
            while abs(term)>Decimal("1e-105"):
                term *= -z*z
                value += term/(2*j+1)
                j += 1
            return value
        return +(16*atan_inverse(5)-4*atan_inverse(239))


def geometric_cdf(beta: float, p: int, h: int) -> dict:
    # A small classical table, not an amplitude oracle. Its h-bit interval
    # endpoints are the actual constants used by the proposed SELECT circuit.
    # Decimal keeps rounding noise well below one endpoint unit in our examples.
    with localcontext() as ctx:
        ctx.prec = 90
        pi=decimal_pi()
        Z=8*pi
        g=Decimal(".01"); E0=Decimal(2); gamma=Decimal(".5")
        betaD=Decimal(".12")+3*g*Z/E0+2*(4*pi*g*gamma**2)/E0+(4*g*gamma*(2*pi*Z).sqrt())/(Decimal(2).sqrt()*E0)
        weights=[betaD**k for k in range(p+1)]
        normalizer=sum(weights)
        total=Decimal(0)
        scale=1<<h
        endpoints=[0]
        for weight in weights[:-1]:
            total += weight/normalizer
            endpoints.append(int((scale*total).to_integral_value(rounding="ROUND_HALF_EVEN")))
        endpoints.append(scale)
        l1=sum(abs(Decimal(endpoints[k+1]-endpoints[k])/scale-weights[k]/normalizer) for k in range(p+1))
    return dict(endpoint_bits=h, integer_endpoints=endpoints,
                observed_L1_weight_error=float(l1),
                analytic_L1_weight_error_upper=2*(p+1)/scale,
                constant_table_scope="Default parameters only; 90-digit Decimal evaluation, not interval-certified compilation")


def instance(L: int, epsilon: float, failure: float) -> dict:
    if L<4 or L&(L-1):
        raise ValueError("L must be a power of two >=4")
    n=L.bit_length()-1
    c=constants()
    prec=select_precision(c, epsilon)
    arithmetic=taylor_plan(n, prec["amplitude_bits"])
    ops=operation_budget(n, prec, arithmetic)
    reps, failbound=confidence_repetitions(failure)
    M=prec["AE_M"]
    preparations=reps*(2*M-1)
    grover=reps*(M-1)
    ae_bits=M.bit_length()-1
    width=ops["analytic_logical_qubit_upper"]+ae_bits
    # Controlled zero/good reflections: multi-control AND ladder and inverse;
    # deliberately overcharge 8*width+32 Toffoli per Grover iteration.
    reflection_ccx=grover*(8*width+32)
    scalar=ops["analytic_controlled_scalar_budget"]
    # Inverse QFT on the AE counter once per repetition.
    ae_phase=reps*ae_bits*(ae_bits-1)//2
    return dict(L=L, V=L**3//2, epsilon=epsilon, target_failure=failure,
                constants=c, precision=prec, arithmetic=arithmetic,
                geometric_weight_table=geometric_cdf(c["beta"],prec["neumann_degree"],prec["geometric_cdf_bits"]),
                **ops,
                amplitude_estimation=dict(repetitions=reps,
                    median_failure_upper=failbound, M=M,
                    inverse_QFT_operator_error=1e-6,
                    inverse_QFT_rotation_error=1e-6/max(1,ae_bits*(ae_bits-1)//2),
                    Hadamard_test_preparations=preparations,
                    Grover_iterations=grover,
                    W_calls_upper=preparations*ops["calls_per_scalar"]["W_upper"],
                    density_half_Green_calls_upper=preparations*ops["calls_per_scalar"]["density_half_Green_upper"]),
                complete_analytic_budget=dict(
                    logical_qubits_upper=width,
                    Toffoli_upper=preparations*scalar["Toffoli_upper"]+reflection_ccx,
                    additional_exact_T_upper=preparations*scalar["additional_exact_T_upper"],
                    unsynthesized_phase_rotations_upper=preparations*scalar["unsynthesized_phase_rotations_upper"]+ae_phase,
                    rotation_synthesis_status="Not run; no compiled T count asserted",
                    Clifford_T_conversion="7*Toffoli + additional_exact_T + sum(actual synthesized phase T counts); clean ancillary Toffoli convention",
                    setup="Polynomial-size arithmetic/phase and small CDF constant compilation; not timed. Effective model is supplied, not derived from nuclei.",
                    output="One scalar and classical median over repetitions; no amplitudes or full SD vector read out"),
                total_scalar_error_budget=prec["tail_error"]+prec["primitive_scalar_error"]+prec["geometric_scalar_error"]+prec["AE_energy_error_bound"])


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", nargs="+",type=int,default=[4,8,16,64,256])
    parser.add_argument("--epsilon",type=float,default=constants()["default_epsilon"])
    parser.add_argument("--failure",type=float,default=.01)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    report=dict(status="analytic_upper_budgets_not_compiled_circuit_counts",
                version=1, coefficient_program="p_x=1/4+bit_reverse_(3log2L-1)((5x+3) mod V)/(2V)",
                assumptions=["constant onsite E0=2", "g=.01, t_o=t_v=.02 bounds, gamma_max=.5",
                    "exact procedural p; no arbitrary table or free QRAM", "Z=8pi rotation-free shell reference",
                    "fresh success ancillas in all projected products", "six extra controls padded on each primitive gate",
                    "failure probability applies to ideal logical execution of fixed compiled access circuits",
                    "no fault-tolerance, connectivity, magic-state throughput, or runtime model"],
                instances=[instance(L,args.epsilon,args.failure) for L in args.L])
    out=json.dumps(report,indent=2)+"\n"
    if args.output:
        args.output.write_text(out)
    else:
        print(out,end="")


if __name__=="__main__":
    main()
