#!/usr/bin/env python3
"""Exact full-SD variational energy/residual from a local-pair trial.

The trial is diagonal-MP2, but the residual includes the complete SD matrix:
direct/exchange, local electron/hole hopping, and K† into singles. Requires
O(V²) memory and O(V² log V) arithmetic, without an O(V⁴) SD vector.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
import os
import platform
from pathlib import Path
import time
import numpy as np


def coefficients(L):
    assert L>=2 and L&(L-1)==0
    shape=(L,L,L//2);V=math.prod(shape);w=int(math.log2(V))
    x=np.arange(V,dtype=np.int64);a=(5*x+3)&(V-1);rev=np.zeros(V,dtype=np.int64)
    for k in range(w):rev|=((a>>k)&1)<<(w-1-k)
    p=.25+rev/(2.*V)
    return shape,p,np.sqrt(p*(1-p))


def pair_trial(L,g=.01,E0=2.,hopping=.02,return_arrays=False):
    begin=time.perf_counter()
    shape,p,gamma=coefficients(L);V=len(p)
    coords=np.array(list(itertools.product(*[range(n)for n in shape])))
    k=2*np.pi*np.arange(L)/L
    lam=4*(np.sin(k[:,None,None]/2)**2+np.sin(k[None,:,None]/2)**2+np.sin(k[None,None,:]/2)**2)
    ilam=np.zeros_like(lam);np.divide(1.,lam,out=ilam,where=lam>1e-14)
    green=(4*np.pi*g*np.fft.ifftn(ilam)).real
    # Displacement kernel between unit internal dipoles, on the cell lattice.
    gz=2*np.arange(L//2)
    base=2*green[:,:,gz]-green[:,:,(gz+1)%L]-green[:,:,(gz-1)%L]
    r=(coords[:,None,:]-coords[None,:,:])%np.array(shape)
    B=base[r[:,:,0],r[:,:,1],r[:,:,2]]
    J=gamma[:,None]*B*gamma[None,:]
    T=J/(2*E0);np.fill_diagonal(T,0.)
    # J T using the primitive translation structure; J itself is nonuniform.
    fft_axes=(0,1,2)
    x=(gamma[:,None]*T).reshape((*shape,V))
    JT=(gamma[:,None]*np.fft.ifftn(np.fft.fftn(x,axes=fft_axes)*np.fft.fftn(base)[...,None],axes=fft_axes).real.reshape(V,V))
    kappa=2*p-1
    self_vo=2*(gamma**2)*green[0,0,0]+(1-2*gamma**2)*green[0,0,1]
    direct=kappa[:,None]*kappa[None,:]*B-self_vo[:,None]-self_vo[None,:]
    interaction=direct*T+JT+JT.T
    np.fill_diagonal(interaction,0.)
    local_pair_residual=-interaction
    # Single local excitation x receives all inserted-pair predecessors y.
    singles_residual=-kappa*np.sum(gamma[:,None]*B*T,axis=0)
    local_pair_normsq=.5*float(np.sum(local_pair_residual**2))
    singles_normsq=float(singles_residual@singles_residual)
    charge_normsq=0.;degree=0
    # Each charged double has exactly one local pair; its spectator pair labels
    # the output uniquely. Electron and hole routes can collide and are added
    # coherently before squaring: fv=-hopping, fo=+hopping.
    for axis,length in enumerate(shape):
        if length==1:continue
        directions=(1,)if length==2 else(-1,1)
        for step in directions:
            degree+=1
            shifted=coords.copy();shifted[:,axis]=(shifted[:,axis]+step)%length
            neighbors=np.ravel_multi_index(shifted.T,shape)
            amplitude=-hopping*(T+T[neighbors,:])
            amplitude[np.arange(V),np.arange(V)]=0.
            amplitude[np.arange(V),neighbors]=0.
            charge_normsq+=float(np.sum(amplitude**2))
    residual_normsq=local_pair_normsq+singles_normsq+charge_normsq
    mp2=.5*float(np.sum(J*T))/V
    variational=mp2-.5*float(np.sum(T*interaction))/V
    G0=float(green[0,0,0]/g)
    Gbar=3*np.pi+np.pi**3/24
    def bound(G):return g*(6*G+4*np.pi+2*np.sqrt(2*np.pi*G))
    Eg=E0-2*degree*hopping
    delta_uniform=Eg-bound(Gbar)
    delta_finite=Eg-bound(G0)
    # Rank preconditioning: D=E0 on singles and 2E0 on doubles.  The physical
    # contraction follows from rank-dependent pair/exchange/hopping counts;
    # unlike the finite gap above, it uses only an all-volume Green bound.
    # This formula requires the constant onsite separation used in this fixture.
    gamma_star=.5
    alpha_J=4*np.pi*g*gamma_star**2
    K_norm_bound=4*g*gamma_star*np.sqrt(2*np.pi*Gbar)
    q_rank_uniform=(6*(2*hopping)+3*g*Gbar+2*alpha_J+K_norm_bound/np.sqrt(2))/E0
    delta_rank_uniform=E0*(1-q_rank_uniform)
    assert delta_uniform>0
    assert delta_rank_uniform>0
    result={'L':L,'V':V,'unrestricted_SD_dimension':V*V+(V*(V-1)//2)**2,
            'parameters':{'g':g,'E0':E0,'hopping':hopping,'max_graph_degree':degree},
            'coefficient_program':'p_x=1/4+bit_reverse_w((5*x+3) mod V)/(2V)',
            'MP2_minus_energy_density':mp2,'variational_minus_energy_density':variational,
            'residual_norm_squared_channels':{'local_pair_doubles':local_pair_normsq,'local_singles':singles_normsq,'charge_transfer_doubles':charge_normsq,'total':residual_normsq},
            'actual_finite_Green_diagonal':G0,'Delta_uniform':delta_uniform,'Delta_paid_finite_Green':delta_finite,
            'rank_preconditioned_uniform_q':q_rank_uniform,'Delta_rank_uniform':delta_rank_uniform,
            'energy_error_certificate_uniform':residual_normsq/(V*delta_uniform),
            'energy_error_certificate_paid_finite_Green':residual_normsq/(V*delta_finite),
            'energy_error_certificate_rank_uniform':residual_normsq/(V*delta_rank_uniform),
            'seconds_total':time.perf_counter()-begin,
            'pair_matrix_float64_bytes':8*V*V,
            'scope':'Full finite SD target certified by exact residual identity; displayed arithmetic evaluated in double precision, not interval arithmetic; no physical material or quantum crossover claim.'}
    if return_arrays:return result,T,local_pair_residual,singles_residual
    return result


def validate():
    from benchmark_lccsd import Model
    checks=[]
    for L in(2,4):
        result,T,rd,rs=pair_trial(L,return_arrays=True)
        model=Model(L)
        trial=model.b/model.denominator
        residual=model.b-model.matvec(trial)
        full_var=float((2*model.b@trial-trial@model.matvec(trial))/model.V)
        residual_norm=float(residual@residual)
        pair_locations=np.arange(model.P)*model.P+np.arange(model.P)
        local_expected=residual[model.ns:][pair_locations]
        local_actual=rd[model.pairs[:,0],model.pairs[:,1]]
        singles_expected=residual[:model.ns].reshape(model.V,model.V)
        item={'L':L,'variational_scalar_error':abs(result['variational_minus_energy_density']-full_var),
              'residual_norm_squared_error':abs(result['residual_norm_squared_channels']['total']-residual_norm),
              'local_pair_component_max_error':float(np.max(abs(local_actual-local_expected))),
              'single_component_max_error':float(np.max(abs(np.diag(singles_expected)-rs))),
              'offdiagonal_single_max':float(np.max(abs(singles_expected-np.diag(np.diag(singles_expected)))))}
        assert max(v for k,v in item.items()if k!='L')<1e-12,item
        checks.append(item)
    return checks


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--L',type=int,nargs='+',default=[4,8,16]);args=parser.parse_args()
    data={'status':'complete_SD_variational_scalar_and_residual_without_SD_vector','numpy':np.__version__,
          'environment':{'python':platform.python_version(),'numpy':np.__version__,
                         'platform':platform.platform(),'machine':platform.machine(),
                         'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS')},
          'structural_validation':validate(),'target_epsilon':5.698089085861592e-9,
          'instances':[pair_trial(L)for L in args.L]}
    for r in data['instances']:
        r['meets_target_with_uniform_gap']=bool(r['energy_error_certificate_uniform']<=data['target_epsilon'])
        r['meets_target_with_paid_finite_gap']=bool(r['energy_error_certificate_paid_finite_Green']<=data['target_epsilon'])
        r['meets_target_with_rank_uniform_gap']=bool(r['energy_error_certificate_rank_uniform']<=data['target_epsilon'])
    (Path(__file__).parent / 'results' / 'pair_variational_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))


if __name__=='__main__':main()
