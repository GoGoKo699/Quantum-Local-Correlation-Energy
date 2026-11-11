#!/usr/bin/env python3
"""Matched classical solvers for the explicit nonuniform Poisson LCCSD model.

No full-electron Hilbert space or dense SD matrix is used in the benchmarks.
An independent full-fermion comparison is used only at L=2 for validation.
The p_x program is a bit-reversed affine permutation, with no random table.
"""
from __future__ import annotations
import argparse
import importlib.util
import itertools
import json
import math
import os
import platform
from pathlib import Path
import time

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, eigsh


def reverse_bits(x: np.ndarray, width: int) -> np.ndarray:
    out = np.zeros_like(x)
    for k in range(width):
        out |= ((x >> k) & 1) << (width - 1 - k)
    return out


class Model:
    def __init__(self, L=4, g=0.01, E0=2.0, hopping=0.02):
        started = time.perf_counter()
        assert L >= 2 and L & (L-1) == 0
        self.L, self.g, self.E0, self.hopping = L, g, E0, hopping
        self.shape = (L, L, L//2)
        self.xyz = np.array(list(itertools.product(range(L), range(L), range(L//2))))
        self.V = V = len(self.xyz)
        self.width = int(math.log2(V))
        ids = np.arange(V, dtype=np.int64)
        # An exact O(log V)-bit coefficient program; no stored disorder seed/list.
        affine = (5*ids+3) & (V-1)
        self.p = 0.25 + reverse_bits(affine, self.width)/(2.0*V)
        self.gamma = np.sqrt(self.p*(1-self.p))
        primitive_k = 2*np.pi*np.arange(L)/L
        lam = 4*(np.sin(primitive_k[:,None,None]/2)**2+
                 np.sin(primitive_k[None,:,None]/2)**2+
                 np.sin(primitive_k[None,None,:]/2)**2)
        inverse_lam = np.zeros_like(lam)
        np.divide(1.0, lam, out=inverse_lam, where=lam>1e-14)
        kernel_started = time.perf_counter()
        self.green = (4*np.pi*g*np.fft.ifftn(inverse_lam)).real
        self.kernel_fft_seconds = time.perf_counter()-kernel_started
        r = self.xyz[:,None,:]-self.xyz[None,:,:]
        r[:,:,2] *= 2
        C = {}
        for u,v in itertools.product(range(2),repeat=2):
            C[u,v] = self.green[r[:,:,0]%L,r[:,:,1]%L,(r[:,:,2]+u-v)%L]
        wo = np.stack([self.p,1-self.p])
        wv = np.stack([1-self.p,self.p])
        def charges(w1,w2):
            return sum(w1[u,:,None]*w2[v,None,:]*C[u,v] for u,v in itertools.product(range(2),repeat=2))
        self.Vvv, self.Voo, self.Vvo = charges(wv,wv),charges(wo,wo),charges(wv,wo)
        self.Jbase = C[0,0]+C[1,1]-C[0,1]-C[1,0]
        self.J = self.gamma[:,None]*self.Jbase*self.gamma[None,:]
        self.Qv = self.gamma[:,None]*((C[1,0]-C[0,0])*wv[0,None,:]+(C[1,1]-C[0,1])*wv[1,None,:])
        self.Qo = self.gamma[:,None]*((C[1,0]-C[0,0])*wo[0,None,:]+(C[1,1]-C[0,1])*wo[1,None,:])
        self.Lcoeff = self.Qv[:,:,None]-self.Qo[:,None,:]  # inserted pair, electron, hole
        # Simple constant nonzero hopping. Local orbital inhomogeneity still breaks translation.
        rows,cols=[],[]
        for x,c in enumerate(self.xyz):
            neighbors=set()
            for axis,length in enumerate(self.shape):
                for sign in (-1,1):
                    d=c.copy(); d[axis]=(d[axis]+sign)%length
                    y=int(np.ravel_multi_index(tuple(d),self.shape))
                    if y!=x:neighbors.add(y)
            for y in sorted(neighbors):rows.append(x);cols.append(y)
        adjacency=sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(V,V))
        self.fo=-0.5*E0*sparse.eye(V,format='csr')+hopping*adjacency
        self.fv=0.5*E0*sparse.eye(V,format='csr')-hopping*adjacency
        self.pairs=np.array(list(itertools.combinations(range(V),2)),dtype=np.int32)
        self.P=len(self.pairs)
        self.pair_index=np.full((V,V),-1,dtype=np.int32)
        for q,(a,b) in enumerate(self.pairs):self.pair_index[a,b]=self.pair_index[b,a]=q
        self.fo2=self.exterior_two(self.fo)
        self.fv2=self.exterior_two(self.fv)
        self.ns=V*V; self.nd=self.P*self.P; self.dim=self.ns+self.nd
        a,b=self.pairs[:,0,None],self.pairs[:,1,None]
        i,j=self.pairs[None,:,0],self.pairs[None,:,1]
        self.direct1=-self.Vvo
        self.direct2=(self.Vvv[a,b]+self.Voo[i,j]-self.Vvo[a,i]-self.Vvo[a,j]-self.Vvo[b,i]-self.Vvo[b,j])
        erows,ecols,evals=[],[],[]
        for y,a,i in itertools.product(range(V),repeat=3):
            if y==a or y==i:continue
            erows.append((y*V+a)*V+i)
            ecols.append(int(self.pair_index[y,a]*self.P+self.pair_index[y,i]))
            evals.append(1 if (y<a)==(y<i) else -1)
        self.E=sparse.csr_matrix((evals,(erows,ecols)),shape=(V**3,self.nd),dtype=float)
        self.b=np.zeros(self.dim)
        b2=self.b[self.ns:].reshape(self.P,self.P)
        b2[np.arange(self.P),np.arange(self.P)]=self.J[self.pairs[:,0],self.pairs[:,1]]
        self.denominator=np.r_[np.full(self.ns,E0),np.full(self.nd,2*E0)]
        self.invroot=1/np.sqrt(self.denominator)
        self.A=LinearOperator((self.dim,self.dim),matvec=self.matvec,rmatvec=self.matvec,dtype=float)
        Gbar=3*np.pi+np.pi**3/24
        self.gamma_bound=g*(6*Gbar+16*np.pi*.25+2*np.sqrt(2*np.pi*Gbar))
        self.row_degree=int(np.diff(adjacency.indptr).max())
        self.Eg=E0-2*hopping*self.row_degree
        self.delta=self.Eg-self.gamma_bound
        self.q_uniform=(2*hopping*self.row_degree+self.gamma_bound)/E0
        self.B0sq=(4*np.pi*g*.25)**2/2
        self.c=self.invroot*self.b
        self.BDsq=float(self.c@self.c/V)
        self.preconditioned_abs_bound=self.compute_absolute_bound()
        self.setup_seconds=time.perf_counter()-started

    def exterior_two(self,f):
        rows,cols,values=[],[],[]
        f=f.tocsc()
        for col,(a,b) in enumerate(self.pairs):
            for old,other,position in ((a,b,0),(b,a,1)):
                for at in range(f.indptr[old],f.indptr[old+1]):
                    new=int(f.indices[at]); value=f.data[at]
                    if new==other:continue
                    first,second=(new,other) if position==0 else (other,new)
                    row=int(self.pair_index[first,second])
                    sign=1 if first<second else -1
                    rows.append(row);cols.append(col);values.append(sign*value)
        return sparse.csr_matrix((values,(rows,cols)),shape=(self.P,self.P))

    def matvec(self,x):
        V=self.V
        s=x[:self.ns].reshape(V,V)
        d=x[self.ns:].reshape(self.P,self.P)
        ys=self.fv@s-(self.fo@s.T).T+self.direct1*s
        ys[np.arange(V),np.arange(V)]+=self.J@np.diag(s)
        yd=self.fv2@d-(self.fo2@d.T).T+self.direct2*d
        incidence=(self.E@d.ravel()).reshape(V,V,V)
        yd+=(self.E.T@(self.J@incidence.reshape(V,-1)).ravel()).reshape(self.P,self.P)
        yd+=(self.E.T@(self.Lcoeff*s[None,:,:]).ravel()).reshape(self.P,self.P)
        ys+=np.sum(self.Lcoeff*incidence,axis=0)
        return np.r_[ys.ravel(),yd.ravel()]

    def preconditioned(self,x):
        return self.invroot*self.matvec(self.invroot*x)

    def M(self,x):
        return x-self.preconditioned(x)

    def compute_absolute_bound(self):
        # Triangle bound on |I-D^-1/2 A D^-1/2|, fully computable from paid finite data.
        V=self.V
        h1v=self.fv-sparse.diags(self.fv.diagonal())
        h1o=self.fo-sparse.diags(self.fo.diagonal())
        h2v=self.fv2-sparse.diags(self.fv2.diagonal())
        h2o=self.fo2-sparse.diags(self.fo2.diagonal())
        absE=abs(self.E)
        s=np.ones((V,V));d=np.ones((self.P,self.P))
        rs=(abs(h1v)@s+(abs(h1o)@s.T).T+abs(self.direct1))/self.E0
        rs[np.arange(V),np.arange(V)]+=np.sum(abs(self.J),axis=1)/self.E0
        rd=(abs(h2v)@d+(abs(h2o)@d.T).T+abs(self.direct2))/(2*self.E0)
        inc=(absE@d.ravel()).reshape(V,V,V)
        rd+=(absE.T@(abs(self.J)@inc.reshape(V,-1)).ravel()).reshape(self.P,self.P)/(2*self.E0)
        rd+=(absE.T@abs(self.Lcoeff).ravel()).reshape(self.P,self.P)/(np.sqrt(2)*self.E0)
        rs+=np.sum(abs(self.Lcoeff)*inc,axis=0)/(np.sqrt(2)*self.E0)
        return float(max(rs.max(),rd.max()))

    def mp2_fft(self):
        # Ordered pairs produce a factor 1/2, denominator is 2E0.
        kernel=self.Jbase[0].reshape(self.shape).copy()**2
        kernel.flat[0]=0.0
        w=(self.gamma**2).reshape(self.shape)
        conv=np.fft.ifftn(np.fft.fftn(kernel)*np.fft.fftn(w)).real
        return float(np.sum(w*conv)/(4*self.E0*self.V))

    def output_memory_bytes(self):
        arrays=[self.p,self.gamma,self.green,self.Vvv,self.Voo,self.Vvo,self.J,self.Jbase,self.Qv,self.Qo,self.Lcoeff,self.pairs,self.pair_index,self.direct1,self.direct2,self.b,self.denominator,self.invroot,self.c]
        return int(sum(x.nbytes for x in arrays)+sum(x.data.nbytes+x.indices.nbytes+x.indptr.nbytes for x in (self.fo,self.fv,self.fo2,self.fv2,self.E)))


def validate_small():
    # Independent second-quantized H in the full fixed-electron sector.
    module_path=Path(__file__).with_name('fermion_reference.py')
    spec=importlib.util.spec_from_file_location('full_fermions',module_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    m=Model(2)
    N=2*m.V
    xyz=np.array(list(itertools.product(range(2),repeat=3)))
    diff=(xyz[:,None,:]-xyz[None,:,:])%2
    C=m.green[diff[:,:,0],diff[:,:,1],diff[:,:,2]]
    orb=np.zeros((N,N))
    for x,p in enumerate(m.p):
        orb[2*x,2*x]=np.sqrt(p);orb[2*x+1,2*x]=np.sqrt(1-p)
        orb[2*x,2*x+1]=-np.sqrt(1-p);orb[2*x+1,2*x+1]=np.sqrt(p)
    rho=np.einsum('xp,xq->xpq',orb,orb)
    eri=np.einsum('xpq,xy,yrs->pqrs',rho,C,rho,optimize=True)
    f=np.zeros((N,N));f[0::2,0::2]=m.fo.toarray();f[1::2,1::2]=m.fv.toarray()
    h=f.copy()
    for i in range(0,N,2):h-=eri[:,:,i,i]-eri[:,i,i,:]
    states,pos,H=module.hamiltonian({'norb':N,'nelec':m.V,'one_body':h,'eri':eri})
    ref=sum(1<<i for i in range(0,N,2));ri=pos[ref]
    excitations=[]
    # Order exactly the single matrix and pair-electron/pair-hole tensor used above.
    for a,i in itertools.product(range(m.V),repeat=2):
        excitations.append(module.act(ref,module.excitation_operations((2*i,),(2*a+1,))))
    for (a,b),(i,j) in itertools.product(m.pairs,m.pairs):
        excitations.append(module.act(ref,module.excitation_operations((2*int(i),2*int(j)),(2*int(a)+1,2*int(b)+1))))
    ix=np.array([pos[s]for s,sign in excitations]); signs=np.array([sign for s,sign in excitations])
    exact=signs[:,None]*H[np.ix_(ix,ix)]*signs[None,:]-H[ri,ri]*np.eye(m.dim)
    actual=np.column_stack([m.matvec(np.eye(m.dim)[:,j])for j in range(m.dim)])
    source=signs*H[ix,ri]
    result={'full_fixed_electron_dimension':len(states),'SD_dimension':m.dim,
            'matrix_max_error':float(np.max(abs(exact-actual))),
            'matrix_spectral_error':float(np.linalg.norm(exact-actual,2)),
            'source_max_error':float(np.max(abs(source-m.b))),
            'symmetry_error':float(np.linalg.norm(actual-actual.T)),
            'MP2_fft_error':abs(m.mp2_fft()-m.BDsq)}
    assert max(result[k]for k in result if 'error'in k)<1e-11,result
    return result


def run(L, tol=1e-10, max_degree=30):
    m=Model(L)
    assert m.delta>0
    mv_started=time.perf_counter();m.matvec(m.b)
    matvec_seconds=time.perf_counter()-mv_started
    count=0
    def callback(x):
        nonlocal count
        count+=1
    started=time.perf_counter()
    pre=LinearOperator((m.dim,m.dim),matvec=lambda x:x/m.denominator,dtype=float)
    solution,info=cg(m.A,m.b,rtol=1e-12,atol=0,M=pre,callback=callback,maxiter=200)
    cg_seconds=time.perf_counter()-started
    assert info==0,info
    residual=m.b-m.matvec(solution)
    exact_proxy=float(m.b@solution/m.V)
    # Variational energy estimate with a rigorous residual error certificate.
    variational=float((2*m.b@solution-solution@m.matvec(solution))/m.V)
    residual_certificate=float(residual@residual/(m.V*m.delta))
    mp2_started=time.perf_counter();mp2=m.mp2_fft();mp2_seconds=time.perf_counter()-mp2_started
    term=m.c.copy();partial=0.;rows=[];trial_sum=np.zeros_like(m.c)
    degree_started=time.perf_counter()
    for p in range(max_degree+1):
        partial+=float(m.c@term/m.V)
        trial_sum+=term
        trial=m.invroot*trial_sum
        trial_residual=m.b-m.matvec(trial)
        trial_variational=float((2*m.b@trial-trial@m.matvec(trial))/m.V)
        trial_certificate=float(trial_residual@trial_residual/(m.V*m.delta))
        rem_uniform=m.BDsq*m.q_uniform**(p+1)/(1-m.q_uniform)
        q=m.preconditioned_abs_bound
        rem_finite=m.BDsq*q**(p+1)/(1-q)if q<1 else None
        rows.append({'degree':p,'scalar':partial,'absolute_error_vs_CG':abs(partial-variational),
                     'uniform_bias_certificate':rem_uniform,'finite_absolute_row_bias_certificate':rem_finite,
                     'trial_variational_scalar':trial_variational,'trial_variational_error_vs_CG':abs(trial_variational-variational),
                     'trial_variational_residual_certificate':trial_certificate})
        if rem_finite is not None and rem_finite<tol and p>=3:break
        term=m.M(term)
    neumann_seconds=time.perf_counter()-degree_started
    # The MP2 control itself is a deterministic FFT, stronger than sampling it.
    weights=(m.gamma**2)
    pair_value=m.J**2/(2*m.E0)
    np.fill_diagonal(pair_value,0)
    root_values=np.sum(pair_value,axis=1)/2
    # Uniform anchor sampling after exact conditional displacement summation.
    source_variance=float(np.var(root_values))
    # A homogeneous medium with the actual average gamma^2 supplies a local control.
    reference=float(weights.mean())
    homogeneous=np.sum(m.Jbase**2-np.diag(np.diag(m.Jbase**2)),axis=1)*reference**2/(4*m.E0)
    control_variance=float(np.var(root_values-homogeneous))
    squared_kernel=m.Jbase**2
    np.fill_diagonal(squared_kernel,0)
    squared_kernel_sum=float(np.sum(squared_kernel[0]))
    pair_probability=squared_kernel/(m.V*squared_kernel_sum)
    scalar_factor=squared_kernel_sum/(4*m.E0)
    pair_estimator=scalar_factor*weights[:,None]*weights[None,:]
    known_weight_mean=11/48-1/(24*m.V**2)
    assert abs(weights.mean()-known_weight_mean)<1e-14
    centered_estimator=scalar_factor*(weights[:,None]-known_weight_mean)*(weights[None,:]-known_weight_mean)
    pair_mean=float(np.sum(pair_probability*pair_estimator))
    centered_mean=float(np.sum(pair_probability*centered_estimator))
    pair_variance=float(np.sum(pair_probability*pair_estimator**2)-pair_mean**2)
    centered_variance=float(np.sum(pair_probability*centered_estimator**2)-centered_mean**2)
    result={'L':L,'V':m.V,'dimensions':{'singles':m.ns,'doubles':m.nd,'total':m.dim},
            'coefficient_program':'p_x=1/4+bit_reverse_w((5*x+3) mod V)/(2V), w=log2 V',
            'parameters':{'g':m.g,'E0':m.E0,'hopping':m.hopping,'max_graph_degree':m.row_degree},
            'certificates':{'Delta_uniform':m.delta,'Gamma_uniform':m.gamma_bound,'q_uniform':m.q_uniform,
                            'q_finite_absolute_rows':m.preconditioned_abs_bound,'B0_squared':m.B0sq,'BD_squared_actual':m.BDsq},
            'reference':{'minus_energy_per_cell':variational,'CG_linear_scalar':exact_proxy,'CG_iterations':count,
                         'residual_l2':float(np.linalg.norm(residual)),'variational_residual_error_bound':residual_certificate},
            'mp2':{'minus_energy_per_cell':mp2,'absolute_bias_vs_SD':abs(mp2-variational),
                   'relative_bias_vs_SD':abs(mp2-variational)/abs(variational),
                   'fft_vs_direct_contraction_error':abs(mp2-m.BDsq),
                   'anchor_estimator_variance':source_variance,'homogeneous_control_variance':control_variance,
                   'pair_importance_estimator_mean':pair_mean,
                   'pair_importance_variance':pair_variance,
                   'known_control_mean':scalar_factor*known_weight_mean**2,
                   'centered_pair_variance':centered_variance,
                   'pair_variance_reduction_factor':pair_variance/centered_variance,
                   'note':'The exact same MP2 number is obtained deterministically by FFT. Anchor variance does not imply that sampling is the best comparator.'},
            'neumann':rows,
            'seconds':{'input_and_operator_setup':m.setup_seconds,'primitive_Green_FFT_only':m.kernel_fft_seconds,
                       'one_structural_SD_matvec':matvec_seconds,'preconditioned_CG':cg_seconds,
                       'MP2_scalar_FFT_after_input_setup':mp2_seconds,'Neumann_sum':neumann_seconds},
            'retained_operator_data_bytes':m.output_memory_bytes(),
            'validation_scope':'Same finite effective model in double precision; residual arithmetic itself is not interval certified. No material, FT timing, or asymptotic advantage claim.'}
    assert abs(mp2-m.BDsq)<1e-12
    assert residual_certificate<tol/100
    for row in rows:
        assert row['absolute_error_vs_CG']<=row['uniform_bias_certificate']+1e-13
        if row['finite_absolute_row_bias_certificate']is not None:
            assert row['absolute_error_vs_CG']<=row['finite_absolute_row_bias_certificate']+1e-13
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--L',type=int,nargs='+',default=[2,4])
    parser.add_argument('--tol',type=float,default=1e-10)
    parser.add_argument('--output',default='benchmark_results.json')
    args=parser.parse_args()
    result={'status':'matched_classical_small_system_benchmark','numpy':np.__version__,'scipy':scipy.__version__,
            'python':platform.python_version(),'system':platform.system(),'machine':platform.machine(),
            'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS','not set'),
            'independent_full_fermion_validation':validate_small(),'requested_scalar_tolerance':args.tol,
            'instances':[run(L,args.tol)for L in args.L]}
    path=Path(__file__).parent / "results" / args.output
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
