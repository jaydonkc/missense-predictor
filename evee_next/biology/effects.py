"""Pure continuous-effect fitting and paired biological discrimination metrics.

No implicit data access, clinical labels, thresholds, or external model calls.
"""
from collections import defaultdict
import numpy as np
from scipy.stats import rankdata
from sklearn.linear_model import Ridge


def defined_ge(a,b):
    return a is not None and b is not None and bool(np.isfinite(a) and np.isfinite(b) and a>=b)


def finite(values, ndim=None):
    a=np.asarray(values,dtype=np.float64)
    if (ndim is not None and a.ndim!=ndim) or not np.isfinite(a).all():
        raise ValueError('Finite numeric array with the expected shape required')
    return a


def transfer_folds(genes):
    g=np.asarray(genes)
    if g.ndim!=1 or set(g)!= {'ASPA','G6PD'}:
        raise ValueError('Exactly the currently permitted ASPA/G6PD sources required')
    return [(str(gene),g!=gene,g==gene) for gene in ('ASPA','G6PD')]


def context_donors(rows):
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        groups[(r['gene'],r['ref_aa'],r['alt_aa'])].append(i)
    donor=np.arange(len(rows));informative=np.zeros(len(rows),dtype=bool)
    for indices in groups.values():
        indices=sorted(indices,key=lambda i:(rows[i]['position0'],rows[i]['id']))
        if len({rows[i]['position0'] for i in indices})!=len(indices):
            raise ValueError('Duplicate source/site/substitution identity')
        if len(indices)>1:
            donor[indices]=np.roll(indices,1);informative[indices]=True
    return donor,informative


def fit_ridge(x,y,alpha=10.):
    x=finite(x,2);y=finite(y,1)
    if len(x)!=len(y) or len(y)<2 or x.shape[1]<1 or not np.isfinite(alpha) or alpha<=0:
        raise ValueError('Invalid fitting dimensions or regularization')
    mean=x.mean(axis=0);scale=x.std(axis=0);scale[scale==0]=1.
    z=(x-mean)/scale
    model=Ridge(alpha=alpha,fit_intercept=True,solver='cholesky').fit(z,y)
    return dict(mean=mean.tolist(),scale=scale.tolist(),coef=model.coef_.tolist(),
                intercept=float(model.intercept_),alpha=float(alpha),fitting_rows=len(y))


def predict_ridge(state,x):
    x=finite(x,2);m=finite(state['mean'],1);s=finite(state['scale'],1);c=finite(state['coef'],1)
    if m.shape!=s.shape or m.shape!=c.shape or x.shape[1]!=len(m) or (s<=0).any():
        raise ValueError('Model feature dimensions/scales differ')
    return finite(((x-m)/s)@c+state['intercept'],1)


def metrics(y,p):
    y=finite(y,1);p=finite(p,1)
    if y.shape!=p.shape or len(y)<2: raise ValueError('Paired observations required')
    def corr(a,b):
        a=a-a.mean();b=b-b.mean();den=np.linalg.norm(a)*np.linalg.norm(b)
        return float(a@b/den) if den>0 else None
    mse=float(np.mean((y-p)**2));var=float(np.var(y))
    return dict(n=len(y),spearman=corr(rankdata(y),rankdata(p)),pearson=corr(y,p),
                rmse=mse**.5,mae=float(np.mean(abs(y-p))),r2=1-mse/var if var>0 else None)


class Pairs:
    """Precomputed within-group pairs; target ties excluded from concordance."""
    def __init__(self,y,groups,blocks=None):
        y=finite(y,1);g=np.asarray(groups);self.n=len(y)
        if g.shape!=y.shape: raise ValueError('Group identity count differs')
        left=[];right=[]
        for label in sorted(set(g.tolist())):
            ix=np.flatnonzero(g==label);a,b=np.triu_indices(len(ix),1)
            left.extend(ix[a]);right.extend(ix[b])
        self.left=np.asarray(left,dtype=int);self.right=np.asarray(right,dtype=int)
        # Stable dense group IDs for pair aggregation, independent of source row order.
        _,codes=np.unique(g,return_inverse=True);self.group=codes[self.left]
        self.group_count=int(codes.max()+1) if len(codes) else 0
        self.target_sign=np.sign(y[self.right]-y[self.left]);self.strict=self.target_sign!=0
        self.blocks=None if blocks is None else np.asarray(blocks,dtype=int)
        if self.blocks is not None and (self.blocks.shape!=y.shape or (self.blocks<0).any()):
            raise ValueError('Invalid block identities')

    def scores(self,p):
        p=finite(p,1)
        if p.shape!=(self.n,):raise ValueError('Prediction row count differs')
        sign=np.sign(p[self.right]-p[self.left])
        return np.where(sign==0,.5,(sign==self.target_sign).astype(float))

    def summarize(self,p,multiplicity=None):
        s=self.scores(p)
        if multiplicity is None:w=np.ones(len(s))
        else:
            m=finite(multiplicity,1)
            if self.blocks is None or len(m)<=self.blocks.max() or (m<0).any() or (m!=np.floor(m)).any():
                raise ValueError('Integer block multiplicities required')
            a=self.blocks[self.left];b=self.blocks[self.right]
            w=m[a]*m[b]
        keep=self.strict
        denominator=np.bincount(self.group[keep],weights=w[keep],minlength=self.group_count)
        numerator=np.bincount(self.group[keep],weights=(s*w)[keep],minlength=self.group_count)
        valid=denominator>0
        return dict(equal_group_concordance=float(np.mean(numerator[valid]/denominator[valid])) if valid.any() else None,
                    informative_groups=int(valid.sum()),strict_pairs=int(w[keep].sum()),
                    target_ties=int(w[~keep].sum()),prediction_ties=int(w[keep&(s==.5)].sum()),
                    correct_pairs=int(w[keep&(s==1)].sum()),wrong_pairs=int(w[keep&(s==0)].sum()))


def concordance(y,p,groups,blocks=None,multiplicity=None):
    return Pairs(y,groups,blocks).summarize(p,multiplicity)


def designs(chemistry,context,protein,rows):
    c=finite(chemistry,2);r=finite(context,2);p=finite(protein,2)
    if len(c)!=len(r) or len(c)!=len(p) or len(rows)!=len(c) or r.shape[1]!=3:
        raise ValueError('Design rows or reference context shape differ')
    donor,informative=context_donors(rows)
    # Directed physical properties and absolute changes; no one-hot expansion here.
    def interaction(z):return (c[:,:10,None]*z[:,None,:2]).reshape(len(c),20)
    base=np.column_stack([c,p]);rd=r[donor]
    return dict(C=c,R=r,K=np.column_stack([c,r,interaction(r)]),P=base,
                J=np.column_stack([base,r,interaction(r)]),
                D=np.column_stack([base,rd,interaction(rd)])),donor,informative
