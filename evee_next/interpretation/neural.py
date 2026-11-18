"""Exact arithmetic accounting for a frozen one-hidden-layer ReLU readout.

The active-region offset is input dependent. Terms reconstruct the calculation;
they are not effects of removing features or claims about biological mechanism.
"""
from dataclasses import dataclass,field
import hashlib,math,re
import numpy as np
from scipy.special import expit

INTERPRETATIONS=frozenset({'task_trained_DNA','biologically_unresolved_embedding','amino_acid_chemistry','availability_indicator','reference_structure_context','sequence_likelihood_not_mechanism','reference_context_x_chemistry'})
def require(ok,message):
 if not ok:raise ValueError(message)
def real_array(value):
 a=np.asarray(value);require(a.dtype.kind in 'fiu','Expected real numeric features/parameters');a=np.asarray(a,dtype=np.float64);require(np.isfinite(a).all(),'Nonfinite features/parameters');return a

def immutable(value):
 a=real_array(value);return np.frombuffer(np.ascontiguousarray(a).tobytes(),dtype=np.float64).reshape(a.shape)

@dataclass(frozen=True)
class FeatureBlock:
 name:str
 start:int
 end:int
 interpretation:str
 def __post_init__(self):
  require(isinstance(self.name,str) and re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,63}',self.name) is not None,'Invalid block name')
  require(type(self.start) is int and type(self.end) is int and 0<=self.start<self.end,'Invalid block range')
  require(self.interpretation in INTERPRETATIONS,'Unsupported biological interpretation')

@dataclass(frozen=True)
class AccountingTerm:
 block:str
 start:int
 end:int
 signed_logit:float
 interpretation:str

@dataclass(frozen=True)
class NeuralAccounting:
 value:float
 logit:float
 cutoff:float
 predicted_class:str
 output_bias:float
 active_hidden_bias:float
 activation_offset:float
 terms:tuple[AccountingTerm,...]
 active_units:tuple[int,...]
 reconstruction_error:float
 feature_map_id:str
 head_parameter_sha256:str
 method:str=field(default='active_relu_region_arithmetic',init=False)
 mechanism_status:str=field(default='unresolved',init=False)
 evaluation_scope:str=field(default='development_only',init=False)
 calibration_status:str=field(default='uncalibrated_sigmoid_score',init=False)
 terms_are_removal_effects:bool=field(default=False,init=False)

@dataclass(frozen=True)
class FrozenReluHead:
 W0:np.ndarray
 b0:np.ndarray
 W1:np.ndarray
 b1:float
 cutoff:float
 blocks:tuple[FeatureBlock,...]
 feature_map_id:str
 def __post_init__(self):
  for name in ['W0','b0','W1']:object.__setattr__(self,name,immutable(getattr(self,name)))
  for name in ['b1','cutoff']:
   scalar=real_array(getattr(self,name));require(scalar.shape==(),'Expected scalar output bias/cutoff');object.__setattr__(self,name,float(scalar))
  require(self.W0.ndim==2 and all(n>0 for n in self.W0.shape),'Invalid first-layer shape')
  require(self.b0.shape==self.W1.shape==(self.W0.shape[1],),'Invalid hidden/output shape')
  require(np.ndim(self.b1)==np.ndim(self.cutoff)==0 and math.isfinite(self.b1) and math.isfinite(self.cutoff) and 0<=self.cutoff<=1,'Invalid output bias/cutoff')
  require(isinstance(self.feature_map_id,str) and re.fullmatch(r'[0-9a-f]{64}',self.feature_map_id) is not None,'A qualified feature-map SHA is required')
  object.__setattr__(self,'blocks',tuple(self.blocks));require(self.blocks and all(isinstance(x,FeatureBlock) for x in self.blocks),'Invalid feature blocks')
  require(len({x.name for x in self.blocks})==len(self.blocks),'Repeated feature block')
  end=0
  for block in self.blocks:require(block.start==end,'Blocks must partition the full input in order');end=block.end
  require(end==self.W0.shape[0],'Blocks do not cover input coordinates')
 @property
 def parameter_sha256(self):
  h=hashlib.sha256()
  for a in [self.W0,self.b0,self.W1,np.asarray(self.b1,dtype=np.float64)]:h.update(str(a.shape).encode());h.update(np.ascontiguousarray(a).tobytes())
  return h.hexdigest()
 def account(self,features,feature_map_id):
  require(feature_map_id==self.feature_map_id,'Feature coordinate map differs from frozen head');x=real_array(features);require(x.shape==(self.W0.shape[0],),'Wrong transformed input shape')
  h=self.b0+x@self.W0;require(np.isfinite(h).all(),'Nonfinite neural preactivation');active=h>0;weighted=self.W1*active;gradient=self.W0@weighted;require(np.isfinite(gradient).all(),'Nonfinite effective coefficient')
  try:
   logit=float(self.b1+math.fsum(float(y*w) for y,w in zip(np.maximum(h,0),self.W1)));hidden=math.fsum(float(b*w) for b,w in zip(self.b0,weighted));offset=float(self.b1+hidden)
   terms=tuple(AccountingTerm(b.name,b.start,b.end,math.fsum(float(y*w) for y,w in zip(x[b.start:b.end],gradient[b.start:b.end])),b.interpretation) for b in self.blocks)
   reconstructed=offset+math.fsum(t.signed_logit for t in terms)
  except (OverflowError,FloatingPointError) as exc:raise ValueError('Nonfinite neural accounting') from exc
  require(all(math.isfinite(v) for v in [logit,hidden,offset,reconstructed]+[t.signed_logit for t in terms]),'Nonfinite neural accounting')
  error=abs(reconstructed-logit);require(error<=1e-10*max(1,abs(logit)),'Numerical accounting failed reconstruction');value=float(expit(logit))
  return NeuralAccounting(value,logit,float(self.cutoff),'pathogenic' if value>=self.cutoff else 'benign',float(self.b1),hidden,offset,terms,tuple(np.flatnonzero(active).tolist()),error,self.feature_map_id,self.parameter_sha256)
