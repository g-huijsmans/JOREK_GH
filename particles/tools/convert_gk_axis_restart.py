"""Explicit conversion of an initial model011 field-import file, never a population restart.
Geometry and all physical Hermite fields are preserved. Only equation topology,
treat_axis and the equation-count metadata change. Requires zero potential.
"""
from pathlib import Path
import argparse,json,hashlib,shutil
import h5py,numpy as np

def convert(src,dst):
 src=Path(src).resolve();dst=Path(dst).resolve()
 if src==dst or dst.exists():raise ValueError('Output must be a new file, distinct from the input.')
 with h5py.File(src,'r')as f:
  if int(f['n_order'][0])!=3:raise ValueError('Only cubic Hermite field imports are supported.')
  if f['treat_axis'][0]!=b'F':raise ValueError('Input is already marked treated; no implicit reconversion.')
  ax=np.flatnonzero(f['axis_node'][:]==b'T');idx=f['index'][:].copy();x=f['x'][:];oldmax=int(idx.max());axis_dof=f['axis_dof'][:][ax]
  if len(ax)<2 or not np.all(np.isin(axis_dof,[2,3])):raise ValueError('Unsupported/missing axis metadata.')
  if np.any(f['constrained'][:][ax]==b'T'):raise ValueError('Constrained axis nodes are unsupported.')
  if np.any(idx<=0)or not np.array_equal(np.unique(idx),np.arange(1,oldmax+1)):raise ValueError('Input must have positive dense equation numbering.')
  coords=x[:,:,0,:].transpose(2,1,0)
  if not np.all(coords[ax,0]==coords[ax[0],0]):raise ValueError('Axis nodes do not coincide exactly.')
  vertices=f['vertex'][:].T-1
  for v in vertices[np.isin(vertices,ax).any(axis=1)]:
   flags=np.isin(v,ax)
   if not np.array_equal(flags,[True,False,False,True]):raise ValueError('Unsupported axis vertex ordering.')
  if np.any(f['values'][1]!=0):raise ValueError('Only zero-potential initial field imports are supported; no evolved state clipping.')
  # Build equivalence classes as the grid generator does, preserving other shared equations.
  representative=np.arange(oldmax+1)
  for dof in range(4):representative[idx[dof,ax]]=idx[dof,ax[0]]
  equivalence=representative[idx];labels=np.unique(equivalence);new=np.searchsorted(labels,equivalence)+1
  assert len(labels)==oldmax-4*(len(ax)-1)
  assert all(len(np.unique(new[d,ax]))==1 for d in range(4))
  eq=int(len(labels)*int(f['n_tor'][0])*int(f['n_var'][0]))
  metadata={'source':str(src),'output':str(dst),'old_spatial_equations':oldmax,'new_spatial_equations':len(labels),'axis_nodes':len(ax),'axis_dof':np.unique(axis_dof).tolist(),'new_file_equation_count':eq,'model011_field_import_recomputes_n_dof':True,'qualification':'Initial field import only; all physical fields retained, no particle/operator restart.'}
 dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
 with h5py.File(dst,'r+')as f:f['index'][...]=new;f['treat_axis'][...]=b'T';f['n_dof'][...]=eq
 changed=[];checked=0
 with h5py.File(src)as a,h5py.File(dst)as b:
  for k in a:
   if not isinstance(a[k],h5py.Dataset):continue
   checked+=1
   same = np.array_equal(a[k][...], b[k][...], equal_nan=True) if a[k].dtype.kind in 'fc' else np.array_equal(a[k][...], b[k][...])
   if not same: changed.append(k)
 assert set(changed)=={'index','treat_axis','n_dof'},changed
 metadata.update(changed_datasets=changed,checked_datasets=checked,geometry_connectivity_and_physical_fields_identical=True)
 dst.with_suffix('.axis_conversion.json').write_text(json.dumps(metadata,indent=2));print(json.dumps(metadata,indent=2));return metadata
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output');a=p.parse_args();convert(a.source,a.output)
