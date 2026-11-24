#!/usr/bin/env python3
"""Compile an isolated native-MPI binary, overriding just the drydep driver.
Run on SW64. Existing source, archives, executable and run data stay intact.
"""
from pathlib import Path
import os,subprocess,shlex,json,hashlib
here=Path(__file__).resolve().parent
base=Path('/home/export/online1/mdt00/shisuan/swlicom1/hk/iamas')
ref=base/'precomm-20261118'
build=here/'build-sw64';build.mkdir(exist_ok=True)
mods=build/'modules';mods.mkdir(exist_ok=True)
env=os.environ.copy();env['PATH']='/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr/bin:/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr_sw/bin:/usr/sw/mpi/mpi_20220608_SEA/bin:'+env['PATH'];env['LC_ALL']='en_US.utf8'
manifest=json.loads((here/'source-manifest.json').read_text())
name='src/atmosphere/chemistry/iamas_atmchem_driver_drydep.F'
if hashlib.sha256((ref/'source-sw64'/name).read_bytes()).hexdigest()!=manifest[name]['original_sha256']:
 raise SystemExit('Reference driver changed; inspect before compiling')
flags=(ref/'build-sw64/src/atmosphere/CMakeFiles/core_atmosphere.dir/flags.make').read_text().splitlines()
vals={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags if ' = ' in l}
f90='/usr/sw/mpi/mpi_20220608_SEA/bin/mpif90';cxx='/usr/sw/mpi/mpi_20220608_SEA/bin/mpicxx'
ff=[x for x in vals['Fortran_FLAGS'] if not x.startswith('-J')]
inc=['-I'+str(mods)]+vals['Fortran_INCLUDES']
schemes=here/'source/src/atmosphere/chemistry/schemes'
commands=[]
def run(args,name):
 commands.append(args)
 with (build/(name+'.log')).open('w') as f:
  p=subprocess.run(args,cwd=build,env=env,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:
  print((build/(name+'.log')).read_text()[-12000:]);raise SystemExit(p.returncode)
run([f90,'-mhost']+vals['Fortran_DEFINES']+inc+ff+['-cpp','-ffree-form','-J'+str(mods),'-c',str(schemes/'module_vertmx_kokkos.F'),'-o',str(build/'vertmx_bridge.o')],'bridge')
run([f90,'-mhost']+vals['Fortran_DEFINES']+inc+ff+['-cpp','-ffree-form','-J'+str(mods),'-c',str(here/'source'/name),'-o',str(build/'drydep.o')],'drydep')
run([cxx,'-mhost','-O2','-g','-std=c++1z','-mieee','-mftz','-msimd','-fopenmp','-ffp-contract=off','-I/tmp/iamas-gcc7-deps/flcl/include','-I/tmp/iamas-gcc7-deps/kokkos/include','-c',str(schemes/'module_vertmx_kokkos.cpp'),'-o',str(build/'vertmx_kernel.o')],'kernel')
cmd=shlex.split((ref/'native-link.command').read_text())
cmd[1:1]=[str(build/n) for n in ['drydep.o','vertmx_bridge.o','vertmx_kernel.o']]
exe=here/'atmosphere-native-vertmx-kokkos'
cmd[cmd.index('-o')+1]=str(exe)
cmd=[x for x in cmd if not x.startswith('-Wl,-Map=')]+['-Wl,-Map='+str(build/'native-link.map')]
commands.append(cmd)
with (build/'link.log').open('w') as f:
 p=subprocess.run(cmd,cwd=ref/'build-sw64/src/atmosphere',env=env,stdout=f,stderr=subprocess.STDOUT)
if p.returncode:print((build/'link.log').read_text()[-12000:]);raise SystemExit(p.returncode)
(build/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(here/'binary-metadata.json').write_text(json.dumps({'sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'reference_binary':str(ref/'atmosphere-native-precomm'),'new_objects':['drydep.o','vertmx_bridge.o','vertmx_kernel.o'],'execution_space':'Kokkos DefaultHostExecutionSpace, existing SERIAL library','modes':{'0':'original Fortran with parent timer','1':'Kokkos batch','2':'Kokkos plus full local field comparison against Fortran'}},indent=2)+'\n')
print('Built',exe,flush=True)
