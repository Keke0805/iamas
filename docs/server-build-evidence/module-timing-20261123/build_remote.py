#!/usr/bin/env python3
"""Build a new timer-only executable using the native GCC7.1/MPI toolchain."""
from pathlib import Path
import hashlib,json,os,shlex,subprocess
here=Path(__file__).resolve().parent;base=here.parent
ref=base/'precomm-20261118';vert=base/'vertmx-kokkos-20261122'
build=here/'build-sw64';build.mkdir(exist_ok=True);mods=build/'modules';mods.mkdir(exist_ok=True)
env=os.environ.copy();env['PATH']='/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr/bin:/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr_sw/bin:/usr/sw/mpi/mpi_20220608_SEA/bin:'+env['PATH'];env['LC_ALL']='en_US.utf8'
manifest=json.loads((here/'source-manifest.json').read_text())
for name,record in manifest.items():
 origin=(vert/'source'/name) if 'driver_drydep' in name else ref/'source-sw64'/name
 assert hashlib.sha256(origin.read_bytes()).hexdigest()==record['original_sha256'],str(origin)+' changed'
 assert hashlib.sha256((here/'source'/name).read_bytes()).hexdigest()==record['instrumented_sha256'],name
assert hashlib.sha256((vert/'atmosphere-native-vertmx-kokkos').read_bytes()).hexdigest()=='da809f71abf317ae8c81a7b3a907faa078f30161fa59df12797b43cdd34ff492'
flags=(ref/'build-sw64/src/atmosphere/CMakeFiles/core_atmosphere.dir/flags.make').read_text().splitlines()
vals={l.split(' = ',1)[0]:shlex.split(l.split(' = ',1)[1]) for l in flags if ' = ' in l}
f90='/usr/sw/mpi/mpi_20220608_SEA/bin/mpif90'
ff=[x for x in vals['Fortran_FLAGS'] if not x.startswith('-J')]
inc=['-I'+str(mods),'-I'+str(vert/'build-sw64/modules')]+vals['Fortran_INCLUDES']
order=['chemistry/schemes/module_mosaic_wetscav.F','chemistry/schemes/module_optical_averaging.F','chemistry/iamas_atmchem_todynamics.F','chemistry/iamas_atmchem_driver_wetscav.F','chemistry/iamas_atmchem_driver_emission.F','chemistry/iamas_atmchem_driver_drydep.F','dynamics/iamas_atm_time_integration.F','physics/iamas_atmphys_driver_microphysics.F','physics/iamas_atmphys_driver.F']
commands=[];objects=[]
for name in order:
 p=here/'source/src/atmosphere'/name;obj=build/(p.stem+'.o');objects.append(str(obj))
 args=[f90,'-mhost']+vals['Fortran_DEFINES']+inc+ff+['-cpp','-ffree-form','-J'+str(mods),'-c',str(p),'-o',str(obj)]
 commands.append(args)
 print('Compiling '+name,flush=True)
 with (build/(p.stem+'.log')).open('w') as f:
  proc=subprocess.run(args,cwd=str(build),env=env,stdout=f,stderr=subprocess.STDOUT)
 if proc.returncode:
  print((build/(p.stem+'.log')).read_text()[-16000:]);raise SystemExit(proc.returncode)
vert_objects=[str(vert/'build-sw64'/n) for n in ['vertmx_bridge.o','vertmx_kernel.o']]
exe=here/'atmosphere-native-module-timing'
cmd=shlex.split((ref/'native-link.command').read_text());cmd[1:1]=objects+vert_objects
cmd[cmd.index('-o')+1]=str(exe)
cmd=[x for x in cmd if not x.startswith('-Wl,-Map=')]+['-Wl,-Map='+str(build/'native-link.map')]
commands.append(cmd)
with (build/'link.log').open('w') as f:
 proc=subprocess.run(cmd,cwd=str(ref/'build-sw64/src/atmosphere'),env=env,stdout=f,stderr=subprocess.STDOUT)
if proc.returncode:print((build/'link.log').read_text()[-16000:]);raise SystemExit(proc.returncode)
(build/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
metadata={'sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'binary':str(exe),'reference_binary':str(vert/'atmosphere-native-vertmx-kokkos'),'new_objects':objects,'reused_vertmx_objects':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in vert_objects},'execution_space':'Existing Kokkos Serial on MPE','numerical_changes':False,'new_timer_labels':len(json.loads((here/'timer-catalog.json').read_text())),'fortran_flags':ff}
(here/'binary-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(metadata,indent=2),flush=True)
