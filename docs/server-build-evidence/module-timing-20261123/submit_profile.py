#!/usr/bin/env python3
"""Submit isolated instrumentation smoke test or an unchanged formal 12h case."""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,sys,time,xml.etree.ElementTree as ET
root=Path(__file__).resolve().parent;base=root.parent
kind=sys.argv[1] if len(sys.argv)>1 else 'smoke'
assert kind in ('smoke','12h')
ref=base/'vertmx-kokkos-20261122/run-vertmx-kokkos-12h-20261122-233310'
exe=root/'atmosphere-native-module-timing'
build=json.loads((root/'binary-metadata.json').read_text(encoding='utf-8'))
assert hashlib.sha256(exe.read_bytes()).hexdigest()==build['sha256']
run=root/('run-'+kind+'-'+time.strftime('%Y%m%d-%H%M%S'));run.mkdir();(run/'log').mkdir()
for p in (base/'data/test_iAMAS').iterdir():
 if p.name in ('atmosphere','log','namelist.atmosphere','streams.atmosphere'):continue
 (run/p.name).symlink_to(p)
for name in ('namelist.atmosphere','streams.atmosphere'):shutil.copy2(str(ref/name),str(run/name))
if kind=='smoke':
 p=run/'namelist.atmosphere';text=p.read_text(encoding='utf-8')
 text,n=re.subn(r'(?m)^(\s*config_run_duration\s*=\s*).*$',lambda m:m[1]+'0_00:06:00',text);assert n==1;p.write_text(text, encoding='utf-8')
 p=run/'streams.atmosphere';tree=ET.fromstring(p.read_text(encoding='utf-8'))
 for stream in tree:
  if stream.get('name') in ('output','diagnostics','restart'):stream.set('output_interval','none')
 p.write_text(ET.tostring(tree,encoding='unicode'), encoding='utf-8')
else:
 for name in ('namelist.atmosphere','streams.atmosphere'):assert (run/name).read_bytes()==(ref/name).read_bytes()
shutil.copy2(str(exe),str(run/'atmosphere'))
env=os.environ.copy();env.update({'IAMAS_PRECOMM_MODE':'3','IAMAS_VERTMX_MODE':'1','OMP_NUM_THREADS':'1','SWCH_USE_HW_COLL':'0','LC_ALL':'en_US.utf8'});env.pop('IAMAS_MONO_KPAD',None)
meta={'kind':kind,'run_directory':str(run),'steps':4 if kind=='smoke' else 480,'simulation_minutes':6 if kind=='smoke' else 720,'mpe':360,'dt':90,'pio_tasks':6,'pio_stride':60,'binary_sha256':build['sha256'],'reference_run':str(ref),'precomm_mode':3,'vertmx_mode':1,'vertical_padding':'none','new_timer_labels':build['new_timer_labels'],'formal_output':kind=='12h','configuration_hashes':{n:hashlib.sha256((run/n).read_bytes()).hexdigest() for n in ('namelist.atmosphere','streams.atmosphere')},'environment':{k:env.get(k,'') for k in ('IAMAS_PRECOMM_MODE','IAMAS_VERTMX_MODE','SWCH_USE_HW_COLL','OMP_NUM_THREADS')},'status':'prepared','submitted_at':time.strftime('%Y-%m-%d %H:%M:%S %z')}
cmd=['bsub','-b','-pr','-K','-timelimit','0-01:00:00' if kind=='smoke' else '0-06:00:00','-q','q_swopt','-J','module_timing_'+kind,'-o','./log/run360.out','-n','360','-cgsp','64','-host_stack','1024','-share_size','13312','-cache_size','128','./atmosphere','-n','namelist.atmosphere','-s','streams.atmosphere']
p=subprocess.run(cmd,cwd=str(run),env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
(run/'submit.log').write_text(p.stdout, encoding='utf-8')
m=re.search(r'Job <(\d+)>',p.stdout)
if p.returncode or not m:raise SystemExit(p.stdout)
meta.update({'job_id':int(m[1]),'submission_command':cmd,'status':'submitted'})
(run/'run-metadata.json').write_text(json.dumps(meta,indent=2)+'\n', encoding='utf-8')
(root/('latest-'+kind+'-run.txt')).write_text(str(run)+'\n', encoding='utf-8')
with (run/'log/collector.log').open('a') as f:
 watcher=subprocess.Popen(['python3','-u',str(root/'watch_profile.py'),str(run)],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,env=env,start_new_session=True)
(run/'collector.pid').write_text(str(watcher.pid)+'\n', encoding='utf-8')
print(json.dumps({'job_id':meta['job_id'],'run':str(run),'watcher_pid':watcher.pid},indent=2),flush=True)
