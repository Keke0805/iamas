#!/usr/bin/env python3
"""Build an isolated SW64 experiment from an existing compatible source tree.
Run on 40.0.2.17 after copying this directory's deploy files to the server.
"""
from pathlib import Path
import argparse,os,subprocess,shlex,shutil,hashlib,json
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,required=True,help='Existing SW64-compatible source, including CMakeLists.txt')
args=parser.parse_args()
here=Path(__file__).resolve().parent
base=Path('/home/export/online1/mdt00/shisuan/swlicom1/hk/iamas')
src=here/'source-sw64';build=here/'build-sw64'
if not args.source.is_dir() or not (args.source/'CMakeLists.txt').is_file():raise SystemExit('Source tree missing')
if src.exists() or build.exists():raise SystemExit('Use a fresh experiment directory; source/build already exist')
# The original successful source contains toolchain compatibility edits. Keep all
# of them while changing only the transport storage and profiling table.
subprocess.run(['patch','--dry-run','-p1','-i',str(here/'precomm.patch')],cwd=args.source,check=True)
shutil.copytree(args.source,src,symlinks=True,ignore=shutil.ignore_patterns('.git','build'))
subprocess.run(['patch','-p1','-i',str(here/'precomm.patch')],cwd=src,check=True)
timer=src/'src/framework/iamas_timer.F'
if 'calls(call ranks' not in timer.read_text():
 subprocess.run(['patch','--dry-run','-p1','-i',str(here/'iamas-profile-timer.patch')],cwd=src,check=True)
 subprocess.run(['patch','-p1','-i',str(here/'iamas-profile-timer.patch')],cwd=src,check=True)
# Reuse the successful GCC7/SWGCC1449 configuration, with isolated source/build.
text=(here/'build-baseline-reference.sh').read_text()
text=text.replace('src=/tmp/iamas-gcc7-work/iamas-v2.6.5','src='+shlex.quote(str(src)))
text=text.replace('build=/tmp/iamas-gcc7-swgomp1449-opt-build','build='+shlex.quote(str(build)))
text='\n'.join(l for l in text.splitlines() if not l.startswith('python3 '))+'\n'
script=here/'build-sw64.sh';script.write_text(text)
with (here/'build-sw64.log').open('w') as log:
 subprocess.run(['bash',str(script)],cwd=here,stdout=log,stderr=subprocess.STDOUT,check=True)
cmd=shlex.split((build/'src/atmosphere/CMakeFiles/atmosphere.dir/link.txt').read_text())
if not cmd[0].endswith('/mpixfort.py'):raise SystemExit('Unexpected link wrapper')
cmd[0]='/usr/sw/mpi/mpi_20220608_SEA/bin/mpif90'
binary=here/'atmosphere-native-precomm'
cmd[cmd.index('-o')+1]=str(binary)
cmd.append('-Wl,-Map='+str(here/'native-link.map'))
env=os.environ.copy()
env['PATH']='/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr/bin:/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr_sw/bin:/usr/sw/mpi/mpi_20220608_SEA/bin:'+env['PATH']
env['LC_ALL']='en_US.utf8'
(here/'native-link.command').write_text(' '.join(shlex.quote(arg) for arg in cmd)+'\n')
with (here/'native-link.log').open('w') as log:
 subprocess.run(cmd,cwd=build/'src/atmosphere',env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
(here/'binary-metadata.json').write_text(json.dumps({'sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'compiler':'GCC7.1.0/SWGCC1449','optimization':'O2; same flags as baseline','vertical_padding':False,'source':str(src),'build':str(build)},indent=2)+'\n')
print('Built',binary)
