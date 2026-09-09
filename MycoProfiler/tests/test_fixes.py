"""Regression tests for the September fixes; external tools are stand-ins."""
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from test_mycoprofiler import TempCase, run_cli, ROOT

class ForceTests(TempCase):
    def test_failed_rerun_preserves_backup_without_stale_success(self):
        args = self.base_argv() + ['-q']
        self.assertEqual(run_cli(args), 0)
        self.assertEqual(run_cli(args + ['--force'], {'MOCK_KRAKEN2_FAIL':'BACTERIA'}), 1)
        self.assertFalse(Path(self.out, 'mycoprofiler_run_manifest.json').exists())
        self.assertFalse(Path(self.out, 'mycoprofiler_summary.tsv').exists())
        backups=list(Path(self.tmp).glob('out.backup-*'))
        self.assertEqual(len(backups), 1)
        self.assertTrue((backups[0]/'mycoprofiler_summary.tsv').exists())
    def test_sample_path_escape_rejected(self):
        self.assertEqual(run_cli(self.base_argv()+['--sample-id','../escape']),1)
        self.assertFalse(Path(self.out).exists())

class DownloadTests(unittest.TestCase):
    def test_interrupted_download_resumes_then_verified_skip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); bin=p/'bin';bin.mkdir()
            curl=bin/'curl'
            curl.write_text('''#!/usr/bin/env python3
import os,sys,pathlib
args=sys.argv[1:]; dest=pathlib.Path(args[args.index('-o')+1])
with open(os.environ['CALLS'],'a') as f:f.write('call\\n')
if os.environ.get('INTERRUPT')=='1':
 dest.write_bytes(b'123');sys.exit(18)
assert dest.read_bytes()==b'123'
dest.write_bytes(b'12345678')
''');curl.chmod(0o755)
            s=Path(ROOT,'scripts/download_bacterial_db.sh').read_text()
            funcs=s[s.index('file_bytes()'):s.index('for tool in')]
            harness=p/'harness.sh';harness.write_text('set -euo pipefail\n'+funcs+'\nfetch example "$DEST" 8\n')
            env=dict(os.environ,PATH=str(bin)+os.pathsep+os.environ['PATH'],DEST=str(p/'hash.k2d'),CALLS=str(p/'calls'))
            r=subprocess.run(['bash',str(harness)],env=dict(env,INTERRUPT='1'),capture_output=True)
            self.assertNotEqual(r.returncode,0);self.assertFalse((p/'hash.k2d').exists())
            self.assertEqual((p/'hash.k2d.part').read_bytes(),b'123')
            for _ in range(2):
                r=subprocess.run(['bash',str(harness)],env=env,capture_output=True)
                self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual((p/'calls').read_text().count('call'),2)
            self.assertEqual((p/'hash.k2d').read_bytes(),b'12345678')
    def test_legacy_partial_file_is_resumed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'hash').write_bytes(b'123')
            s=Path(ROOT,'scripts/download_bacterial_db.sh').read_text()
            funcs=s[s.index('file_bytes()'):s.index('for tool in')]
            harness=p/'test.sh';harness.write_text('set -euo pipefail\n'+funcs+'''
curl() { local dest=""; while [[ $# -gt 0 ]]; do if [[ "$1" == '-o' ]]; then dest="$2"; shift; fi; shift; done; [[ "$(cat "$dest")" == 123 ]]; printf 12345678 > "$dest"; }
fetch example "$1" 8
''')
            r=subprocess.run(['bash',str(harness),str(p/'hash')],capture_output=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual((p/'hash').read_bytes(),b'12345678')

class BuildTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.bin=self.p/'bin';self.bin.mkdir();self.cat=self.p/'cat';self.cat.mkdir()
        (self.cat/'A.fa').write_text('>contig\n'+'ACGT'*30+'\n')
        columns=['identifier','gtdb_domain','gtdb_phylum','gtdb_class','gtdb_order','gtdb_family','gtdb_genus','gtdb_species']
        with (self.p/'meta.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(columns);w.writerow(['A','d__Bacteria','p__P','c__C','o__O','f__F','g__G','s__S'])
        builder=self.bin/'kraken2-build';builder.write_text('''#!/usr/bin/env python3
import sys,os,pathlib,json
args=sys.argv[1:]
with open(os.environ['CALLS'],'a') as f:f.write(json.dumps(args)+'\\n')
if '--build' in args:
 p=pathlib.Path(args[args.index('--db')+1])
 for n in ['hash.k2d','opts.k2d','taxo.k2d']:(p/n).write_bytes(b'valid-mock')
''');builder.chmod(0o755)
        inspect=self.bin/'kraken2-inspect';inspect.write_text('''#!/usr/bin/env python3
import os,sys
if os.environ.get('BAD_INSPECT'):sys.exit(2)
print('100.00\\t42\\t0\\tR\\t1\\troot')
''');inspect.chmod(0o755)
        self.env=dict(os.environ,PATH=str(self.bin)+os.pathsep+os.environ['PATH'],CALLS=str(self.p/'calls'))
    def tearDown(self):self.tmp.cleanup()
    def run_build(self,*extra,**env):
        return subprocess.run([sys.executable,str(Path(ROOT,'scripts/build_smgc_kraken2_db.py')),'--catalogue',str(self.cat),'--metadata',str(self.p/'meta.csv'),'--output',str(self.p/'db'),*extra],env=dict(self.env,**env),capture_output=True)
    def test_masking_flag_and_provenance(self):
        import json
        r=self.run_build();self.assertEqual(r.returncode,0,r.stderr)
        cmds=[json.loads(x) for x in (self.p/'calls').read_text().splitlines()]
        self.assertIn('--no-masking',cmds[0]);self.assertNotIn('--no-masking',cmds[1])
        self.assertTrue((self.p/'db/build_provenance.json').exists())
    def test_masking_enabled(self):
        r=self.run_build('--mask-low-complexity');self.assertEqual(r.returncode,0,r.stderr)
        self.assertNotIn('--no-masking',(self.p/'calls').read_text())
    def test_missing_genome_stops_before_build(self):
        (self.cat/'A.fa').unlink();r=self.run_build()
        self.assertNotEqual(r.returncode,0);self.assertFalse((self.p/'calls').exists())
    def test_inspection_failure_is_fatal(self):
        r=self.run_build(BAD_INSPECT='1');self.assertNotEqual(r.returncode,0)
        self.assertFalse((self.p/'db/build_provenance.json').exists())
    def test_force_build_starts_fresh_and_backs_up(self):
        self.assertEqual(self.run_build().returncode,0)
        self.assertEqual(self.run_build('--force').returncode,0)
        self.assertEqual(len(list(self.p.glob('db.backup-*'))),1)

if __name__=='__main__': unittest.main()
