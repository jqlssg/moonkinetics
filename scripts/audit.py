#!/usr/bin/env python3
"""Source/history inventory; 1500 physical lines including ordinary whitespace."""
import pathlib, subprocess, re, json
ROOT = pathlib.Path(__file__).resolve().parents[1]
def git(*args):
    return subprocess.check_output(['git']+list(args),cwd=str(ROOT),encoding='utf8')
def main():
    sources=sorted(f for f in ROOT.glob('*.mbt') if not f.name.endswith(('_test.mbt','_wbtest.mbt')))
    rows=[{'file':f.name,'physical_lines':len(f.read_text(encoding='utf8').splitlines())} for f in sources]
    total=sum(r['physical_lines'] for r in rows)
    assert total>=1500, total
    commits=[]
    for sha in git('rev-list','--reverse','HEAD').splitlines():
        paths=git('diff-tree','--root','--no-commit-id','--name-only','-r',sha).splitlines()
        production=[p for p in paths if '/' not in p and p.endswith('.mbt') and not p.endswith(('_test.mbt','_wbtest.mbt'))]
        reference='scripts/reference.py' in paths and 'cmd/reference/main.mbt' in paths
        workflow=any(p in paths for p in ('cmd/biochemical/main.mbt','cmd/population/main.mbt','cmd/epidemic/main.mbt'))
        # These end-to-end integrations have numerical/invariant assertions, not print-only glue.
        # This is a candidate inventory; purpose/retained result must also be reviewed.
        if production or reference or workflow:
            commits.append({'sha':sha,'subject':git('show','-s','--format=%s',sha).strip(),'files':paths})
    assert len(commits)>=20, len(commits)
    readme=(ROOT/'README.md').read_text(encoding='utf8')
    quick=(ROOT/'cmd/quickstart/main.mbt').read_text(encoding='utf8').strip()
    assert quick in readme, 'README example drift'
    for target in re.findall(r'\]\(([^)]+)\)',readme):
        if '://' not in target and not target.startswith('#'):
            assert (ROOT/target.split('#')[0]).exists(), target
    assert 'jqlssg/moonkinetics@0.1.0' in readme
    print(json.dumps({'method':'Root handwritten non-test .mbt physical lines, normal whitespace included; cmd/generated/examples/vendor excluded', 'files':rows,'production_physical_lines':total,'engineering_commit_candidates':len(commits),'commits':commits},ensure_ascii=False,indent=2))
if __name__=='__main__': main()
