"""Independent analytic and arbitrary-integer reference; no copied package code."""
import argparse, bisect, json, math, pathlib, shutil, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]

def near(a, b, tol=2e-10):
    if isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), (a, b)
        for x, y in zip(a, b): near(x, y, tol)
    else:
        assert math.isfinite(a) and abs(a-b) <= tol * max(1.0, abs(b)), (a,b)

def reference_ssa(seed):
    # Arbitrary integer modulo is independent of the package's Schrage method.
    def uniform():
        nonlocal seed
        seed = (16807 * seed) % 2147483647
        return seed / 2147483647
    state, t = [12,0], 0.0
    states, times, channels = [state[:]], [], []
    log_likelihood, score = 0.0, [0.0,0.0]
    while t < 4.0:
        a = [0.2 * math.comb(state[0], 2) * 2 / 3 if state[0]>=2 else 0.0, 0.7 * state[1]]
        total = math.fsum(a)
        dt = -math.log(uniform()) / total
        if dt > 4-t:
            log_likelihood -= total*(4-t)
            for j,k in enumerate([0.2,0.7]): score[j] -= (a[j]/k)*(4-t)
            break
        j = bisect.bisect_right([a[0],total], uniform()*total)
        log_likelihood += math.log(a[j])-total*dt
        for c,k in enumerate([0.2,0.7]): score[c] -= a[c]/k*dt
        score[j] += 1/[0.2,0.7][j]
        t += dt
        state = [state[0]-2,state[1]+1] if j==0 else [state[0]+2,state[1]-1]
        channels.append(j); times.append(t); states.append(state[:])
    return states,times,channels,log_likelihood,score

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--target',default='js'); args=parser.parse_args()
    moon=shutil.which('moon'); assert moon, 'moon not in PATH'
    proc=subprocess.run([moon,'run','cmd/reference','--target',args.target],cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf8',check=True)
    records=[json.loads(line) for line in proc.stdout.splitlines() if line.startswith('{')]
    kinds={}
    for v in records:
        kind=v['kind']; kinds[kind]=kinds.get(kind,0)+1
        if kind=='rate':
            r,s,n,k,V=v['input'],v['output'],v['counts'],v['rate'],v['volume']
            c=[x/V for x in n]; order=sum(r)
            # Combinations times factorial explicitly checks the convention.
            a=k*V**(1-order)*math.prod(math.comb(ni,ri)*math.factorial(ri) if ni>=ri else 0 for ni,ri in zip(n,r))
            f=k*math.prod(ci**ri for ci,ri in zip(c,r))
            delta=[si-ri for si,ri in zip(s,r)]
            J=[[di*k*r[j]*math.prod(c[l]**(r[l]-(1 if l==j else 0)) for l in range(2)) if r[j] else 0 for j in range(2)] for di in delta]
            near(v['propensity'],a); near(v['flux'],f); near(v['concentrations'],c)
            near(v['drift'],[x*f for x in delta]); near(v['jacobian'],J)
            near(v['diffusion'],[[di*dj*a for dj in delta] for di in delta])
        elif kind=='binomial':
            n,t=v['n'],v['time']; p=1/3+2/3*math.exp(-3*t)
            expected=[math.comb(n,i)*p**i*(1-p)**(n-i) for i in range(n+1)]
            near(v['masses'],expected,3e-12)
            assert 0<=v['tail']<=1e-12
            assert abs(math.fsum(v['masses'])-1)<3e-12
        elif kind=='ssa':
            states,times,channels,ll,score=reference_ssa(v['seed'])
            assert v['states']==states and v['path']['channels']==channels
            near(v['path']['times'],times); near(v['log_likelihood'],ll); near(v['score'],score)
        else: raise AssertionError(kind)
    assert kinds=={'rate':80,'binomial':25,'ssa':4}, kinds
    print('independent reference: 80 kinetic vectors, 25 full binomial distributions, 4 exact-seed event paths; all passed')
if __name__=='__main__': main()
