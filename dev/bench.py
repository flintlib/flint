import subprocess

def run(s):
    print("running", s)
    s = subprocess.check_output(s.split()).decode()
    for line in s.splitlines():
        print("--- ", line)
    return s

results = []

def add_result(prefix, time):
    global results
    results.append((prefix, time))

def print_results():
    global results
    print()
    print("==================== DETAILS ====================")
    total = 0.0
    groups = {}
    group_count = {}
    for prefix, time in results:
        total += time
        descr = " ".join(prefix)
        print(descr.ljust(40), f'{time:6.3f}')
        group = prefix[0]
        if group in groups:
            groups[group] += time
        else:
            groups[group] = time
        group_count[group] = group_count.get(group, 0) + 1
    print("==================== BY GROUP ===================")
    for group in groups:
        descr = "%s (%i benchmarks)" % (group, group_count[group])
        time = groups[group]
        print(descr.ljust(40), f'{time:6.3f}')
    print("====================== TOTAL ====================")
    descr = "all (%s benchmarks)" % len(results)
    print(descr.ljust(40), f'{total:6.3f}')
    print()

def arb_bench():
    s = run("build/examples/pi 10000000")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "pi 10000000"], time)

    s = run("build/examples/pi 10000000 -threads 8")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "pi 10000000 -threads 8"], time)

    s = run("build/examples/bernoulli 1000000 -threads 8")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "bernoulli 1000000 -threads 8"], time)

    s = run("build/examples/partitions 1000000000000 -threads 8 -quiet")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "partitions 1e12 -threads 8"], time)

    s = run("build/examples/zeta_zeros -count 1000 -prec 333 -threads 8")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "zeta_zeros 1000 -prec 333 -threads 8"], time)

    s = run("build/examples/keiper_li 3000 -threads 8")
    for line in s.splitlines():
        if "total: " in line:
            time = float(line.split()[-1])
    add_result(["arb", "keiper_li 3000 -threads 8"], time)

    s = run("build/examples/integrals -i 4 -prec 1000")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "integrals -i 4 -prec 1000"], time)

    s = run("build/examples/hilbert_matrix 120")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "hilbert_mat 120"], time)

    s = run("build/examples/hilbert_matrix -eig 50")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "hilbert_matrix -eig 50"], time)

    s = run("build/examples/poly_roots a 256")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "poly_roots a 256"], time)

    s = run("build/examples/class_poly -1000000")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            time = float(line.split()[-1])
    add_result(["arb", "class_poly -1000000"], time)

def fmpz_bench():
    s = run("build/fmpz/profile/p-is_prime_bench")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            size = " ".join(line.split()[:3])
            time = float(line.split()[-1])
            add_result(["fmpz", "is_prime", size], time)

    s = run("build/fmpz/profile/p-is_probabprime_bench")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            size = " ".join(line.split()[:3])
            time = float(line.split()[-1])
            add_result(["fmpz", "is_probabprime", size], time)

    s = run("build/fmpz_factor/profile/p-factor_bench")
    for line in s.splitlines():
        if "cpu/wall(s):" in line:
            size = " ".join(line.split()[:3])
            time = float(line.split()[-1])
            add_result(["fmpz", "factor", size], time)



def fmpz_poly_bench():
    s = run("build/fmpz_poly_factor/profile/p-factor_hard")
    for line in s.splitlines():
        if " ms" in line:
            poly = line.split()[0]
            time = float(line.split()[-2]) * 0.001
            add_result(["fmpz_poly", "factor", poly], time)

def fmpz_mpoly_bench():
    s = run("build/fmpz_mpoly_factor/profile/p-factor")
    for line in s.splitlines():
        if "variables" in line:
            vars = line.split()[1]
        if "total_time:" in line:
            time = float(line.split()[-1]) * 0.001
            add_result(["fmpz_mpoly", "factor %s variables" % vars], time)

def nmod_mpoly_bench():
    s = run("build/nmod_mpoly_factor/profile/p-factor")
    for line in s.splitlines():
        if "---" in line:
            info = " ".join(line.split()[1:-1])
            info = info.replace("variables", "vars").replace("characteristic", "p =")
        if "total_time:" in line:
            time = float(line.split()[-1]) * 0.001
            add_result(["nmod_mpoly", "factor %s" % info], time)

def calcium_bench():
    s = run("build/examples/huge_expr -ca")
    for line in s.splitlines():
        if "Total: cpu/wall(s):" in line:
            time = float(line.split()[-1])
            add_result(["calcium", "huge_expr"], time)

fmpz_bench()
fmpz_poly_bench()
nmod_mpoly_bench()
fmpz_mpoly_bench()
arb_bench()
calcium_bench()

print_results()

