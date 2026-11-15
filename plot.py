import matplotlib.pyplot as plt

x = []
depth = []
for line in open("coverage.csv").readlines()[1:]:
    start, d = line.split(",")
    x.append(int(start) / 1000000)
    depth.append(float(d))

mean = sum(depth) / len(depth)
plt.figure(figsize=(12, 4))
plt.plot(x, depth)
plt.axhline(mean, color="gray", linestyle="--", label="mean " + str(round(mean)) + "x")
plt.xlabel("position in genome, Mb")
plt.ylabel("mean coverage (10 kb windows)")
plt.title("Genome coverage")
plt.legend()
plt.savefig("coverage.png")

pos = []
mismatch = []
insertion = []
deletion = []
for line in open("errors.csv").readlines()[1:]:
    p, bases, m, i, d = line.split(",")
    pos.append(int(p))
    mismatch.append(100 * int(m) / int(bases))
    insertion.append(100 * int(i) / int(bases))
    deletion.append(100 * int(d) / int(bases))

plt.figure(figsize=(10, 4))
plt.plot(pos, mismatch, label="mismatch")
plt.plot(pos, insertion, label="insertion")
plt.plot(pos, deletion, label="deletion")
plt.xlabel("position in read")
plt.ylabel("errors, % of bases")
plt.title("Errors by position in read")
plt.legend()
plt.savefig("errors.png")
