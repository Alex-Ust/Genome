#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#define MAX_GENOME 6000000
#define MAX_LEN 128
#define K 16 // длина (seed)
#define ERR_PCT 8
#define MAX_ERR (MAX_LEN * ERR_PCT / 100)
#define MIN_LEN 30
#define MAX_HITS 1000
#define MAX_BEST 64
#define HASH_SIZE 4000037
#define STEP 64
#define BIN 10000
#define INF 1000000

#define UNIQUE 0
#define MULTI 1
#define UN_SHORT 2
#define UN_N 3
#define UN_NOSEED 4
#define UN_ERRORS 5
#define NSTATUS 6

const char *status_name[NSTATUS] = {
    "unique", "multi", "too short after trim", "too many N", "no seed found", "too many errors"
};
const char *letters = "ACGT";

typedef struct {
    int status;
    int strand; // 0 - strand, 1 - reverse
    int pos;
    int span;
    int dist;
    int len;
    int nops;
    char ops[2 * MAX_LEN]; // M fit, X change, I insert, D delete
} hit_t;

char genome[MAX_GENOME + 2];
int glen;

int bucket[HASH_SIZE + 1];
int pos[MAX_GENOME];
int fill_at[HASH_SIZE];

int sa[MAX_GENOME + 1];
char bwt[MAX_GENOME + 1];
int occ[MAX_GENOME / STEP + 2][4];
int cnt[4];

int use_bwt = 0;

int dp[MAX_LEN + 1][MAX_LEN + 2 * MAX_ERR + 1];
int found[MAX_HITS];
int cand[MAX_LEN / K * MAX_HITS];

long n_status[NSTATUS];
long n_trimmed;
long n_strand[2];
long n_dist[MAX_ERR + 1];
long bases[MAX_LEN], mism[MAX_LEN], ins[MAX_LEN], del[MAX_LEN];
int cov[MAX_GENOME + 1];

double seconds(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int code(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return -1;
}

char complement(char c) {
    if (c == 'A') return 'T';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    if (c == 'T') return 'A';
    return 'N';
}

void revcomp(const char *s, int n, char *out) {
    for (int i = 0; i < n; i++)
        out[i] = complement(s[n - 1 - i]);
    out[n] = '\0';
}

int cmp_int(const void *a, const void *b) {
    int x = *(int *)a;
    int y = *(int *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

void read_genome(const char *path) {
    FILE *f = fopen(path, "r");
    if (f == NULL) { printf(" cannot open %s \n", path); exit(1); }
    char line[4096];
    glen = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') continue;
        for (int i = 0; line[i] != '\0' && line[i] != '\n' && line[i] != '\r'; i++) {
            if (glen == MAX_GENOME) { printf("genome is too long \n"); exit(1); }
            genome[glen] = toupper(line[i]);
            glen++;
        }
    }
    fclose(f);
    genome[glen] = '$';
    genome[glen + 1] = '\0';
}


// ======================== HASH =========================

int kmer_bucket(const char *s) {
    unsigned long x = 0;
    for (int i = 0; i < K; i++) {
        int c = code(s[i]);
        if (c < 0) return -1;
        x = x * 4 + c;
    }
    return x % HASH_SIZE;
}

void build_hash(void) {
    int n = glen - K + 1;
    for (int b = 0; b <= HASH_SIZE; b++) bucket[b] = 0;

    for (int i = 0; i < n; i++) {
        int b = kmer_bucket(genome + i);
        if (b >= 0) bucket[b + 1]++;
    }
    for (int b = 0; b < HASH_SIZE; b++) bucket[b + 1] += bucket[b];

    for (int b = 0; b < HASH_SIZE; b++) fill_at[b] = bucket[b];
    for (int i = 0; i < n; i++) {
        int b = kmer_bucket(genome + i);
        if (b >= 0) {
            pos[fill_at[b]] = i;
            fill_at[b]++;
        }
    }
}

int find_hash(const char *s, int *out) {
    int b = kmer_bucket(s);
    if (b < 0) return 0;
    int n = 0;
    for (int i = bucket[b]; i < bucket[b + 1]; i++) {
        if (strncmp(genome + pos[i], s, K) != 0) continue;
        if (n == MAX_HITS) return -1;
        out[n] = pos[i];
        n++;
    }
    return n;
}


// ======================== BWT =========================

int cmp_suffix(const void *a, const void *b) {
    return strcmp(genome + *(int *)a, genome + *(int *)b);
}

void build_bwt(void) {
    int n = glen + 1;

    for (int i = 0; i < n; i++) sa[i] = i;
    qsort(sa, n, sizeof(int), cmp_suffix);

    for (int i = 0; i < n; i++) {
        if (sa[i] == 0) bwt[i] = '$';
        else bwt[i] = genome[sa[i] - 1];
    }

    int count[4] = {0, 0, 0, 0};
    for (int i = 0; i <= n; i++) {
        if (i % STEP == 0)
            for (int c = 0; c < 4; c++) occ[i / STEP][c] = count[c];
        if (i < n && code(bwt[i]) >= 0) count[code(bwt[i])]++;
    }
    cnt[0] = 1;
    for (int c = 1; c < 4; c++) cnt[c] = cnt[c - 1] + count[c - 1];
}

int occ_count(int c, int i) {
    int r = occ[i / STEP][c];
    for (int j = i / STEP * STEP; j < i; j++)
        if (bwt[j] == letters[c]) r++;
    return r;
}

int find_bwt(const char *s, int *out) {
    int lo = 0, hi = glen + 1;
    for (int i = K - 1; i >= 0; i--) {
        int c = code(s[i]);
        if (c < 0) return 0;
        lo = cnt[c] + occ_count(c, lo);
        hi = cnt[c] + occ_count(c, hi);
        if (lo >= hi) return 0;
    }
    if (hi - lo > MAX_HITS) return -1;
    for (int i = lo; i < hi; i++) out[i - lo] = sa[i];
    return hi - lo;
}

int find_seed(const char *s, int *out) {
    if (use_bwt) return find_bwt(s, out);
    return find_hash(s, out);
}

int align(const char *r, int n, const char *g, int e, hit_t *h) {
    for (int j = 0; j <= 2 * e; j++) dp[0][j] = 0;

    for (int i = 1; i <= n; i++) {
        int row_min = INF;
        for (int j = i; j <= i + 2 * e; j++) {
            int v = dp[i - 1][j - 1];
            if (r[i - 1] != g[j - 1]) v++;
            if (j < i + 2 * e && dp[i - 1][j] + 1 < v) v = dp[i - 1][j] + 1;
            if (j > i && dp[i][j - 1] + 1 < v) v = dp[i][j - 1] + 1;
            dp[i][j] = v;
            if (v < row_min) row_min = v;
        }
        if (row_min > e) return INF; // next wii be worse
    }

    int best = INF, end = 0;
    for (int j = n; j <= n + 2 * e; j++) {
        if (dp[n][j] <= best) {
            best = dp[n][j];
            end = j;
        }
    }
    if (best > e) return INF;

    char back[2 * MAX_LEN];
    int i = n, j = end, k = 0;
    while (i > 0) {
        int sub = dp[i - 1][j - 1];
        if (r[i - 1] != g[j - 1]) sub++;
        if (dp[i][j] == sub) {
            if (r[i - 1] == g[j - 1]) back[k] = 'M';
            else back[k] = 'X';
            i--;
            j--;
        } else if (j < i + 2 * e && dp[i][j] == dp[i - 1][j] + 1) {
            back[k] = 'I';
            i--;
        } else {
            back[k] = 'D';
            j--;
        }
        k++;
    }
    for (int t = 0; t < k; t++) h->ops[t] = back[k - 1 - t];
    h->nops = k;
    h->dist = best;
    h->pos = j;
    h->span = end - j;
    return best;
}


// One read mapping:

void map_read(const char *seq, const char *qual, int len, hit_t *h) {
    int n = len;
    while (n > 0 && qual[n - 1] <= '#') n--;
    h->len = n;
    if (n < MIN_LEN) {
        h->status = UN_SHORT;
        return;
    }

    int e = n * ERR_PCT / 100;
    int nn = 0;
    for (int i = 0; i < n; i++)
        if (seq[i] == 'N') nn++;
    if (nn > e) {
        h->status = UN_N;
        return;
    }

    char rc[MAX_LEN + 1];
    revcomp(seq, n, rc);

    int had_seed = 0, best = INF, nbest = 0;
    int best_pos[MAX_BEST], best_strand[MAX_BEST];

    for (int s = 0; s < 2; s++) {
        const char *r = seq;
        if (s == 1) r = rc;

        int nc = 0;
        for (int k = 0; k + K <= n; k += K) {
            int m = find_seed(r + k, found);
            if (m != 0) had_seed = 1;
            for (int i = 0; i < m; i++) {
                cand[nc] = found[i] - k;
                nc++;
            }
        }
        qsort(cand, nc, sizeof(int), cmp_int);

        int last = -INF;
        for (int i = 0; i < nc; i++) {
            if (cand[i] - last <= e) continue;
            last = cand[i];
            int from = cand[i] - e;
            if (from < 0) from = 0;
            if (from + n + 2 * e > glen) from = glen - n - 2 * e;

            hit_t a;
            int d = align(r, n, genome + from, e, &a);
            if (d == INF || d > best) continue;
            int p = from + a.pos;
            if (d < best) {
                best = d;
                nbest = 0;
            }

            int dup = 0;
            for (int b = 0; b < nbest && b < MAX_BEST; b++)
                if (best_strand[b] == s && abs(best_pos[b] - p) <= e) dup = 1;
            if (dup) continue;

            if (nbest < MAX_BEST) {
                best_pos[nbest] = p;
                best_strand[nbest] = s;
            }
            nbest++;

            if (rand() % nbest == 0) {
                *h = a;
                h->pos = p;
                h->strand = s;
                h->len = n;
            }
        }
    }

    if (nbest == 0 && had_seed) h->status = UN_ERRORS;
    else if (nbest == 0) h->status = UN_NOSEED;
    else if (nbest == 1) h->status = UNIQUE;
    else h->status = MULTI;
}

// ======================== STATS =========================

void add_stats(hit_t *h, int raw_len) {
    n_status[h->status]++;
    if (h->len < raw_len) n_trimmed++;
    if (h->status != UNIQUE && h->status != MULTI) return;

    n_strand[h->strand]++;
    n_dist[h->dist]++;

    int ri = 0;
    for (int k = 0; k < h->nops; k++) {
        int p = ri;
        if (h->strand == 1) p = h->len - 1 - ri;
        char op = h->ops[k];
        if (op == 'X') mism[p]++;
        if (op == 'I') ins[p]++;
        if (op == 'D') del[p]++;
        if (op != 'D') {
            bases[p]++;
            ri++;
        }
    }

    cov[h->pos]++;
    cov[h->pos + h->span]--;
}

double pct(long a, long b) {
    if (b == 0) return 0;
    return 100.0 * a / b;
}

void report(long total) {
    long mapped = n_status[UNIQUE] + n_status[MULTI];

    printf("\n");
    printf("Total reads: %ld\n", total);
    printf("Trimmed reads: %ld (%.2f%%)\n", n_trimmed, pct(n_trimmed, total));
    printf("Mapped reads: %ld (%.2f%%)\n", mapped, pct(mapped, total));
    printf("Unique: %ld (%.2f%%)\n", n_status[UNIQUE], pct(n_status[UNIQUE], total));
    printf("Multi-mapped: %ld (%.2f%%)\n", n_status[MULTI], pct(n_status[MULTI], total));
    printf("Unmapped reads: %ld (%.2f%%)\n", total - mapped, pct(total - mapped, total));

    printf("\nUnmapped because:\n");
    for (int s = UN_SHORT; s < NSTATUS; s++)
        printf("%s: %ld (%.2f%%)\n", status_name[s], n_status[s], pct(n_status[s], total));

    printf("\n");
    printf("Forward strand: %ld (%.2f%%)\n", n_strand[0], pct(n_strand[0], mapped));
    printf("Reverse strand: %ld (%.2f%%)\n", n_strand[1], pct(n_strand[1], mapped));

    printf("\nErrors in read, number of reads:\n");
    for (int d = 0; d <= MAX_ERR; d++)
        if (n_dist[d] > 0) printf("%d %ld (%.2f%%)\n", d, n_dist[d], pct(n_dist[d], mapped));

    long all = 0, x = 0, in = 0, de = 0;
    for (int i = 0; i < MAX_LEN; i++) {
        all += bases[i];
        x += mism[i];
        in += ins[i];
        de += del[i];
    }
    long wrong = x + in + de;
    printf("\n");
    printf("Identity: %.3f%%\n", 100 - pct(wrong, all));
    printf("Mismatches: %.1f%% of errors\n", pct(x, wrong));
    printf("Insertions: %.1f%% of errors\n", pct(in, wrong));
    printf("Deletions: %.1f%% of errors\n", pct(de, wrong));

    FILE *ef = fopen("errors.csv", "w");
    fprintf(ef, "pos,bases,mismatch,insertion,deletion\n");
    for (int i = 0; i < MAX_LEN; i++)
        if (bases[i] > 0) fprintf(ef, "%d,%ld,%ld,%ld,%ld\n", i + 1, bases[i], mism[i], ins[i], del[i]);
    fclose(ef);

    // +1/-1 -> глубина покрытия в каждой точке
    for (int i = 1; i < glen; i++) cov[i] += cov[i - 1];

    // сколько букв генома покрыто 0, 1, ..., 9 ридами и 10+
    long sum = 0;
    long by_depth[11] = {0};
    for (int i = 0; i < glen; i++) {
        sum += cov[i];
        if (cov[i] < 10) by_depth[cov[i]]++;
        else by_depth[10]++;
    }
    printf("\n");
    printf("Mean coverage: %.1f\n", (double)sum / glen);
    printf("\nCoverage, number of genome positions:\n");
    for (int d = 0; d < 10; d++)
        printf("%d %ld (%.3f%%)\n", d, by_depth[d], pct(by_depth[d], glen));
    printf("10+ %ld (%.3f%%)\n", by_depth[10], pct(by_depth[10], glen));
    printf("\nRegions with zero coverage:\n");

    int start = -1;
    for (int i = 0; i <= glen; i++) {
        int empty = i < glen && cov[i] == 0;
        if (empty && start == -1) start = i;
        if (!empty && start != -1) {
            printf("%d - %d, length %d\n", start + 1, i, i - start);
            start = -1;
        }
    }

    FILE *cf = fopen("coverage.csv", "w");
    fprintf(cf, "start,mean_depth\n");
    for (int b = 0; b + BIN <= glen; b += BIN) {
        long s = 0;
        for (int k = b; k < b + BIN; k++) s += cov[k];
        fprintf(cf, "%d,%.1f\n", b, (double)s / BIN);
    }
    fclose(cf);
}


int main(int argc, char **argv) {
    if (argc < 4 || (strcmp(argv[1], "hash") != 0 && strcmp(argv[1], "bwt") != 0)) {
        printf("usage: %s hash|bwt genome.fna reads.fastq [max_reads]\n", argv[0]);
        return 2;
    }
    if (strcmp(argv[1], "bwt") == 0) use_bwt = 1;
    long max_reads = -1;
    if (argc > 4) max_reads = atol(argv[4]);
    srand(1);

    double t0 = seconds();
    read_genome(argv[2]);
    long mem;
    if (use_bwt) {
        build_bwt();
        mem = (glen + 1) * (sizeof(int) + 1) + (glen / STEP + 1) * 4 * sizeof(int);
    } else {
        build_hash();
        mem = (HASH_SIZE + 1 + glen) * sizeof(int);
    }
    double t_index = seconds() - t0;

    FILE *f = fopen(argv[3], "r");
    if (f == NULL) { printf("cannot open %s\n", argv[3]); return 1; }

    char name[4096], seq[4096], plus[4096], qual[4096];
    long total = 0;
    double t1 = seconds();
    while (total != max_reads) {
        if (!fgets(name, sizeof(name), f)) break;
        if (!fgets(seq, sizeof(seq), f)) break;
        if (!fgets(plus, sizeof(plus), f)) break;
        if (!fgets(qual, sizeof(qual), f)) break;

        int len = 0;
        while (seq[len] != '\0' && seq[len] != '\n' && seq[len] != '\r') len++;
        if (len > MAX_LEN) len = MAX_LEN;

        hit_t h;
        map_read(seq, qual, len, &h);
        add_stats(&h, len);
        total++;
        if (total % 1000000 == 0) fprintf(stderr, "\r%ld reads", total);
    }
    double t_map = seconds() - t1;
    fprintf(stderr, "\n");
    fclose(f);

    printf("Index: %s (k = %d)\n", argv[1], K);
    printf("Genome length: %d\n", glen);
    printf("Index build time: %.2f sec\n", t_index);
    printf("Index memory: %.1f MB\n", mem / 1048576.0);
    printf("Mapping time: %.1f sec\n", t_map);
    report(total);

    return 0;

}
