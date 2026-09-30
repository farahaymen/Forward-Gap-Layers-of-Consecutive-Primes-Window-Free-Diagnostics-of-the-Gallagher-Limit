
/*
============================================================================
* fgp_sieve.c -- Forward-Gap Partition census of consecutive prime gaps.
*
* Bit-packed segmented sieve of Eratosthenes over odd integers. Streams the
* primes in increasing order and, for each prime p_n (n >= 2), accumulates the
* prime into the forward-gap layer L(g;N) indexed by its OUTGOING gap
* g = d_n = p_{n+1} - p_n. At each of a fixed schedule of checkpoints N_k the
* cumulative layer occupancies |L(g;N_k)| and layer prime-sums are snapshotted.
*
* Output (CSV, stdout):
* header line: k,N,X,g,count,sum
* one row per (checkpoint, gap) pair with nonzero count, where X = p_{N_k}.
*
* Build: gcc -O3 -march=native -funroll-loops -o fgp_sieve fgp_sieve.c -lm
* Run: ./fgp_sieve 10000000000 > fgp_layers.csv 2> fgp_sieve.log
*
* Memory: ~O(sqrt(LIMIT)) for the base primes plus one 256 KiB segment window,
* i.e. a few MiB in total, independent of LIMIT.
*
* Author: (manuscript author). Released under CC0 / public domain.
*
============================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
/* --- tunables ------------------------------------------------------ */
#define SEG_BYTES (1u << 18)   /* 256 KiB window: stays in L2 cache */
#define SEG_BITS  (SEG_BYTES * 8u) /* one bit per ODD integer */
#define SEG_SPAN  ((uint64_t)SEG_BITS * 2u)
#define GMAX      512           /* max gap tracked; the maximal gap below 10^10 is 354 */
#define MAXCK     64
static unsigned char *seg;      /* 1 = composite */
/* bit i of the segment stands for the odd integer base + 2*i */
#define SETBIT(a,i) ((a)[(i) >> 3] |= (unsigned char)(1u << ((i) & 7)))
#define GETBIT(a,i) ((a)[(i) >> 3] & (unsigned char)(1u << ((i) & 7)))
/* --- checkpoint schedule ------------------------------------------- */
static uint64_t ckN[MAXCK];
static int nck = 0;
/* snapshots */
static uint64_t snapCount[MAXCK][GMAX];
static uint64_t snapSum [MAXCK][GMAX];
static unsigned __int128 snapSqs[MAXCK][GMAX];
static uint64_t snapX [MAXCK];
/* running accumulators */
static uint64_t cnt[GMAX];
static uint64_t sum[GMAX];
static unsigned __int128 sqs[GMAX];   /* sum of p^2: needs 128 bits (p^2 ~ 1e20) */
/* print an unsigned __int128 in decimal */
static void fmt_u128(unsigned __int128 v, char *out)
{
    char tmp[40]; int i = 0;
    if (!v) { out[0] = '0'; out[1] = 0; return; }
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    for (int j = 0; j < i; j++) out[j] = tmp[i - 1 - j];
    out[i] = 0;
}
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}
/* Build a 20-point geometric checkpoint schedule in the prime INDEX N,
 * from N=2e5 up to N=pi(LIMIT), with the three scales quoted in the
 * original manuscript (10^6, 9x10^6, 10^8) forced into the grid. */
static void build_schedule(uint64_t piMax)
{
    const int K = 20;
    const double lo = 2.0e5;
    const double hi = (double)piMax;
    const double r = pow(hi / lo, 1.0 / (K - 1));
    double v = lo;
    for (int k = 0; k < K; k++) { ckN[k] = (uint64_t)(v + 0.5); v *= r; }
    ckN[K-1] = piMax;
    /* force-include the manuscript's scales by snapping the nearest grid point */
    const uint64_t forced[3] = {1000000ULL, 9000000ULL, 100000000ULL};
    for (int f = 0; f < 3; f++) {
        int best = 0; double bd = 1e300;
        for (int k = 0; k < K - 1; k++) {
            double d = fabs(log((double)ckN[k]) - log((double)forced[f]));
            if (d < bd) { bd = d; best = k; }
        }
        ckN[best] = forced[f];
    }
    qsort(ckN, K, sizeof(uint64_t), cmp_u64);
    /* force-snapping can collide two grid points; keep the schedule strictly
     * increasing so no checkpoint is emitted twice. */
    int m = 0;
    for (int k = 0; k < K; k++)
        if (k == 0 || ckN[k] != ckN[k-1]) ckN[m++] = ckN[k];
    if (m != K)
        fprintf(stderr, "note: %d duplicate checkpoint(s) removed\n", K - m);
    nck = m;
}
int main(int argc, char **argv)
{
    uint64_t LIMIT = (argc > 1) ? strtoull(argv[1], NULL, 10) : 10000000000ULL;
    /* pi(LIMIT) is not known in advance; the schedule needs an upper index.
     * We pass it explicitly as argv[2] when known, else use li(x) as estimate. */
    uint64_t piMax = (argc > 2) ? strtoull(argv[2], NULL, 10) : 0;
    if (!piMax) {
        double x = (double)LIMIT, lx = log(x);
        piMax = (uint64_t)(x / lx * (1 + 1/lx + 2/(lx*lx))); /* rough li(x) */
    }
    /* Optional checkpoint/restart: argv[3] = state file, argv[4] = time budget
     * in seconds of CPU time (clock(), not wall-clock). The run stops at a
     * segment boundary once the budget is exhausted, saves state, and exits with code 2; re-invoking with the same
     * arguments resumes exactly where it stopped. This makes very long runs
     * possible in environments with a per-process wall-clock limit, and does
     * not affect the result: segment boundaries are deterministic. */
    const char *statefile = (argc > 3) ? argv[3] : NULL;
    double budget = (argc > 4) ? atof(argv[4]) : 1e18;
    build_schedule(piMax);
    clock_t t0 = clock();
    /* --- base primes up to sqrt(LIMIT) ----------------------------- */
    uint64_t root = 0; while ((root + 1) * (root + 1) <= LIMIT) root++;
    unsigned char *small = calloc(root + 1, 1);
    uint64_t nbase = 0;
    for (uint64_t i = 3; i <= root; i += 2)
        if (!small[i]) { nbase++; for (uint64_t j = i * i; j <= root; j += 2 * i) small[j] = 1; }
    uint32_t *bp = malloc(nbase * sizeof(uint32_t));   /* base primes (odd) */
    uint64_t *next = malloc(nbase * sizeof(uint64_t)); /* next multiple */
    { uint64_t t = 0;
      for (uint64_t i = 3; i <= root; i += 2) if (!small[i]) bp[t++] = (uint32_t)i; }
    free(small);
    for (uint64_t i = 0; i < nbase; i++) next[i] = (uint64_t)bp[i] * bp[i];
    fprintf(stderr, "base primes below %llu: %llu\n",
            (unsigned long long)root, (unsigned long long)nbase);
    seg = malloc(SEG_BYTES);
    /* --- stream the primes ----------------------------------------- */
    uint64_t n = 1;              /* index of `prev` in the prime sequence */
    uint64_t prev = 2;           /* p_1 = 2 */
    uint64_t maxgap = 0, maxgap_at = 0;
    int ck = 0;
    uint64_t nprimes = 1;
    uint64_t start = 3;
    if (statefile) {  /* try to resume */
        FILE *f = fopen(statefile, "rb");
        if (f) {
            if (fread(&start,sizeof start,1,f) == 1) {
                int rd = 1;
                rd &= fread(&n,sizeof n,1,f) == 1; rd &= fread(&prev,sizeof prev,1,f) == 1;
                rd &= fread(&maxgap,sizeof maxgap,1,f) == 1; rd &= fread(&maxgap_at,sizeof maxgap_at,1,f) == 1;
                rd &= fread(&nprimes,sizeof nprimes,1,f) == 1; rd &= fread(&ck,sizeof ck,1,f) == 1;
                if (!rd) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(cnt,sizeof cnt,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(sum,sizeof sum,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(sqs,sizeof sqs,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(snapCount,sizeof snapCount,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(snapSum,sizeof snapSum,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(snapSqs,sizeof snapSqs,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                if (fread(snapX,sizeof snapX,1,f) != 1) { fprintf(stderr,"FATAL: truncated state file\n"); return 4; }
                fprintf(stderr,"resumed at %llu (n=%llu, ck=%d)\n",
                        (unsigned long long)start,(unsigned long long)n,ck);
            }
            fclose(f);
        }
    }
    for (uint64_t base = start; base <= LIMIT; base += SEG_SPAN) {
        if (statefile && (double)(clock()-t0)/CLOCKS_PER_SEC > budget) {
            FILE *f = fopen(statefile,"wb");
            fwrite(&base,sizeof base,1,f); fwrite(&n,sizeof n,1,f);
            fwrite(&prev,sizeof prev,1,f);
            fwrite(&maxgap,sizeof maxgap,1,f); fwrite(&maxgap_at,sizeof maxgap_at,1,f);
            fwrite(&nprimes,sizeof nprimes,1,f); fwrite(&ck,sizeof ck,1,f);
            fwrite(cnt,sizeof cnt,1,f); fwrite(sum,sizeof sum,1,f);
            fwrite(sqs,sizeof sqs,1,f);
            fwrite(snapCount,sizeof snapCount,1,f);
            fwrite(snapSum,sizeof snapSum,1,f);
            fwrite(snapSqs,sizeof snapSqs,1,f);
            fwrite(snapX,sizeof snapX,1,f);
            fclose(f);
            fprintf(stderr,"PAUSED at %llu (%.1f%% of range, n=%llu)\n",
                    (unsigned long long)base, 100.0*base/LIMIT,
                    (unsigned long long)n);
            return 2;
        }
        memset(seg, 0, SEG_BYTES);
        uint64_t hi = base + SEG_SPAN - 2;
        if (hi > LIMIT) hi = LIMIT;
        uint64_t nbits = (hi - base) / 2 + 1;
        for (uint64_t i = 0; i < nbase; i++) {
            uint64_t p = bp[i];
            if (p * p > hi) break;
            uint64_t s = next[i];
            if (s < base) {  /* first odd multiple >= base */
                s = base + (p - base % p) % p;
                if (!(s & 1)) s += p;
            }
            for (; s <= hi; s += 2 * p) SETBIT(seg, (s - base) >> 1);
            next[i] = s;
        }
        /* walk the survivors in order */
        uint64_t nbytes = (nbits + 7) >> 3;
        for (uint64_t b = 0; b < nbytes; b++) {
            unsigned char v = (unsigned char)~seg[b];
            while (v) {
                int t = __builtin_ctz((unsigned)v);
                v &= (unsigned char)(v - 1);
                uint64_t idx = b * 8 + (uint64_t)t;
                if (idx >= nbits) break;
                uint64_t q = base + 2 * idx;    /* q = p_{n+1} */
                uint64_t g = q - prev;
                if (n >= 2) {
                    if (g >= GMAX) {   /* silent truncation would break the partition */
                        fprintf(stderr, "FATAL: gap %llu at prime %llu exceeds GMAX=%d;"
                                " recompile with a larger GMAX\n",
                                (unsigned long long)g, (unsigned long long)prev, GMAX);
                        return 3;
                    }
                    cnt[g]++; sum[g] += prev;
                    sqs[g] += (unsigned __int128)prev * prev;
                }
                if (n >= 2 && g > maxgap) { maxgap = g; maxgap_at = prev; }
                /* snapshot once the layer contribution of index n is recorded */
                while (ck < nck && n == ckN[ck]) {
                    memcpy(snapCount[ck], cnt, sizeof(cnt));
                    memcpy(snapSum [ck], sum, sizeof(sum));
                    memcpy(snapSqs [ck], sqs, sizeof(sqs));
                    snapX[ck] = prev;   /* X = p_N */
                    fprintf(stderr, " checkpoint %2d: N=%llu X=p_N=%llu (%.1fs)\n",
                            ck, (unsigned long long)ckN[ck],
                            (unsigned long long)prev,
                            (double)(clock() - t0) / CLOCKS_PER_SEC);
                    ck++;
                }
                prev = q; n++; nprimes++;
            }
        }
    }
    fprintf(stderr, "pi(%llu) = %llu\nlargest prime <= LIMIT: %llu\n"
                    "max gap = %llu after %llu\nelapsed %.1fs\n",
                    (unsigned long long)LIMIT, (unsigned long long)nprimes,
                    (unsigned long long)prev, (unsigned long long)maxgap,
                    (unsigned long long)maxgap_at,
                    (double)(clock() - t0) / CLOCKS_PER_SEC);
    if (ck < nck)
        fprintf(stderr, "WARNING: only %d of %d checkpoints reached; "
                "pi(LIMIT) argument was probably too large\n", ck, nck);
    printf("k,N,X,g,count,sum,sumsq\n");
    char sq[40];
    for (int k = 0; k < ck; k++)
        for (int g = 2; g < GMAX; g += 2)
            if (snapCount[k][g]) {
                fmt_u128(snapSqs[k][g], sq);
                printf("%d,%llu,%llu,%d,%llu,%llu,%s\n", k,
                       (unsigned long long)ckN[k], (unsigned long long)snapX[k],
                       g, (unsigned long long)snapCount[k][g],
                       (unsigned long long)snapSum[k][g], sq);
            }
    return 0;
}
