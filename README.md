# libgref

Libgref provides k-mer hashing index on string graphs on nucleotide sequences. It is implemented in pure C99 for Unix-like environments.

## Usage

The gref object takes three states regarding to the mutabilty and search capability. The first `gref_pool_t` object, which is generated on init, is a mutable sequence pool. Any sequences, links between sequences, and SNPs can be added to the pool in any order. The sequences are identified with its name, represented in an array of `char` and its length. The pool object can be converted to a `gref_acv_t` archive object, which can enumerate all the k-mers in the graph. The `gref_idx_t` object, which is converted from the archive object, provides k-mer matching (searching) functionality. The k-mers are represented in an array of `char` or a 2-bit-encoded sequence packed (in little-endian) in unsigned 64-bit variable.
``
```
/* initialize seqeunce pool */
gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
	.seed_length = 14));

/* append sequences and links */
gref_append_seq(pool, ...);
gref_append_link(pool, ...);

/* archive (build index on links) */
gref_acv_t *acv = gref_freeze_pool(pool);

/* enumerate kmers */
gref_iter_t *iter = gref_iter_init(acv);
uint64_t kmer;
while((kmer = gref_iter_next(acv)) != GREF_KMER_TERM) {
	/* do something on kmer */
}
gref_iter_clean(iter);

/* build index on kmers */
gref_idx_t *idx = gref_build_index(acv);

/* search kmer on graph */
struct gref_match_res_s match = gref_match(
	idx, "ACCTCCTTGCGCTT");

/* cleanup */
gref_clean(idx);
```

## Functions

### Init, destroy and conversion

#### gref\_init\_pool

Initialize empty pool.

```
gref_pool_t *gref_init_pool(
	gref_params_t const *params);
```

#### gref\_freeze\_pool

Build index on links to enable enumeration of kmers. (`pool` -> `acv` conversion)

```
gref_acv_t *gref_freeze_pool(
	gref_pool_t *pool)
```

#### gref\_melt\_archive

Destroy the index on links. (`acv` -> `pool` conversion)

```
gref_pool_t *gref_melt_archive(
	gref_acv_t *acv);
```

#### gref\_build\_index

Build index on kmers. (`acv` -> `idx` conversion)

```
gref_idx_t *gref_build_index(
	gref_acv_t *acv);
```

#### gref\_disable\_index

Destroy index on kmers. (`idx` -> `acv` conversion)

```
gref_acv_t *gref_disable_index(
	gref_idx_t *idx);
```

#### gref\_clean

Cleanup an object. `gref` can be any of `gref_pool_t`, `gref_acv_t`, or `gref_idx_t`.

```
void gref_clean(
	void *gref);
```

### Segment and link handling

#### gref\_append\_segment

Append a segment (a contiguous sequence) to `pool`. Returns 0 if succeeded.

```
int gref_append_segment(
	gref_pool_t *pool,
	char const *name,
	int32_t name_len,
	uint8_t const *seq,
	int64_t seq_len);
```

#### gref\_append\_link

Append a link to `pool`. Returns 0 if succeeded.

```
int gref_append_link(
	gref_pool_t *pool,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori);
```

#### gref\_append\_snp

TBD

#### gref\_split\_segment

TBD

### kmer enumeration

#### gref\_iter\_init

Initialize iterator.

```
gref_iter_t *gref_iter_init(
	gref_acv_t *acv);
```

#### gref\_iter\_next

Get the next kmer in the archive.

```
gref_kmer_tuple_t gref_iter_next(
	gref_iter_t *iter);
```

#### gref\_iter\_clean	

Cleanup.

```
void gref_iter_clean(
	gref_iter_t *iter);
```

### kmer search

#### gref\_match, gref\_match\_2bitpacked

Search kmer in the indexed graph.

```
struct gref_match_res_s gref_match(
	gref_idx_t const *idx,
	uint8_t const *seq);
struct gref_match_res_s gref_match_2bitpacked(
	gref_idx_t const *idx,
	uint64_t seq);
```

### Miscellaneous

#### gref\_get\_section

Get section info by section id.

```
struct gref_section_s const *gref_get_section(
	gref_idx_t const *ref,
	uint32_t id);
```

#### gref\_get\_name

Get section name by section id. 

```
struct gref_idx_str_s gref_get_name(
	gref_idx_t const *ref,
	uint32_t id);
```

#### gref\_get\_ptr

Returns a pointer to the sequence array.

```
uint8_t const *gref_get_ptr(
	gref_idx_t const *ref);
```

#### gref\_get\_total\_len

Total sequence length in the object.

```
int64_t gref_get_total_len(
	gref_idx_t const *ref);
```

## License

MIT
