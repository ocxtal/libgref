
/**
 * @file gref.c
 *
 * @brief reference sequence indexer and searcher.
 *
 * @author Hajime Suzuki
 * @date 2016/3/25
 * @license MIT
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "hmap/hmap.h"
#include "psort/psort.h"
#include "zf/zf.h"
#include "gref.h"
#include "kvec.h"
#include "sassert.h"
#include "log.h"

#define UNITTEST_UNIQUE_ID			50

#ifdef TEST
/* use auto-generated main function to run tests */
#define UNITTEST 					1
#define UNITTEST_ALIAS_MAIN			1
#endif

#include "unittest.h"

/* inline directive */
#define _force_inline				inline

/* max, min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MAX3(x, y, z)				( MAX2(MAX2(x, y), z) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )
#define MIN3(x, y, z)				( MIN2(MIN2(x, y), z) )

/* encode and decode id */
#define _rev(_d)					( 0x01 ^ (_d) )
#define _encode_id(_x, _d)			( ((_x)<<1) | (0x01 & (_d)) )
#define _decode_id(_x)				( (_x)>>1 )
#define _decode_dir(_x)				( (_x) & 0x01 )


/**
 * structs and typedefs
 */
typedef kvec_t(struct gref_kmer_tuple_s) kvec_kmer_tuple_t;

/**
 * @struct gref_gid_pair_s
 */
struct gref_gid_pair_s {
	int32_t from;
	int32_t to;
};

/**
 * @struct gref_seq_interval_s
 */
struct gref_seq_interval_s {
	uint64_t base;
	uint64_t tail;
};

/**
 * @struct gref_section_intl_s
 * @brief sizeof(gref_section_intl_s) == 48
 */
struct gref_section_intl_s {
	hmap_header_t header;

	/* forward link index */
	uint32_t fw_link_idx_base;

	/* splitted section link (used to get original name) */
	uint32_t base_id;

	/* gref_section_s compatible */
	struct gref_section_s sec;

	/* reverse link index */
	uint32_t rv_link_idx_base;
	uint32_t reserved2;
};
_static_assert(sizeof(hmap_header_t) == sizeof(struct gref_section_s));
_static_assert(sizeof(struct gref_section_intl_s) == 48);

/**
 * @struct gref_section_half_s
 * @brief former half of the section_intl_s
 */
struct gref_section_half_s {
	/* gref_section_s compatible */
	uint32_t reserved1[2];
	uint64_t reserved2;

	/* section table */
	uint32_t link_idx_base;
	uint32_t reserved3;
};
_static_assert(sizeof(struct gref_section_half_s) == 24);

/**
 * @enum gref_type
 * @breif gref->type
 */
enum gref_type {
	GREF_POOL 	= 1,
	GREF_ACV 	= 2,
	GREF_IDX 	= 3
};

/**
 * @struct gref_s
 * @brief aliased to gref_pool_t, gref_acv_t, and gref_idx_t
 */
struct gref_s {
	/* name -> section mapping */
	hmap_t *hmap;					/* name -> section_info hashmap */
	uint32_t tail_id;

	/* status */
	int8_t type;
	int8_t kmer_available;
	uint8_t reserved2[2];

	/* internal params */
	int64_t iter_init_stack_size;

	/* params */
	struct gref_params_s params;

	/* sequence container */
	kvec_t(uint8_t) seq;

	/* link info container */
	kvec_t(struct gref_gid_pair_s) link;

	/* kv_ptr(link) and link_table shares pointer */
	uint64_t mask;
	int64_t link_table_size;
	uint32_t *link_table;

	/* kmer index container */
	int64_t *kmer_idx_table;
	int64_t kmer_table_size;
	struct gref_gid_pos_s *kmer_table;

	/* sequence encoder */
	struct gref_seq_interval_s (*append_seq)(
		struct gref_s *gref,
		uint8_t const *seq,
		int64_t len);
};
_static_assert(sizeof(struct gref_params_s) == 16);

/**
 * @fn gref_encode_2bit
 * @brief mapping IUPAC amb. to 2bit encoding
 */
static _force_inline
uint8_t gref_encode_2bit(
	int c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x00, C = 0x01, G = 0x02, T = 0x03
	};
	uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('N')] = A,		/* treat 'N' as 'A' */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn gref_encode_4bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t gref_encode_4bit(
	uint8_t c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x01, C = 0x02, G = 0x04, T = 0x08
	};
	uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('R')] = A | G,
		[_b('Y')] = C | T,
		[_b('S')] = G | C,
		[_b('W')] = A | T,
		[_b('K')] = G | T,
		[_b('M')] = A | C,
		[_b('B')] = C | G | T,
		[_b('D')] = A | G | T,
		[_b('H')] = A | C | T,
		[_b('V')] = A | C | G,
		[_b('N')] = 0,		/* treat 'N' as a gap */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b(c)]);

	#undef _b
}

/**
 * @fn gref_copy_seq_ascii
 */
static
struct gref_seq_interval_s gref_copy_seq_ascii(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	uint64_t base = kv_size(gref->seq);
	kv_reserve(gref->seq, base + len);

	/* append */
	for(int64_t i = 0; i < len; i++) {
		kv_at(gref->seq, base + i) = gref_encode_4bit(seq[i]);
	}

	/* resize array */
	kv_size(gref->seq) += len;
	return((struct gref_seq_interval_s){
		.base = base,
		.tail = base + len
	});	
}

/**
 * @fn gref_copy_seq_4bit
 */
static
struct gref_seq_interval_s gref_copy_seq_4bit(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	uint64_t base = kv_size(gref->seq);
	kv_pushm(gref->seq, seq, len);
	return((struct gref_seq_interval_s){
		.base = base,
		.tail = base + len
	});
}

/**
 * @fn gref_nocopy_seq_4bit
 */
static
struct gref_seq_interval_s gref_nocopy_seq_4bit(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	return((struct gref_seq_interval_s){
		.base = (uint64_t)seq,
		.tail = (uint64_t)seq + len
	});
}

/* init / destroy pool */
/**
 * @fn gref_init_pool
 */
gref_pool_t *gref_init_pool(
	gref_params_t const *params)
{
	/* check sanity of params */
	if(params == NULL) {
		return(NULL);
	}
	struct gref_params_s p = *params;

	/* restore defaults */
	#define restore(param, def)		{ (param) = ((uint64_t)(param) == 0) ? (def) : (param); }
	
	restore(p.k, 14);
	restore(p.hash_size, 1024);
	restore(p.seq_format, GREF_ASCII);
	restore(p.copy_mode, GREF_COPY);
	restore(p.index_mode, GREF_INDEX_HASH);
	restore(p.num_threads, 0);

	#undef restore

	/* check sanity */
	if((uint32_t)p.k > 32) { return(NULL); }
	if((uint8_t)p.seq_format > GREF_4BIT) { return(NULL); }
	if((uint8_t)p.copy_mode > GREF_NOCOPY) { return(NULL); }
	if((uint8_t)p.index_mode > GREF_INDEX_ITER) { return(NULL); }

	/* malloc mem */
	struct gref_s *pool = (struct gref_s *)malloc(sizeof(struct gref_s));
	if(pool == NULL) {
		return(NULL);
	}
	memset(pool, 0, sizeof(struct gref_s));

	/* init hmap */
	pool->hmap = hmap_init(p.hash_size, sizeof(struct gref_section_intl_s));
	if(pool->hmap == NULL) {
		goto _gref_init_pool_error_handler;
	}
	pool->tail_id = 0;
	pool->type = GREF_POOL;

	/* calc iterator buffer size */
	pool->iter_init_stack_size = MAX2(1024, (int64_t)pow(3.0, p.k * 0.5));

	/* init seq vector */
	if(p.copy_mode != GREF_NOCOPY) {
		kv_init(pool->seq);
	} else {
		pool->seq.a = NULL;
	}

	/* init link vector */
	kv_init(pool->link);

	/* init seq encoder */
	struct gref_seq_interval_s (*table[][3])(
		struct gref_s *gref,
		uint8_t const *seq,
		int64_t len) = {
		[GREF_ASCII] = {
			[GREF_COPY] = gref_copy_seq_ascii,
			[GREF_4BIT] = NULL
		},
		[GREF_4BIT] = {
			[GREF_COPY] = gref_copy_seq_4bit,
			[GREF_NOCOPY] = gref_nocopy_seq_4bit
		}
	};
	if((pool->append_seq = table[p.seq_format][p.copy_mode]) == NULL) {
		goto _gref_init_pool_error_handler;
	}

	/* copy params */
	pool->params = p;
	return((gref_pool_t *)pool);

_gref_init_pool_error_handler:;
	if(pool != NULL) {
		hmap_clean(pool->hmap);
		if(kv_ptr(pool->seq) != NULL) { kv_destroy(pool->seq); }
		if(kv_ptr(pool->link) != NULL) { kv_destroy(pool->link); }
		free(pool);
	}
	return(NULL);
}

/**
 * @fn gref_clean
 */
void gref_clean(
	gref_t *_gref)
{
	struct gref_s *gref = (struct gref_s *)_gref;

	if(gref != NULL) {
		/* cleanup, cleanup... */
		hmap_clean(gref->hmap);
		if(kv_ptr(gref->seq) != NULL) { kv_destroy(gref->seq); }
		if(gref->link_table != NULL) { free(gref->link_table); }
		if(gref->kmer_idx_table != NULL) { free(gref->kmer_idx_table); }
		if(gref->kmer_table != NULL) { free(gref->kmer_table); }
	}
	return;
}


/* pool modify operation */
/**
 * @fn gref_append_segment
 */
int gref_append_segment(
	gref_pool_t *_pool,
	char const *name,
	int32_t name_len,
	uint8_t const *seq,
	int64_t seq_len)
{
	struct gref_s *pool = (struct gref_s *)_pool;
	debug("append segment");

	/* add sequence at the tail of the seq buffer */
	struct gref_seq_interval_s iv = pool->append_seq(pool, seq, seq_len);

	/* append the first section */
	uint64_t const max_sec_len = 0x80000000;
	uint64_t len = MIN2(iv.tail - iv.base, max_sec_len);

	uint32_t id = hmap_get_id(pool->hmap, name, name_len);
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, id);

	/* update tail_id */
	pool->tail_id = MAX2(pool->tail_id, id + 1);

	/* store section info */
	sec->base_id = id;
	sec->fw_link_idx_base = 0;
	sec->rv_link_idx_base = 0;
	sec->sec = (struct gref_section_s){
		.id = id,
		.len = len,
		.base = iv.base
	};
	return(0);
}

/**
 * @fn gref_append_link
 */
int gref_append_link(
	gref_pool_t *_pool,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori)
{
	struct gref_s *pool = (struct gref_s *)_pool;
	debug("append link");

	/* get ids */
	uint32_t src_id = hmap_get_id(pool->hmap, src, src_len);
	uint32_t dst_id = hmap_get_id(pool->hmap, dst, dst_len);

	/* add forward link */
	kv_push(pool->link, ((struct gref_gid_pair_s){
		.from = _encode_id(src_id, src_ori),
		.to = _encode_id(dst_id, dst_ori)
	}));

	/* add reverse link */
	kv_push(pool->link, ((struct gref_gid_pair_s){
		.from = _encode_id(dst_id, _rev(dst_ori)),
		.to = _encode_id(src_id, _rev(src_ori))
	}));

	/* update tail_id */
	pool->tail_id = MAX3(pool->tail_id, src_id, dst_id);
	return(0);
}

/**
 * @fn gref_split_section
 * @brief split base section and give new name (splitted) to the latter section.
 */
int gref_split_section(
	gref_pool_t *_pool,
	char const *base,
	int32_t base_len,
	int64_t pos,
	char const *splitted,
	int32_t splitted_len)
{
	return(0);
}


/* build link table (pool -> acv conversion) */

/**
 * @fn gref_add_tail_section
 */
static _force_inline
void gref_add_tail_section(
	struct gref_s *pool)
{
	uint32_t tail_id = pool->tail_id + 1;
	if(hmap_get_count(pool->hmap) >= tail_id) {
		/* sentinel already exists */
		return;
	}

	/* push sentinel with unique name */
	char const *template = "tail_sentinel_";
	int64_t len = strlen(template);

	char buf[256];
	strcpy(buf, template);
	uint32_t id = (uint32_t)-1;
	do {
		buf[len++] = '0';
		buf[len] = '\0';

		/* push sentinel to section array */
		id = hmap_get_id(pool->hmap, buf, len);
	} while(id != tail_id && len < 256);

	/* set info */
	struct gref_section_intl_s *tail_sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, tail_id);
	tail_sec->base_id = tail_id;
	tail_sec->sec = (struct gref_section_s){
		.id = tail_id,
		.len = 0,
		.base = 0
	};
	return;
}

/**
 * @fn gref_build_link_idx_table
 * @brief build link_idx table. gref_add_tail_section must be called beforehand.
 */
static _force_inline
int gref_build_link_idx_table(
	struct gref_s *pool)
{
	int64_t link_idx_table_size = 2 * (pool->tail_id + 1);
	int64_t link_table_size = kv_size(pool->link);

	/* store link table size */
	pool->link_table_size = link_table_size;

	/* sort by src, build src->dst mapping */
	debug("sort src->dst mapping, size(%llu)", kv_size(pool->link));
	if(psort_half(kv_ptr(pool->link), kv_size(pool->link),
		sizeof(struct gref_gid_pair_s), pool->params.num_threads) != 0) {

		/* sort failed */
		return(-1);
	}

	/* init index */
	uint32_t prev_gid = 0;
	struct gref_section_half_s *sec_half =
		(struct gref_section_half_s *)hmap_get_object(pool->hmap, 0);

	/* store link index */
	sec_half[prev_gid].link_idx_base = 0;		/* head is always zero */
	for(int64_t i = 0; i < link_table_size; i++) {
		uint32_t gid = kv_at(pool->link, i).from;
		debug("i(%lld), gid(%u), prev_gid(%u)", i, gid, prev_gid);

		if(prev_gid == gid) { continue; }

		debug("index table for gid(%u) ends at %lld, next_gid(%u)", prev_gid, i, gid);

		/* sequence update detected */
		for(int64_t j = prev_gid + 1; j < gid + 1; j++) {
			debug("fill gaps j(%lld)", j);
			sec_half[j].link_idx_base = i;
		}
		prev_gid = gid;
	}

	/* store tail */
	for(int64_t j = prev_gid + 1; j < link_idx_table_size + 1; j++) {
		debug("fill gaps j(%lld)", j);
		sec_half[j].link_idx_base = link_table_size;
	}
	return(0);
}

/**
 * @fn gref_shrink_link_table
 * @brief shrink link table. gref_build_link_idx_table must be called beforehand.
 */
static _force_inline
int gref_shrink_link_table(
	struct gref_s *pool)
{
	int64_t link_table_size = pool->link_table_size;
	uint32_t *link_table = (uint32_t *)kv_ptr(pool->link);

	/* pack */
	for(int64_t i = 0; i < link_table_size; i++) {
		link_table[i] = kv_at(pool->link, i).to;
	}
	kv_resize(pool->link, link_table_size / 2);
	
	/* store info */
	pool->link_table = (uint32_t *)kv_ptr(pool->link);
	return((pool->link_table == NULL) ? -1 : 0);
}

/**
 * @fn gref_expand_link_table
 */
static _force_inline
int gref_expand_link_table(
	struct gref_s *acv)
{
	/* resize mem */
	int64_t link_table_size = acv->link_table_size;
	kv_resize(acv->link, link_table_size);
	uint32_t *link = (uint32_t *)kv_ptr(acv->link);

	if(link == NULL) { return(-1); }

	/* load section ptr */
	int64_t sec_cnt = acv->tail_id + 1;
	struct gref_section_half_s *sec_half =
		(struct gref_section_half_s *)hmap_get_object(acv->hmap, 0);

	/* expand, iterate from tail to head */
	for(int64_t i = 2 * sec_cnt - 1; i >= 0; i--) {
		for(int64_t j = sec_half[i + 1].link_idx_base - 1;
			j >= sec_half[i].link_idx_base;
			j--) {

			kv_at(acv->link, j) = (struct gref_gid_pair_s){
				.from = i,
				.to = link[j]
			};
		}
	}
	return(0);
}

/**
 * @fn gref_freeze_pool
 */
gref_acv_t *gref_freeze_pool(
	gref_pool_t *pool)
{
	struct gref_s *gref = (struct gref_s *)pool;

	if(gref == NULL || gref->type != GREF_POOL) {
		goto _gref_freeze_pool_error_handler;
	}

	/* push tail sentinel */
	gref_add_tail_section(gref);

	/* build link array */
	if(gref_build_link_idx_table(gref) != 0) {
		/* sort failed */
		goto _gref_freeze_pool_error_handler;
	}

	/* shrink table */
	if(gref_shrink_link_table(gref) != 0) {
		goto _gref_freeze_pool_error_handler;
	}

	/* change type */
	gref->type = GREF_ACV;
	return((gref_acv_t *)gref);

_gref_freeze_pool_error_handler:;
	gref_clean((void *)gref);
	return(NULL);
}

/**
 * @fn gref_melt_archive
 */
gref_pool_t *gref_melt_archive(
	gref_acv_t *acv)
{
	struct gref_s *gref = (struct gref_s *)acv;

	if(gref == NULL || gref->type != GREF_ACV) {
		goto _gref_melt_archive_error_handler;
	}

	/* expand table */
	if(gref_expand_link_table(gref) != 0) {
		goto _gref_melt_archive_error_handler;
	}

	/* change type */
	gref->type = GREF_POOL;
	return((gref_pool_t *)gref);

_gref_melt_archive_error_handler:;
	gref_clean((void *)gref);
	return(NULL);
}


/* kmer enumeration */
/**
 * @struct gref_iter_stack_s
 * @brief aliased to gref_iter_t
 */
struct gref_iter_stack_s {
	/* previous stack */
	struct gref_iter_stack_s *prev_stack;

	/* global params */
	uint32_t global_rem_len;
	uint32_t shift_len;

	/* current section info */
	uint32_t sec_gid;
	uint32_t link_idx;

	/* sequence info */
	uint8_t (*fetch)(struct gref_iter_stack_s *stack);
	uint8_t const *seq_base;
	uint32_t rem_len;
	uint32_t pos;

	/* kmer table */
	uint64_t cnt_arr;
	uint32_t kmer_idx;
	uint32_t kmer_occ;
	uint64_t kmer[];
};

/**
 * @struct gref_iter_s
 */
#define GREF_ITER_INTL_MEM_ARR_LEN			( 5 )
struct gref_iter_s {
	/* global info */
	uint32_t base_gid;
	uint32_t tail_gid;
	uint32_t reserved1;
	uint8_t seed_len;
	uint8_t shift_len;
	uint16_t reserved3;

	uint8_t const *seq;
	uint32_t const *link_table;
	struct gref_section_intl_s const *sec;

	/* stack mem array */
	struct gref_iter_stack_s *stack, *root_stack;
	void *mem_arr[GREF_ITER_INTL_MEM_ARR_LEN];
};
_static_assert(sizeof(struct gref_iter_s) == 96);

/**
 * @fn gref_iter_add_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_add_stack(
	struct gref_iter_s *iter,
	struct gref_iter_stack_s *stack)
{
	/* fixme: check remaining memory size and malloc if there is no room */
	struct gref_iter_stack_s *new_stack = (struct gref_iter_stack_s *)(
		(uint64_t *)(stack + 1) + stack->kmer_occ);

	debug("stack(%p), kmer_occ(%u), new_stack(%p)",
		stack, stack->kmer_occ, new_stack);

	/* link to previous stack */
	new_stack->prev_stack = stack;

	/* copy from previous stack */
	new_stack->global_rem_len = stack->global_rem_len;

	/* constant */
	new_stack->shift_len = stack->shift_len;

	/* copy buffer */
	new_stack->pos = stack->pos;
	new_stack->cnt_arr = stack->cnt_arr;
	new_stack->kmer_idx = 0;
	new_stack->kmer_occ = stack->kmer_occ;
	memcpy(&new_stack->kmer, &stack->kmer, stack->kmer_occ * sizeof(uint64_t));
	return(new_stack);
}

/**
 * @fn gref_iter_remove_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_remove_stack(
	struct gref_iter_stack_s *stack)
{
	return(stack->prev_stack);
}

/**
 * @fn gref_iter_append_base
 */
static _force_inline
int gref_iter_append_base(
	struct gref_iter_stack_s *stack,
	uint8_t c)
{
	/* conversion tables */
	uint8_t const popcnt_table[] = {
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 0	/* ignore 0x0f */
	};
	enum bases {
		A = 0x00, C = 0x01, G = 0x02, T = 0x03
	};
	uint8_t const encode_2bit[][3] = {
		{},
		{ A },
		{ C },
		{ A, C },
		{ G },
		{ A, G },
		{ C, G },
		{ A, C, G },
		{ T },
		{ A, T },
		{ C, T },
		{ A, C, T },
		{ G, T },
		{ A, G, T },
		{ C, G, T },
		{},
	};

	uint64_t pcnt = popcnt_table[c];
	stack->cnt_arr = (stack->cnt_arr<<2) | pcnt;
	uint64_t occ = stack->kmer_occ;

	/* branch */
	switch(3 - pcnt) {
		case 0: memcpy(&stack->kmer[2 * occ], stack->kmer, sizeof(uint64_t) * occ);
		case 1: memcpy(&stack->kmer[occ], stack->kmer, sizeof(uint64_t) * occ);
		/* fall through */
		default: break;		/* return(-1); */
	}

	/* append to vector */
	for(int64_t j = 0; j < pcnt; j++) {
		for(int64_t k = 0; k < occ; k++) {
			stack->kmer[j * occ + k] =
				(stack->kmer[j * occ + k]>>2) | (encode_2bit[c][j]<<stack->shift_len);
			debug("%lld, %lld, %lld, %x, %llx",
				j, k, j * occ + k, encode_2bit[c][j], stack->kmer[j * occ + k]);
		}
	}

	/* update occ */
	occ *= pcnt;

	/* merge (shrink buffer) */
	uint64_t shrink_skip = 0x03 & (stack->cnt_arr >> (stack->shift_len + 2));
	debug("cnt_arr(%llx), occ(%llu), shrink_skip(%llu)",
		stack->cnt_arr, occ, shrink_skip);
	if(shrink_skip > 1) {
		occ /= shrink_skip;
		for(int64_t j = 0; j < occ; j++) {
			stack->kmer[j] = stack->kmer[j * shrink_skip];
		}
	}

	/* write back occupancy */
	stack->kmer_idx = 0;
	stack->kmer_occ = occ;
	return(0);
}

/**
 * @fn gref_iter_calc_seq_base
 * @brief len: length to fetch from base
 */
static _force_inline
uint8_t const *gref_iter_calc_seq_base(
	struct gref_iter_s *p,
	uint32_t gid,
	uint32_t len)
{
	debug("gid(%u), id(%u)", gid, _decode_id(gid));
	int64_t pos = p->sec[_decode_id(gid)].sec.base;

	for(int64_t i = 0; i < p->sec[_decode_id(gid)].sec.len; i++) {
		debug("c(%u)", p->seq[pos + i]);
	}

	pos += (_decode_dir(gid) == 0)
		? (len - 1)									/* forward */
		: (p->sec[_decode_id(gid)].sec.len - len);	/* reverse */
	debug("base(%lld), dir(%x), len_ofs(%u), pos(%lld)",
		p->sec[_decode_id(gid)].sec.base,
		_decode_dir(gid),
		(_decode_dir(gid) == 0) ? (len - 1) : (p->sec[_decode_id(gid)].sec.len - len),
		pos);
	return(p->seq + pos);
}

/**
 * @fn gref_iter_foward_fetch, gref_iter_reverse_fetch
 */
static
uint8_t gref_iter_foward_fetch(
	struct gref_iter_stack_s *stack)
{
	debug("rem_len(%u), pos(%u), c(%x)", stack->rem_len - 1, stack->pos + 1, stack->seq_base[-((int64_t)(stack->rem_len - 1))]);
	stack->pos++;
	return(stack->seq_base[-((int64_t)(--stack->rem_len))]);
}
static
uint8_t gref_iter_reverse_fetch(
	struct gref_iter_stack_s *stack)
{
	debug("rem_len(%u), pos(%u), c(%x)", stack->rem_len - 1, stack->pos + 1, stack->seq_base[-((int64_t)(stack->rem_len - 1))]);
	stack->pos++;
	return(stack->seq_base[(int64_t)(--stack->rem_len)]);
}

/**
 * @fn gref_iter_fetch
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_fetch(
	struct gref_iter_s *iter,
	struct gref_iter_stack_s *stack)
{
	if(stack->rem_len > 0) {
		/* fetch seq */
		gref_iter_append_base(stack, stack->fetch(stack));
		return(stack);
	} else if(stack->rem_len == 0) {
		/* return if no more seq remains */
		if(stack->global_rem_len == 0) {
			debug("reached tail");
			stack = stack->prev_stack;
		}

		/* remove stack if no more link remains */
		debug("stack(%p), gid(%u), link_idx(%u), rv_link_idx_base(%u)",
			stack, stack->sec_gid, stack->link_idx, iter->sec[_decode_id(stack->sec_gid)].rv_link_idx_base);
		while(stack->link_idx == iter->sec[_decode_id(stack->sec_gid)].rv_link_idx_base) {
			debug("stack(%p), gid(%u), link_idx(%u), rv_link_idx_base(%u)",
				stack, stack->sec_gid, stack->link_idx, iter->sec[_decode_id(stack->sec_gid)].rv_link_idx_base);
			if((stack = gref_iter_remove_stack(stack)) == NULL) {
				/* root (reached the end of enumeration) */
				debug("reached NULL");
				return(NULL);
			}
		}

		/* load next gid */
		uint32_t gid = iter->link_table[stack->link_idx++];

		/* add stack for the next section */
		stack = gref_iter_add_stack(iter, stack);
		debug("stack added");

		/* init link info */
		stack->sec_gid = gid;
		stack->link_idx = iter->sec[_decode_id(gid)].fw_link_idx_base;
		// debug("link_table(%p), link_idx(%u)", iter->link_table, stack->link_idx);

		/* init seq info */
		stack->fetch = (_decode_dir(gid) == 0)
			? gref_iter_foward_fetch
			: gref_iter_reverse_fetch;
		stack->rem_len = MIN2(stack->global_rem_len, iter->sec[_decode_id(gid)].sec.len);
		stack->seq_base = gref_iter_calc_seq_base(iter, gid, stack->rem_len);

		/* adjust global_rem_len */
		stack->global_rem_len -= stack->rem_len;

		debug("gid(%u), pos(%u), rem_len(%u), global_rem_len(%u)",
			gid, stack->pos, stack->rem_len, stack->global_rem_len);

		/* fetch seq */
		gref_iter_append_base(stack, stack->fetch(stack));
		return(stack);
	}

	/* never reaches here */
	return(NULL);
}

/**
 * @fn gref_iter_init_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_init_stack(
	struct gref_iter_s *iter,
	struct gref_iter_stack_s *stack)
{
	uint32_t gid = iter->base_gid;
	uint32_t len = iter->sec[_decode_id(gid)].sec.len;

	/* init stack */
	stack->prev_stack = NULL;
	// stack->global_rem_len = len;
	stack->global_rem_len = iter->seed_len - 1;
	stack->shift_len = iter->shift_len;

	/* current section info */
	stack->sec_gid = gid;
	stack->link_idx = iter->sec[_decode_id(gid)].fw_link_idx_base;

	/* seq info */
	stack->fetch = gref_iter_foward_fetch;
	stack->seq_base = gref_iter_calc_seq_base(iter, gid, len);
	stack->rem_len = len;
	stack->pos = 0;

	/* init buffer and counter */
	stack->cnt_arr = 0;
	stack->kmer_idx = 0;
	stack->kmer_occ = 1;
	stack->kmer[0] = 0;

	for(int64_t i = 0; i < iter->seed_len; i++) {
		stack = gref_iter_fetch(iter, stack);
	}
	return(stack);
}


/**
 * @fn gref_iter_init
 */
gref_iter_t *gref_iter_init(
	gref_acv_t *acv)
{
	struct gref_s *gref = (struct gref_s *)acv;

	debug("init_stack_size(%lld)", gref->iter_init_stack_size);

	/* malloc mem */
	struct gref_iter_s *iter = (struct gref_iter_s *)malloc(
		sizeof(struct gref_iter_s) + sizeof(uint64_t) * gref->iter_init_stack_size);
	if(iter == NULL) {
		return(NULL);
	}

	/* init param container */
	memset(iter->mem_arr, 0, sizeof(void *) * GREF_ITER_INTL_MEM_ARR_LEN);

	/* iterate from section 0 */
	iter->base_gid = 0;
	iter->tail_gid = _encode_id(acv->tail_id, 0);
	
	/* set params */
	iter->seed_len = acv->params.k;
	iter->shift_len = 2 * (acv->params.k - 1);
	iter->seq = kv_ptr(acv->seq);
	iter->link_table = acv->link_table;
	iter->sec = (struct gref_section_intl_s const *)hmap_get_object(acv->hmap, 0);

	/* init stack */
	iter->stack = iter->root_stack = gref_iter_init_stack(iter,
		(struct gref_iter_stack_s *)(iter + 1));

	debug("stack(%p), iter(%p)", iter->stack, iter);
	return((gref_iter_t *)iter);
}

/**
 * @fn gref_iter_next
 */
struct gref_kmer_tuple_s gref_iter_next(
	gref_iter_t *_iter)
{
	struct gref_iter_s *iter = (struct gref_iter_s *)_iter;
	struct gref_iter_stack_s *stack = iter->stack;

	#define return_kmer(_iter, _stack) { \
		debug("return kmer(%llx), gid(%u), pos(%u)", \
			(_stack)->kmer[(_stack)->kmer_idx], \
			(_stack)->sec_gid, \
			(_stack)->pos); \
		return((struct gref_kmer_tuple_s){ \
			.kmer = (_stack)->kmer[(_stack)->kmer_idx++], \
			.pos = (struct gref_gid_pos_s){ \
				.gid = (_iter)->base_gid, \
				.pos = (_stack)->pos - (_iter)->seed_len \
			} \
		}); \
	}

	debug("stack(%p), kmer_idx(%u), kmer_occ(%u)", stack, stack->kmer_idx, stack->kmer_occ);
	if(stack->kmer_idx < stack->kmer_occ) {
		return_kmer(iter, stack);
	} else if((iter->stack = stack = gref_iter_fetch(iter, stack)) != NULL) {
		return_kmer(iter, stack);
	} else if((iter->base_gid += 2) < iter->tail_gid) {
		debug("base_gid(%u), tail_gid(%u)", iter->base_gid, iter->tail_gid);
		iter->stack = stack = gref_iter_init_stack(iter, iter->root_stack);
		return_kmer(iter, stack);
	}

	#undef return_kmer

	debug("terminal");
	/* reached the end */
	return((struct gref_kmer_tuple_s){
		.kmer = GREF_ITER_KMER_TERM,
		.pos = (struct gref_gid_pos_s){
			.gid = (uint32_t)-1,
			.pos = 0
		}
	});
}

/**
 * @fn gref_iter_clean
 */
void gref_iter_clean(
	gref_iter_t *_iter)
{
	struct gref_iter_s *iter = (struct gref_iter_s *)_iter;
	struct gref_iter_stack_s *stack = iter->stack;
	debug("stack(%p), iter(%p)", stack, iter);

	if(stack != NULL) {
		for(int64_t i = 0; i < GREF_ITER_INTL_MEM_ARR_LEN; i++) {
			if(iter->mem_arr[i] != NULL) { free(iter->mem_arr[i]); }
		}
		free((void *)iter);
	}
	return;
}


/* build kmer index (acv -> idx conversion) */
/**
 * @fn gref_build_kmer_idx_table
 */
static _force_inline
int64_t *gref_build_kmer_idx_table(
	struct gref_s *acv,
	struct gref_kmer_tuple_s *arr,
	int64_t size)
{
	kvec_t(int64_t) kmer_idx;
	kv_init(kmer_idx);

	uint64_t kmer_idx_size = 0x01 << (2 * acv->params.k);
	kv_reserve(kmer_idx, kmer_idx_size);

	uint64_t prev_kmer = 0;
	kv_push(kmer_idx, prev_kmer);
	for(int64_t i = 0; i < size; i++) {
		uint64_t kmer = arr[i].kmer;
		debug("i(%lld), kmer(%llx), id(%u), pos(%u), prev_kmer(%llx)",
			i, kmer, arr[i].pos.gid, arr[i].pos.pos, prev_kmer);

		if(prev_kmer == kmer) { continue; }

		debug("index table for kmer(%llx) ends at %lld, next_kmer(%llx)", prev_kmer, i, kmer);

		/* sequence update detected */
		for(uint64_t j = prev_kmer + 1; j < kmer + 1; j++) {
			debug("fill gaps j(%llx)", j);
			kv_push(kmer_idx, i);
		}
		prev_kmer = kmer;
	}
	for(uint64_t j = prev_kmer; j < kmer_idx_size; j++) {
		debug("fill tail gaps j(%llx)", j);
		kv_push(kmer_idx, size);
	}
	return(kv_ptr(kmer_idx));
}

/**
 * @fn gref_shrink_kmer_table
 */
static _force_inline
struct gref_gid_pos_s *gref_shrink_kmer_table(
	struct gref_s *acv,
	struct gref_kmer_tuple_s *kmer_table,
	int64_t kmer_table_size)
{
	struct gref_gid_pos_s *packed_pos = (struct gref_gid_pos_s *)kmer_table;
	
	for(int64_t i = 0; i < kmer_table_size; i++) {
		packed_pos[i] = kmer_table[i].pos;
	}

	return(realloc(kmer_table, sizeof(struct gref_gid_pos_s) * kmer_table_size));
}

/**
 * @fn gref_build_index
 */
gref_idx_t *gref_build_index(
	gref_acv_t *acv)
{
	struct gref_s *gref = (struct gref_s *)acv;

	if(gref == NULL || gref->type != GREF_ACV) {
		goto _gref_build_index_error_handler;
	}

	/* enumerate kmers and pack into vector */
	kvec_t(struct gref_kmer_tuple_s) v;
	kv_init(v);

	struct gref_iter_s *iter = gref_iter_init(acv);
	struct gref_kmer_tuple_s t;
	while((t = gref_iter_next(iter)).kmer != GREF_ITER_KMER_TERM) {
		kv_push(v, t);
	}
	gref_iter_clean(iter);

	/* sort kmers */
	if(psort_half(kv_ptr(v), kv_size(v),
		sizeof(struct gref_kmer_tuple_s), acv->params.num_threads) != 0) {
		if(kv_ptr(v) != NULL) { kv_destroy(v); }
		goto _gref_build_index_error_handler;
	}

	/* build index of kmer table */
	gref->kmer_idx_table = gref_build_kmer_idx_table(acv, kv_ptr(v), kv_size(v));
	if(gref->kmer_idx_table == NULL) {
		goto _gref_build_index_error_handler;
	}

	/* shrink table */
	gref->kmer_table_size = kv_size(v);
	gref->kmer_table = gref_shrink_kmer_table(acv, kv_ptr(v), kv_size(v));
	if(gref->kmer_table == NULL) {
		goto _gref_build_index_error_handler;
	}

	/* store misc constants for kmer matching */
	gref->mask = (0x01<<(2 * gref->params.k)) - 1;

	/* change state */
	gref->type = GREF_IDX;
	return((gref_idx_t *)gref);

_gref_build_index_error_handler:;
	gref_clean((gref_t *)gref);
	return(NULL);
}

/**
 * @fn gref_disable_index
 */
gref_acv_t *gref_disable_index(
	gref_idx_t *idx)
{
	struct gref_s *gref = (struct gref_s *)idx;

	if(gref == NULL || gref->type != GREF_IDX) {
		gref_clean((gref_t *)gref);
		return(NULL);
	}

	/* cleanup kmer table */
	if(idx->kmer_idx_table != NULL) {
		free(idx->kmer_idx_table);
		idx->kmer_idx_table = NULL;
	}
	if(idx->kmer_table != NULL) {
		idx->kmer_table = NULL;
	}
	idx->kmer_table_size = 0;

	/* change state */
	gref->type = GREF_ACV;
	return((gref_acv_t *)gref);
}

/**
 * @fn gref_match_2bitpacked
 */
struct gref_match_res_s gref_match_2bitpacked(
	gref_idx_t const *_gref,
	uint64_t seq)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	seq &= gref->mask;
	int64_t base = gref->kmer_idx_table[seq];
	int64_t tail = gref->kmer_idx_table[seq + 1];

	debug("seq(%llx), mask(%llx), base(%lld), tail(%lld)",
		seq, gref->mask, base, tail);
	return((struct gref_match_res_s){
		.ptr = &gref->kmer_table[base],
		.len = tail - base
	});
}

/**
 * @fn gref_match
 * @brief seq length must be equal to k.
 */
struct gref_match_res_s gref_match(
	gref_idx_t const *_gref,
	uint8_t const *seq)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	int64_t const seed_len = gref->params.k;
	int64_t const shift_len = 2 * (seed_len - 1);

	uint64_t packed_seq = 0;
	for(int64_t i = 0; i < seed_len; i++) {
		packed_seq = (packed_seq>>2) | (gref_encode_2bit(seq[i])<<shift_len);
	}
	return(gref_match_2bitpacked((gref_t const *)gref, packed_seq));
}

/* misc */
#if 0
/**
 * @fn gref_dump_index
 */
int gref_dump_index(
	gref_t const *gref,
	zf_t *outfp)
{
g	return(0);
}

/**
 * @fn gref_load_index
 */
gref_t *gref_load_index(
	zf_t *infp)
{
	return(NULL);
}
#endif

/**
 * @fn gref_get_section_cnt
 */
int64_t gref_get_section_cnt(
	gref_t const *_gref)
{
	struct gref_s *gref = (struct gref_s *)_gref;
	return((int64_t)gref->tail_id);
}

/**
 * @fn gref_get_section
 */
struct gref_section_s const *gref_get_section(
	gref_t const *_gref,
	uint32_t id)
{
	struct gref_s *gref = (struct gref_s *)_gref;

	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(gref->hmap, id);
	return((struct gref_section_s const *)&sec->sec);
}

/**
 * @fn gref_get_name
 */
struct gref_idx_str_s gref_get_name(
	gref_t const *_gref,
	uint32_t id)
{
	struct gref_s *gref = (struct gref_s *)_gref;
	struct hmap_key_s key = hmap_get_key(gref->hmap, id);
	return((struct gref_idx_str_s){
		.str = key.str,
		.len = key.len
	});
}

/**
 * @fn gref_get_ptr
 */
uint8_t const *gref_get_ptr(
	gref_t const *_gref)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	return((uint8_t const *)kv_ptr(gref->seq));
}

/**
 * @fn gref_get_total_len
 */
int64_t gref_get_total_len(
	gref_t const *_gref)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	return(kv_size(gref->seq));
}


/**
 * unittests
 */
#ifdef TEST

unittest_config(
	.name = "ref",
	.depends_on = { "psort", "hmap" }
);

#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

#define _pack(x) ({ \
	uint64_t _packed_seq = 0; \
	int64_t _len = strlen(x); \
	uint8_t _shift_len = 2 * (_len - 1); \
	for(int64_t i = 0; i < _len; i++) { \
		_packed_seq = (_packed_seq>>2) | (gref_encode_2bit((x)[i])<<_shift_len); \
	} \
	_packed_seq; \
})

/* make pool context */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 3));

	assert(pool != NULL);

	gref_clean(pool);
}

/* add segment */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));

	int ret = gref_append_segment(pool, _str("sec0"), _seq("AARA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_segment(pool, _str("sec1"), _seq("MAAA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_segment(pool, _str("sec2"), _seq("ACGT"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	/* section count */
	assert(gref_get_section_cnt(pool) == 3, "%lld", gref_get_section_cnt(pool));

	/* total len */
	assert(gref_get_total_len(pool) == 12, "len(%lld)", gref_get_total_len(pool));

	gref_clean(pool);
}

/* archive */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_cnt(acv) == 3, "%lld", gref_get_section_cnt(acv));

	/* total len */
	assert(gref_get_total_len(acv) == 13, "len(%lld)", gref_get_total_len(acv));

	gref_clean(acv);
}

/* seed iteration */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	
	/* build archive */
	gref_acv_t *acv = gref_freeze_pool(pool);

	/* create iterator */
	gref_iter_t *iter = gref_iter_init(acv);
	assert(iter != NULL, "%p", iter);

	/* enumerate */
	#define _f(_i)					gref_iter_next(_i)
	#define _check_kmer(_t, _k, _id, _pos) ( \
		   (_t).kmer == _pack(_k) \
		&& (_t).pos.gid == _encode_id(_id, 0) \
		&& (_t).pos.pos == (_pos) \
	)
	#define _print_kmer(_t) \
		"kmer(%llx), sec(%u), pos(%u)", \
		(_t).kmer, _decode_id((_t).pos.gid), (_t).pos.pos

	struct gref_kmer_tuple_s t;

	/* sec0 */
	t = _f(iter); assert(_check_kmer(t, "GGA", 0, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGG", 0, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAA", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGA", 0, 1), _print_kmer(t));
	
	/* sec0-sec1 */
	t = _f(iter); assert(_check_kmer(t, "AAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AAC", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAC", 0, 2), _print_kmer(t));
	
	/* sec0-sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AAA", 0, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACA", 0, 3), _print_kmer(t));

	/* sec0-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AAC", 0, 3), _print_kmer(t));
	
	/* sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AAC", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CAC", 1, 0), _print_kmer(t));

	/* sec2 */
	t = _f(iter); assert(_check_kmer(t, "ACA", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACC", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACG", 2, 0), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "CAA", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCA", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGA", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CAC", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCC", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGC", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CAG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGG", 2, 1), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CAG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAG", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GCG", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGG", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGG", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGT", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGT", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGT", 2, 3), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "GTG", 2, 4), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "TGT", 2, 5), _print_kmer(t));

	t = _f(iter); assert(t.kmer == GREF_ITER_KMER_TERM, "%llx, %llx", t.kmer, GREF_ITER_KMER_TERM);

	#undef _f
	#undef _check_kmer
	#undef _print_kmer

	gref_iter_clean(iter);
	gref_clean(acv);
}

/* build index */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */	
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	gref_idx_t *idx = gref_build_index(pool);
	assert(idx != NULL, "idx(%p)", idx);

	/* pointer to seq */
	uint8_t const *ptr = gref_get_ptr(idx);
	assert(ptr != NULL, "ptr(%p)", ptr);

	/* section count */
	assert(gref_get_section_cnt(idx) == 3, "%lld", gref_get_section_cnt(idx));

	/* total len */
	assert(gref_get_total_len(idx) == 16, "len(%lld)", gref_get_total_len(idx));

	gref_clean(idx);
}

/* get_section */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_acv_t *acv = gref_freeze_pool(pool);
	gref_idx_t *idx = gref_build_index(pool);

	/* section id is given in ascending order from 0 */
	assert(gref_get_section(idx, 0) != NULL, "%p", gref_get_section(idx, 0));
	assert(gref_get_section(idx, 0)->id == 0, "id(%u)", gref_get_section(idx, 0)->id);
	assert(gref_get_section(idx, 0)->len == 4, "len(%u)", gref_get_section(idx, 0)->len);
	assert(gref_get_section(idx, 0)->base == 0, "base(%llu)", gref_get_section(idx, 0)->base);

	/* section 1 */
	assert(gref_get_section(idx, 1) != NULL, "%p", gref_get_section(idx, 1));
	assert(gref_get_section(idx, 1)->id == 1, "id(%u)", gref_get_section(idx, 1)->id);
	assert(gref_get_section(idx, 1)->len == 4, "len(%u)", gref_get_section(idx, 1)->len);
	assert(gref_get_section(idx, 1)->base == 4, "base(%llu)", gref_get_section(idx, 1)->base);

	/* section 2 */
	assert(gref_get_section(idx, 2) != NULL, "%p", gref_get_section(idx, 2));
	assert(gref_get_section(idx, 2)->id == 2, "id(%u)", gref_get_section(idx, 2)->id);
	assert(gref_get_section(idx, 2)->len == 8, "len(%u)", gref_get_section(idx, 2)->len);
	assert(gref_get_section(idx, 2)->base == 8, "base(%llu)", gref_get_section(idx, 2)->base);

	gref_clean(idx);
}

/* get_name */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_acv_t *acv = gref_freeze_pool(pool);
	gref_idx_t *idx = gref_build_index(pool);

	/* section id is given in ascending order from 0 */
	assert(gref_get_name(idx, 0).len == 4, "%d", gref_get_name(idx, 0).len);
	assert(strcmp(gref_get_name(idx, 0).str, "sec0") == 0, "%s", gref_get_name(idx, 0).str);

	/* section 1 */
	assert(gref_get_name(idx, 1).len == 4, "%d", gref_get_name(idx, 1).len);
	assert(strcmp(gref_get_name(idx, 1).str, "sec1") == 0, "%s", gref_get_name(idx, 1).str);

	/* section 2 */
	assert(gref_get_name(idx, 2).len == 4, "%d", gref_get_name(idx, 2).len);
	assert(strcmp(gref_get_name(idx, 2).str, "sec2") == 0, "%s", gref_get_name(idx, 2).str);

	gref_clean(idx);
}

/* match */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(.k = 3));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_acv_t *acv = gref_freeze_pool(pool);
	gref_idx_t *idx = gref_build_index(pool);

	/* without ambiguous bases */
	struct gref_match_res_s r = gref_match(idx, (uint8_t const *)"GTG");
	assert(r.ptr != NULL, "%p", r.ptr);
	assert(r.len == 1, "%lld", r.len);

	/* check pos */
	assert(r.ptr[0].pos == 4, "%u", r.ptr[0].pos);

	/* check section */
	struct gref_section_s const *sec = gref_get_section(idx, gref_id(r.ptr[0].gid));
	assert(sec->id == 2, "id(%u)", sec->id);
	assert(sec->len == 8, "len(%u)", sec->len);
	assert(sec->base == 8, "base(%llu)", sec->base);


	/* with ambiguous bases */
	r = gref_match(idx, (uint8_t const *)"GGG");
	assert(r.ptr != NULL, "%p", r.ptr);
	assert(r.len == 3, "%lld", r.len);

	/* check pos */
	assert(r.ptr[0].pos == 0, "%u", r.ptr[0].pos);

	/* check section */
	sec = gref_get_section(idx, gref_id(r.ptr[0].gid));
	assert(sec->id == 0, "id(%u)", sec->id);
	assert(sec->len == 4, "len(%u)", sec->len);
	assert(sec->base == 0, "base(%llu)", sec->base);

	/* check pos */
	assert(r.ptr[1].pos == 1, "%u", r.ptr[1].pos);

	/* check section */
	sec = gref_get_section(idx, gref_id(r.ptr[1].gid));
	assert(sec->id == 1, "id(%u)", sec->id);
	assert(sec->len == 4, "len(%u)", sec->len);
	assert(sec->base == 4, "base(%llu)", sec->base);

	/* check pos */
	assert(r.ptr[2].pos == 2, "%u", r.ptr[2].pos);

	/* check section */
	sec = gref_get_section(idx, gref_id(r.ptr[2].gid));
	assert(sec->id == 2, "id(%u)", sec->id);
	assert(sec->len == 8, "len(%u)", sec->len);
	assert(sec->base == 8, "base(%llu)", sec->base);

	gref_clean(idx);
}

#endif

/**
 * end of gref.c
 */
