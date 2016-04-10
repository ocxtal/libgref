
/**
 * @file ref.c
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
#include "ref.h"
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

#define _force_inline				inline

/* max, min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )

/* encode and decode id */
#define _rev(_d)					( 0x01 ^ (_d) )
#define _encode_id(_x, _d)			( ((_x)<<1) | (0x01 & (_d)) )
#define _decode_id(_x)				( (_x)>>1 )
#define _decode_dir(_x)				( (_x) & 0x01 )


/**
 * structs and typedefs
 */

/**
 * @struct ref_gid_pair_s
 */
struct ref_gid_pair_s {
	int32_t from;
	int32_t to;
};

/**
 * @struct ref_seq_interval_s
 */
struct ref_seq_interval_s {
	uint64_t base;
	uint64_t tail;
};

/**
 * @struct ref_section_intl_s
 * @brief sizeof(ref_section_intl_s) == 48
 */
struct ref_section_intl_s {
	hmap_header_t header;

	/* forward link index */
	uint32_t fw_link_idx_base;

	/* splitted section link (used to get original name) */
	uint32_t base_id;

	/* ref_section_s compatible */
	struct ref_section_s sec;

	/* reverse link index */
	uint32_t rv_link_idx_base;
	uint32_t reserved2;
};
_static_assert(sizeof(struct ref_section_intl_s) == 48);

/**
 * @struct ref_section_intl_half_s
 * @brief former half of the section_intl_s
 */
struct ref_section_intl_half_s {
	/* ref_section_s compatible */
	uint32_t reserved1[2];
	uint64_t reserved2;

	/* section table */
	uint32_t link_idx_base;
	uint32_t reserved3;
};
_static_assert(sizeof(struct ref_section_intl_half_s) == 24);

/**
 * @struct ref_hash_tuple_s
 */
struct ref_hash_tuple_s {
	uint64_t kmer;
	struct ref_gid_pos_s p;
};
typedef kvec_t(struct ref_hash_tuple_s) kvec_tuple_t;

/**
 * @struct ref_prec_s
 * @brief reference index precursor
 */
struct ref_prec_s {
	/* name -> section mapping */
	hmap_t *hmap;					/* name -> section_info hashmap */

	/* sequence container */
	kvec_t(uint64_t) seq;			/* 4bit packed seq */
	uint64_t seq_rem;

	/* link info container */
	kvec_t(struct ref_gid_pair_s) link;
	// uint64_t next_id;

	/* reserved */
	void *reserved1;
	uint64_t reserved2;
	void *reserved3;

	/* params */
	struct ref_params_s params;
};

/**
 * @struct ref_s
 */
struct ref_s {
	/* name -> section mapping */
	hmap_t *hmap;					/* name -> section_info hashmap */

	/* sequence container */
	kvec_t(uint64_t) seq;
	uint64_t seq_len;

	/* link info container */
	uint64_t mask;
	int64_t link_table_size;
	uint32_t *link_table;			/* fw and rv link array */

	int64_t *kmer_idx_table;
	int64_t kmer_table_size;
	struct ref_gid_pos_s *kmer_table;

	/* params */
	struct ref_params_s params;
};
_static_assert(sizeof(struct ref_prec_s) == sizeof(struct ref_s));
_static_assert_offset(struct ref_prec_s, hmap, struct ref_s, hmap, 0);
_static_assert_offset(struct ref_prec_s, seq, struct ref_s, seq, 0);
_static_assert_offset(struct ref_prec_s, params, struct ref_s, params, 0);

/**
 * @fn ref_prec_init
 */
ref_prec_t *ref_prec_init(
	ref_params_t const *params)
{
	/* check sanity of params */
	if(params == NULL) {
		return(NULL);
	}

	if(params->seed_length > 32) {
		return(NULL);
	}

	/* malloc mem */
	struct ref_prec_s *prec = (struct ref_prec_s *)malloc(sizeof(struct ref_prec_s));
	if(prec == NULL) {
		return(NULL);
	}

	/* fill default params and malloc mems */
	prec->hmap = hmap_init(1024, sizeof(struct ref_section_intl_s));

	kv_init(prec->seq);
	kv_push(prec->seq, 0);
	prec->seq_rem = 64;
	kv_init(prec->link);

	/* copy params */
	prec->params = *params;


	return((ref_prec_t *)prec);
}

/**
 * @fn ref_prec_clean
 */
void ref_prec_clean(
	ref_prec_t *_prec)
{
	struct ref_prec_s *prec = (struct ref_prec_s *)_prec;

	if(prec != NULL) {
		/* cleanup, cleanup... */
		hmap_clean(prec->hmap);
		kv_destroy(prec->seq);
		kv_destroy(prec->link);
	}
	return;
}

/**
 * @fn ref_encode_2bit
 * @brief mapping IUPAC amb. to 2bit encoding
 */
static _force_inline
uint8_t ref_encode_2bit(
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
 * @fn ref_encode_4bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t ref_encode_4bit(
	char c)
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
 * @fn ref_append_sequence
 * @brief append nucleotide sequence, must be null-terminated, accept IUPAC ambiguous encoding
 */
static _force_inline
struct ref_seq_interval_s ref_append_sequence(
	struct ref_prec_s *prec,
	char const *seq,
	int64_t len)
{
	/* load context onto registers */
	uint64_t rem = prec->seq_rem;		/* remaining length inside the array */
	uint64_t arr = kv_at(prec->seq, kv_size(prec->seq) - 1)<<rem;

	debug("rem(%llu), arr(%llx)", rem, arr);

	/* calc start pos */
	uint64_t base = (64 * kv_size(prec->seq) - rem) / 4;

	/* push until ptr reaches the end of the seq */
	for(int64_t i = 0; i < len; i++) {
		arr = (arr>>4) | ((uint64_t)ref_encode_4bit(seq[i])<<(64 - 4));
		debug("rem(%llu), arr(%llx)", rem, arr);
		if((rem -= 4) == 0) {
			kv_at(prec->seq, kv_size(prec->seq) - 1) = arr;
			kv_push(prec->seq, 0);
			arr = 0;
			rem = 64;
		}
	}

	/* write back contexts */
	prec->seq_rem = rem;
	kv_at(prec->seq, kv_size(prec->seq) - 1) = arr>>rem;
	debug("rem(%llu), arr(%llx)", rem, arr>>rem);

	/* calc tail */
	uint64_t tail = (64 * kv_size(prec->seq) - rem) / 4;
	debug("base(%llu), tail(%llu)", base, tail);

	return((struct ref_seq_interval_s){
		.base = base,
		.tail = tail
	});
}

/**
 * @fn ref_append_segment
 */
int ref_append_segment(
	ref_prec_t *_prec,
	char const *name,
	int32_t name_len,
	char const *seq,
	int64_t seq_len)
{
	struct ref_prec_s *prec = (struct ref_prec_s *)_prec;
	debug("append segment");

	/* add sequence at the tail of the seq buffer */
	struct ref_seq_interval_s iv = ref_append_sequence(prec, seq, seq_len);

	/* append the first section */
	uint64_t const max_sec_len = 0x80000000;
	uint64_t len = MIN2(iv.tail - iv.base, max_sec_len);

	uint32_t id = hmap_get_id(prec->hmap, name, name_len);
	struct ref_section_intl_s *sec =
		(struct ref_section_intl_s *)hmap_get_object(prec->hmap, id);

	/* store section info */
	sec->base_id = id;
	sec->sec = (struct ref_section_s){
		.id = id,
		.len = len,
		.base = iv.base
	};
	return(0);
}

/**
 * @fn ref_append_link
 */
int ref_append_link(
	ref_prec_t *_prec,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori)
{
	struct ref_prec_s *prec = (struct ref_prec_s *)_prec;
	debug("append link");

	/* get ids */
	uint32_t src_id = hmap_get_id(prec->hmap, src, src_len);
	uint32_t dst_id = hmap_get_id(prec->hmap, dst, dst_len);

	/* add forward link */
	kv_push(prec->link, ((struct ref_gid_pair_s){
		.from = _encode_id(src_id, src_ori),
		.to = _encode_id(dst_id, dst_ori)
	}));

	/* add reverse link */
	kv_push(prec->link, ((struct ref_gid_pair_s){
		.from = _encode_id(dst_id, _rev(dst_ori)),
		.to = _encode_id(src_id, _rev(src_ori))
	}));
	return(0);
}

/**
 * @fn ref_get_base
 */
static _force_inline
uint8_t ref_get_base(
	struct ref_prec_s const *prec,
	struct ref_section_s const *sec,
	uint32_t dir,
	uint32_t pos)
{
	uint64_t p = sec->base;
	uint64_t *ptr = kv_ptr(prec->seq);

	p += (dir == 0) ? pos : ((sec->len - 1) - pos);
	return((uint8_t)(0x0f & (ptr[p/16]>>((p & 0x0f) * 4))));
}

/**
 * @struct ref_pack_kmer_work_s
 */
struct ref_pack_kmer_work_s {
	uint64_t curr, cnt_arr;
};

/**
 * @fn ref_pack_kmer_sec_update
 */
static _force_inline
struct ref_pack_kmer_work_s ref_pack_kmer_sec_update(
	struct ref_prec_s *prec,
	struct ref_pack_kmer_work_s w,
	uint64_t *buf,
	uint8_t c)
{
	/* seed-related consts */
	uint64_t const seed_len = prec->params.seed_length;
	uint64_t const shift_len = 2 * (seed_len - 1);
	uint64_t const mask = (0x01<<(2 * seed_len)) - 1;

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
	w.cnt_arr = (w.cnt_arr<<2) | pcnt;

	/* branch */
	switch(3 - pcnt) {
		case 0: memcpy(&buf[2 * w.curr], buf, sizeof(uint64_t) * w.curr);
		case 1: memcpy(&buf[w.curr], buf, sizeof(uint64_t) * w.curr);
		/* fall through */
		default: break;
	}

	/* append to vector */
	for(int64_t j = 0; j < pcnt; j++) {
		for(int64_t k = 0; k < w.curr; k++) {
			// buf[j * w.curr + k] = mask &
				// ((buf[j * w.curr + k]<<2) | encode_2bit[c][j]);
			buf[j * w.curr + k] =
				(buf[j * w.curr + k]>>2) | (encode_2bit[c][j]<<shift_len);
			debug("%lld, %lld, %lld, %x, %llx",
				j, k, j * w.curr + k,
				encode_2bit[c][j], buf[j * w.curr + k]);
		}
	}

	/* update curr */
	w.curr *= pcnt;

	/* merge (shrink buffer) */
	uint64_t shrink_skip = 0x03 & (w.cnt_arr >> shift_len);
	debug("cnt_arr(%llx), w.curr(%llu), shrink_skip(%llu)",
		w.cnt_arr, w.curr, shrink_skip);
	if(shrink_skip > 1) {
		w.curr /= shrink_skip;
		for(int64_t j = 0; j < w.curr; j++) {
			buf[j] = buf[j * shrink_skip];
		}
	}

	return(w);
}

/**
 * @fn ref_pack_kmer_sec_push
 */
static _force_inline
void ref_pack_kmer_sec_push(
	struct ref_prec_s *prec,
	struct ref_pack_kmer_work_s w,
	uint64_t *buf,
	kvec_tuple_t *tuple_vec,
	uint32_t sec_id,
	uint32_t pos)
{
	/* section-related consts */
	uint32_t const fw_gid = _encode_id(sec_id, 0);
	// uint32_t const rv_gid = _encode_id(sec_id, 1);

	for(int64_t i = 0; i < w.curr; i++) {
		/* push forward */
		kv_push(*tuple_vec, ((struct ref_hash_tuple_s){
			.kmer = buf[i],
			.p.gid = fw_gid,
			.p.pos = pos
		}));
	}
	return;
}

/**
 * @fn ref_pack_kmer_sec
 */
static _force_inline
void ref_pack_kmer_sec(
	struct ref_prec_s *prec,
	kvec_tuple_t *tuple_vec,
	uint64_t *buf,
	uint32_t *link_idx,
	uint32_t sec_id)
{
	/* seed-related consts */
	uint64_t const seed_len = prec->params.seed_length;
	uint64_t const prefetch_len = seed_len - 1;

	/* working variables */
	struct ref_pack_kmer_work_s w = {
		.curr = 1,
		.cnt_arr = 0
	};

	/* section */
	struct ref_section_intl_s const *sec =
		(struct ref_section_intl_s const *)hmap_get_object(prec->hmap, sec_id);
	debug("seed_len(%llu), prefetch_len(%llu)",
		seed_len, prefetch_len);

	/* prefetch */
	buf[0] = 0;
	for(int64_t i = 0; i < prefetch_len; i++) {
		debug("i(%lld), c(%x)", i, ref_get_base(prec, &sec->sec, 0, i));
		w = ref_pack_kmer_sec_update(prec, w, buf,
			ref_get_base(prec, &sec->sec, 0, i));
	}

	/* body */
	for(int64_t i = prefetch_len; i < sec->sec.len; i++) {
		debug("i(%lld), c(%x)", i, ref_get_base(prec, &sec->sec, 0, i));
		w = ref_pack_kmer_sec_update(prec, w, buf,
			ref_get_base(prec, &sec->sec, 0, i));
		ref_pack_kmer_sec_push(prec, w, buf, tuple_vec,
			sec_id, i - prefetch_len);
	}

	/* tail */
	debug("(%u, %u)", sec->fw_link_idx_base, sec->rv_link_idx_base);
	for(int64_t j = sec->fw_link_idx_base; j < sec->rv_link_idx_base; j++) {
		uint32_t next_sec_id = _decode_id(link_idx[j]);
		uint32_t next_sec_dir = _decode_dir(link_idx[j]);

		/* copy curr and cnt_arr */
		struct ref_pack_kmer_work_s tw = w;

		struct ref_section_intl_s const *next_sec =
			(struct ref_section_intl_s const *)hmap_get_object(prec->hmap, next_sec_id);

		debug("link_idx(%lld), link to %u", j, next_sec_id);
		for(int64_t i = 0; i < prefetch_len; i++) {
			debug("i(%lld), c(%x)", i, ref_get_base(prec, &next_sec->sec, next_sec_dir, i));
			tw = ref_pack_kmer_sec_update(prec, tw, buf,
				ref_get_base(prec, &next_sec->sec, next_sec_dir, i));
			ref_pack_kmer_sec_push(prec, tw, buf, tuple_vec,
				sec_id, i + sec->sec.len - prefetch_len);
		}
	}
	return;
}

/**
 * @fn ref_pack_kmer
 * @brief pack kmer, sec must be sorted by local id
 */
static _force_inline
void ref_pack_kmer(
	struct ref_prec_s *prec,
	uint32_t tail_id,
	kvec_tuple_t *tuple_vec,
	uint32_t *link_idx)
{
	/* init buffer */
	uint64_t size = (uint64_t)pow(3.0, prec->params.seed_length * 0.5);
	uint64_t *buf = (uint64_t *)malloc(MAX2(size, 1024));

	debug("tail_id(%u)", tail_id);
	for(int64_t i = 0; i < tail_id; i++) {
		debug("pack_kmer id(%lld)", i);
		ref_pack_kmer_sec(prec, tuple_vec, buf, link_idx, i);
	}
	free(buf);
	return;
}

/**
 * @fn ref_build_kmer_idx_table
 */
static _force_inline
int64_t *ref_build_kmer_idx_table(
	struct ref_prec_s *prec,
	kvec_tuple_t *tuple_vec)
{
	kvec_t(int64_t) kmer_idx;
	kv_init(kmer_idx);

	uint64_t kmer_idx_size = 0x01 << (2 * prec->params.seed_length);
	kv_reserve(kmer_idx, kmer_idx_size);

	uint64_t prev_kmer = 0;
	kv_push(kmer_idx, prev_kmer);
	for(int64_t i = 0; i < kv_size(*tuple_vec); i++) {
		uint64_t kmer = kv_at(*tuple_vec, i).kmer;
		debug("i(%lld), kmer(%llx), id(%u), pos(%u), prev_kmer(%llx)",
			i, kmer, kv_at(*tuple_vec, i).p.gid, kv_at(*tuple_vec, i).p.pos, prev_kmer);

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
		kv_push(kmer_idx, kv_size(*tuple_vec));
	}
	return(kv_ptr(kmer_idx));
}

/**
 * @fn ref_build_kmer_shrink_table
 */
struct ref_kmer_table_s {
	int64_t *kmer_idx_table;
	int64_t kmer_table_size;
	struct ref_gid_pos_s *kmer_table;
};
static _force_inline
struct ref_kmer_table_s ref_build_kmer_shrink_table(
	struct ref_prec_s *prec,
	kvec_tuple_t *tuple_vec)
{
	struct ref_gid_pos_s *packed_pos = (struct ref_gid_pos_s *)kv_ptr(*tuple_vec);
	for(int64_t i = 0; i < kv_size(*tuple_vec); i++) {
		packed_pos[i] = kv_at(*tuple_vec, i).p;
	}

	return((struct ref_kmer_table_s){
		.kmer_table_size = kv_size(*tuple_vec),
		.kmer_table = packed_pos
	});
}

/**
 * @fn ref_build_kmer_table
 */
static _force_inline
struct ref_kmer_table_s ref_build_kmer_table(
	struct ref_prec_s *prec,
	uint32_t tail_id,
	uint32_t *link_idx)
{
	/* init kmer tuple vector */
	kvec_tuple_t v;
	kv_init(v);

	debug("build kmer table");

	/* pack kmer */
	ref_pack_kmer(prec, tail_id, &v, link_idx);

	/* sort vector by kmer */
	psort_half(
		kv_ptr(v),
		kv_size(v),
		sizeof(struct ref_hash_tuple_s),
		0);

	/* build index of kmer table */
	int64_t *kmer_idx_table = ref_build_kmer_idx_table(prec, &v);

	/* shrink kmer table */
	struct ref_kmer_table_s t = ref_build_kmer_shrink_table(prec, &v);
	t.kmer_idx_table = kmer_idx_table;
	return(t);
}

/**
 * @fn ref_build_link_idx_table
 */
static _force_inline
void ref_build_link_idx_table(
	struct ref_prec_s *prec,
	int64_t link_table_size,
	int64_t link_idx_table_size)
{
	/* sort by src, build src->dst mapping */
	debug("sort src->dst mapping, size(%llu)", kv_size(prec->link));
	psort_half(
		kv_ptr(prec->link),
		kv_size(prec->link),
		sizeof(struct ref_gid_pair_s),
		0);

	debug("forward list");
	for(int64_t i = 0; i < kv_size(prec->link); i++) {
		debug("(%u-%u, %u-%u)",
			_decode_id(kv_at(prec->link, i).from),
			_decode_dir(kv_at(prec->link, i).from),
			_decode_id(kv_at(prec->link, i).to),
			_decode_dir(kv_at(prec->link, i).to));
	}

	/* store forward info */
	uint32_t prev_gid = 0;
	struct ref_section_intl_half_s *sec_half =
		(struct ref_section_intl_half_s *)hmap_get_object(prec->hmap, 0);
	sec_half[prev_gid].link_idx_base = 0;
	for(int64_t i = 0; i < link_table_size; i++) {
		uint32_t gid = kv_at(prec->link, i).from;
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
	for(int64_t j = prev_gid + 1; j < link_idx_table_size + 1; j++) {
		debug("fill gaps j(%lld)", j);
		sec_half[j].link_idx_base = link_table_size;
	}
	return;
}

/**
 * @fn ref_build_link_shrink_table
 */
struct ref_link_table_s {
	uint32_t *gid_arr;
	int64_t size;
};
static _force_inline
struct ref_link_table_s ref_build_link_shrink_table(
	struct ref_prec_s *prec,
	int64_t link_table_size)
{
	/* pack */
	uint32_t *packed_gid = (uint32_t *)kv_ptr(prec->link);
	for(int64_t i = 0; i < link_table_size; i++) {
		packed_gid[i] = kv_at(prec->link, i).to;
	}
	for(int64_t i = 0; i < link_table_size; i++) {
		debug("%lld, %u", i, packed_gid[i]);
	}

	kv_resize(prec->link, link_table_size);
	return((struct ref_link_table_s){
		.gid_arr = (uint32_t *)kv_ptr(prec->link),
		.size = link_table_size
	});
}

/**
 * @fn ref_build_link_table
 */
static _force_inline
struct ref_link_table_s ref_build_link_table(
	struct ref_prec_s *prec,
	uint32_t tail_id)
{
	/* build id -> index on link table mapping */
	int64_t link_idx_table_size = 2 * tail_id;
	int64_t link_table_size = kv_size(prec->link);
	debug("build link_table, gid_size(%llu), size(%llu)",
		link_idx_table_size, link_table_size);

	/* build link_idx_table */
	ref_build_link_idx_table(prec, link_table_size, link_idx_table_size);

	/* shrink link_table */
	return(ref_build_link_shrink_table(prec, link_table_size));
}

/**
 * @fn ref_build_index_add_tail_sentinel
 */
static _force_inline
uint32_t ref_build_index_add_tail_sentinel(
	struct ref_prec_s *prec)
{
	char const *template = "tail_sentinel_";
	int64_t len = strlen(template);

	char buf[256];
	strcpy(buf, template);

	uint32_t id = (uint32_t)-1;
	uint32_t tail_id = hmap_get_count(prec->hmap);
	do {
		buf[len++] = '0';
		buf[len] = '\0';

		/* push sentinel to section array */
		id = hmap_get_id(prec->hmap, buf, len);
	} while(id != tail_id && len < 256);
	return(tail_id);
}

/**
 * @fn ref_build_index
 */
ref_t *ref_build_index(
	ref_prec_t *_prec)
{
	struct ref_prec_s *prec = (struct ref_prec_s *)_prec;
	if(prec == NULL) {
		goto _ref_build_index_error_handler;
	}

	dump(kv_ptr(prec->seq), 28 / 2);

	/* push tail sentinel */
	uint32_t tail_id = ref_build_index_add_tail_sentinel(prec);
	struct ref_section_intl_s *tail_sec =
		(struct ref_section_intl_s *)hmap_get_object(prec->hmap, tail_id);
	tail_sec->base_id = tail_id;
	tail_sec->sec = (struct ref_section_s){
		.id = tail_id,
		.len = 0,
		.base = 0
	};

	/* build link array */
	struct ref_link_table_s l = ref_build_link_table(prec, tail_id);
	if(l.gid_arr == NULL) {
		debug("gid_arr(%p)", l.gid_arr);
		goto _ref_build_index_error_handler;
	}

	/* build kmer array */
	struct ref_kmer_table_s k = ref_build_kmer_table(prec, tail_id, l.gid_arr);
	if(k.kmer_idx_table == NULL || k.kmer_table == NULL) {
		debug("kmer_idx_table(%p), kmer_table(%p)", k.kmer_idx_table, k.kmer_table);
		goto _ref_build_index_error_handler;
	}

	/* cast prec to ref */
	ref_t *ref = (ref_t *)prec;

	/* store misc */
	ref->mask = (0x01<<(2 * prec->params.seed_length)) - 1;
	ref->seq_len = (64 * kv_size(prec->seq) - prec->seq_rem) / 4;

	/* store link info */
	ref->link_table = l.gid_arr;
	ref->link_table_size = l.size;

	/* store kmer info */
	ref->kmer_idx_table = k.kmer_idx_table;
	ref->kmer_table_size = k.kmer_table_size;
	ref->kmer_table = k.kmer_table;
	return(ref);

_ref_build_index_error_handler:;
	if(prec != NULL) {
		hmap_clean(prec->hmap);
		if(kv_ptr(prec->seq) != NULL) { kv_destroy(prec->seq); }

		/* l.gid_arr and prec.link shares the memory */
		if(l.gid_arr != NULL) { free(l.gid_arr); }
	}
	if(k.kmer_idx_table != NULL) { free(k.kmer_idx_table); }
	if(k.kmer_table != NULL) { free(k.kmer_table); }
	return(NULL);
}

/**
 * @fn ref_clean
 */
void ref_clean(
	ref_t *_ref)
{
	struct ref_s *ref = (struct ref_s *)_ref;

	if(ref != NULL) {
		/* cleanup, cleanup... */
		hmap_clean(ref->hmap);
		if(kv_ptr(ref->seq) != NULL) { kv_destroy(ref->seq); }
		if(ref->link_table != NULL) { free(ref->link_table); }
		if(ref->kmer_idx_table != NULL) { free(ref->kmer_idx_table); }
		if(ref->kmer_table != NULL) { free(ref->kmer_table); }
	}
	return;
}

/**
 * @fn ref_dump_index
 */
int ref_dump_index(
	ref_t const *ref,
	zf_t *outfp)
{
	return(0);
}

/**
 * @fn ref_load_index
 */
ref_t *ref_load_index(
	zf_t *infp)
{
	return(NULL);
}

/**
 * @fn ref_get_section
 */
struct ref_section_s const *ref_get_section(
	ref_t const *_ref,
	uint32_t id)
{
	struct ref_s *ref = (struct ref_s *)_ref;

	struct ref_section_intl_s *sec =
		(struct ref_section_intl_s *)hmap_get_object(ref->hmap, id);
	return((struct ref_section_s const *)&sec->sec);
}

/**
 * @fn ref_get_name
 */
struct ref_str_s ref_get_name(
	ref_t const *_ref,
	uint32_t id)
{
	struct ref_s *ref = (struct ref_s *)_ref;
	struct hmap_key_s key = hmap_get_key(ref->hmap, id);
	return((struct ref_str_s){
		.str = key.str,
		.len = key.len
	});
}

/**
 * @fn ref_get_ptr
 */
uint8_t const *ref_get_ptr(
	ref_t const *_ref)
{
	struct ref_s const *ref = (struct ref_s const *)_ref;
	return((uint8_t const *)kv_ptr(ref->seq));
}

/**
 * @fn ref_get_total_len
 */
int64_t ref_get_total_len(
	ref_t const *_ref)
{
	struct ref_s const *ref = (struct ref_s const *)_ref;
	return(ref->seq_len);
}

/**
 * @fn ref_match_2bitpacked
 */
struct ref_match_res_s ref_match_2bitpacked(
	ref_t const *_ref,
	uint64_t seq)
{
	struct ref_s const *ref = (struct ref_s const *)_ref;
	seq &= ref->mask;
	int64_t base = ref->kmer_idx_table[seq];
	int64_t tail = ref->kmer_idx_table[seq + 1];

	debug("seq(%llx), mask(%llx), base(%lld), tail(%lld)",
		seq, ref->mask, base, tail);
	
	struct ref_section_intl_s *sec =
		(struct ref_section_intl_s *)hmap_get_object(ref->hmap, _decode_id(ref->kmer_table[base].gid));
	debug("id(%u), base(%llu), len(%u)",
		sec->sec.id, sec->sec.base, sec->sec.len);
	return((struct ref_match_res_s){
		.ptr = &ref->kmer_table[base],
		.len = tail - base
	});
}

/**
 * @fn ref_match
 * @brief seq length must be equal to seed_length.
 */
struct ref_match_res_s ref_match(
	ref_t const *_ref,
	char const *seq)
{
	struct ref_s const *ref = (struct ref_s const *)_ref;
	int64_t const seed_len = ref->params.seed_length;
	int64_t const shift_len = 2 * (seed_len - 1);

	uint64_t packed_seq = 0;
	for(int64_t i = 0; i < seed_len; i++) {
		packed_seq = (packed_seq>>2) | (ref_encode_2bit(seq[i])<<shift_len);
	}
	return(ref_match_2bitpacked((ref_t const *)ref, packed_seq));
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

/* make prec context */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(
		.seed_length = 3));

	assert(prec != NULL);

	ref_prec_clean(prec);
}

/* add segment */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(.seed_length = 3));

	int ret = ref_append_segment(prec, _str("sec0"), _str("AARA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = ref_append_segment(prec, _str("sec1"), _str("MAAA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = ref_append_link(prec, _str("sec0"), 0, _str("sec1"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = ref_append_link(prec, _str("sec1"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = ref_append_segment(prec, _str("sec2"), _str("ACGT"));
	assert(ret == 0, "ret(%d)", ret);

	ret = ref_append_link(prec, _str("sec0"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ref_prec_clean(prec);
}

/* build index */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(.seed_length = 3));

	/* append */
	ref_append_segment(prec, _str("sec0"), _str("GGRA"));
	ref_append_segment(prec, _str("sec1"), _str("MGGG"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec1"), 0);
	ref_append_link(prec, _str("sec1"), 0, _str("sec2"), 0);
	ref_append_segment(prec, _str("sec2"), _str("ACVVGTGT"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */
	ref_t *ref = ref_build_index(prec);
	assert(ref != NULL, "ref(%p)", ref);

	/* pointer to seq */
	uint8_t const *ptr = ref_get_ptr(ref);
	assert(ptr != NULL, "ptr(%p)", ptr);

	/* total len */
	assert(ref_get_total_len(ref) == 16, "len(%lld)", ref_get_total_len(ref));

	ref_clean(ref);
}

/* get_section */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(.seed_length = 3));
	ref_append_segment(prec, _str("sec0"), _str("GGRA"));
	ref_append_segment(prec, _str("sec1"), _str("MGGG"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec1"), 0);
	ref_append_link(prec, _str("sec1"), 0, _str("sec2"), 0);
	ref_append_segment(prec, _str("sec2"), _str("ACVVGTGT"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec2"), 0);
	ref_t *ref = ref_build_index(prec);

	/* section id is given in ascending order from 0 */
	assert(ref_get_section(ref, 0) != NULL, "%p", ref_get_section(ref, 0));
	assert(ref_get_section(ref, 0)->id == 0, "id(%u)", ref_get_section(ref, 0)->id);
	assert(ref_get_section(ref, 0)->len == 4, "len(%u)", ref_get_section(ref, 0)->len);
	assert(ref_get_section(ref, 0)->base == 0, "base(%llu)", ref_get_section(ref, 0)->base);

	/* section 1 */
	assert(ref_get_section(ref, 1) != NULL, "%p", ref_get_section(ref, 1));
	assert(ref_get_section(ref, 1)->id == 1, "id(%u)", ref_get_section(ref, 1)->id);
	assert(ref_get_section(ref, 1)->len == 4, "len(%u)", ref_get_section(ref, 1)->len);
	assert(ref_get_section(ref, 1)->base == 4, "base(%llu)", ref_get_section(ref, 1)->base);

	/* section 2 */
	assert(ref_get_section(ref, 2) != NULL, "%p", ref_get_section(ref, 2));
	assert(ref_get_section(ref, 2)->id == 2, "id(%u)", ref_get_section(ref, 2)->id);
	assert(ref_get_section(ref, 2)->len == 8, "len(%u)", ref_get_section(ref, 2)->len);
	assert(ref_get_section(ref, 2)->base == 8, "base(%llu)", ref_get_section(ref, 2)->base);

	ref_clean(ref);
}

/* get_name */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(.seed_length = 3));
	ref_append_segment(prec, _str("sec0"), _str("GGRA"));
	ref_append_segment(prec, _str("sec1"), _str("MGGG"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec1"), 0);
	ref_append_link(prec, _str("sec1"), 0, _str("sec2"), 0);
	ref_append_segment(prec, _str("sec2"), _str("ACVVGTGT"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec2"), 0);
	ref_t *ref = ref_build_index(prec);

	/* section id is given in ascending order from 0 */
	assert(ref_get_name(ref, 0).len == 4, "%d", ref_get_name(ref, 0).len);
	assert(strcmp(ref_get_name(ref, 0).str, "sec0") == 0, "%s", ref_get_name(ref, 0).str);

	/* section 1 */
	assert(ref_get_name(ref, 1).len == 4, "%d", ref_get_name(ref, 1).len);
	assert(strcmp(ref_get_name(ref, 1).str, "sec1") == 0, "%s", ref_get_name(ref, 1).str);

	/* section 2 */
	assert(ref_get_name(ref, 2).len == 4, "%d", ref_get_name(ref, 2).len);
	assert(strcmp(ref_get_name(ref, 2).str, "sec2") == 0, "%s", ref_get_name(ref, 2).str);

	ref_clean(ref);
}

/* match */
unittest()
{
	ref_prec_t *prec = ref_prec_init(REF_PARAMS(.seed_length = 3));
	ref_append_segment(prec, _str("sec0"), _str("GGRA"));
	ref_append_segment(prec, _str("sec1"), _str("MGGG"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec1"), 0);
	ref_append_link(prec, _str("sec1"), 0, _str("sec2"), 0);
	ref_append_segment(prec, _str("sec2"), _str("ACVVGTGT"));
	ref_append_link(prec, _str("sec0"), 0, _str("sec2"), 0);
	ref_t *ref = ref_build_index(prec);

	/* without ambiguous bases */
	struct ref_match_res_s r = ref_match(ref, "GTG");
	assert(r.ptr != NULL, "%p", r.ptr);
	assert(r.len == 1, "%lld", r.len);

	/* check pos */
	assert(r.ptr[0].pos == 4, "%u", r.ptr[0].pos);

	/* check section */
	struct ref_section_s const *sec = ref_get_section(ref, ref_id(r.ptr[0].gid));
	assert(sec->id == 2, "id(%u)", sec->id);
	assert(sec->len == 8, "len(%u)", sec->len);
	assert(sec->base == 8, "base(%llu)", sec->base);


	/* with ambiguous bases */
	r = ref_match(ref, "GGG");
	assert(r.ptr != NULL, "%p", r.ptr);
	assert(r.len == 2, "%lld", r.len);

	/* check pos */
	assert(r.ptr[0].pos == 0, "%u", r.ptr[0].pos);

	/* check section */
	sec = ref_get_section(ref, ref_id(r.ptr[0].gid));
	assert(sec->id == 0, "id(%u)", sec->id);
	assert(sec->len == 4, "len(%u)", sec->len);
	assert(sec->base == 0, "base(%llu)", sec->base);

	/* check pos */
	assert(r.ptr[1].pos == 1, "%u", r.ptr[1].pos);

	/* check section */
	sec = ref_get_section(ref, ref_id(r.ptr[1].gid));
	assert(sec->id == 1, "id(%u)", sec->id);
	assert(sec->len == 4, "len(%u)", sec->len);
	assert(sec->base == 4, "base(%llu)", sec->base);
}

#endif

/**
 * end of ref.c
 */
