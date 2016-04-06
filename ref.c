
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
 * MurmurHash3 string hashing function,
 * extracted from https://github.com/aappleby/smhasher
 * modified to make the functions static and to return
 * hashed value directly.
 */
static _force_inline
uint32_t rotl32(uint32_t x,int8_t r)
{
	return((x << r) | (x >> (32 - r)));
}

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

static _force_inline
uint32_t getblock32(const uint32_t *p, int i)
{
	return(p[i]);
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static _force_inline
uint32_t fmix32(uint32_t h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return(h);
}

//-----------------------------------------------------------------------------
static _force_inline
uint32_t MurmurHash3_x86_32(
	const void *key,
	int32_t len,
	uint32_t seed)
{
	const uint8_t * data = (const uint8_t*)key;
	const int nblocks = len / 4;

	uint32_t h1 = seed;

	const uint32_t c1 = 0xcc9e2d51;
	const uint32_t c2 = 0x1b873593;

	//----------
	// body

	const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

	for(int i = -nblocks; i; i++)
	{
		uint32_t k1 = getblock32(blocks,i);

		k1 *= c1;
		k1 = rotl32(k1,15);
		k1 *= c2;

		h1 ^= k1;
		h1 = rotl32(h1,13); 
		h1 = h1*5+0xe6546b64;
	}

	//----------
	// tail

	const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

	uint32_t k1 = 0;

	switch(len & 3)
	{
	case 3: k1 ^= tail[2] << 16;
	case 2: k1 ^= tail[1] << 8;
	case 1: k1 ^= tail[0];
		k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
	};

	//----------
	// finalization

	h1 ^= len;

	h1 = fmix32(h1);

	return(h1);
}
/* end of MurmurHash3.cpp */

/**
 * MurmurHash3 wrapper functions
 */
static _force_inline
uint32_t hash_string(char const *str, int32_t len)
{
	return(MurmurHash3_x86_32(
		(void const *)str,
		(int)len,
		(uint32_t)0xcafebabe));
}
static _force_inline
uint32_t hash_rehash(uint32_t val)
{
	return(MurmurHash3_x86_32(
		(void const *)&val,
		4,
		(uint32_t)0xcafebabe));
}

/**
 * @struct hmap_s
 */
struct hmap_s {
	uint32_t mask;
	uint32_t reserved;
	uint32_t *t;
};

/**
 * hashmap manipulating functions
 */
static _force_inline
struct hmap_s *hmap_init(int64_t size)
{
	struct hmap_s *hmap = malloc(
		sizeof(struct hmap_s) + sizeof(uint32_t) * size);
	*hmap = (struct hmap_s){
		.mask = size - 1,
		.t = (uint32_t *)(hmap + 1)
	};
	memset(hmap->t, 0xff, sizeof(uint32_t) * size);
	return(hmap);
}
static _force_inline
void hmap_clean(struct hmap_s *hash)
{
	free((void *)hash);
	return;
}
static _force_inline
void hmap_add(struct hmap_s *hash, uint32_t hash_val, uint32_t val)
{
	debug("h(%x, %x), v(%u)", hash_val, hash_val & hash->mask, val);
	hash->t[hash->mask & hash_val] = val;
	return;
}
static _force_inline
uint32_t hmap_get(struct hmap_s *hash, uint32_t hash_val)
{
	debug("h(%x, %x)", hash_val, hash_val & hash->mask);
	return(hash->t[hash->mask & hash_val]);
}
static _force_inline
void hmap_resize(struct hmap_s *hash, uint64_t size)
{
	return;
}

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
	/* ref_section_s compatible */
	uint32_t id;
	uint32_t len;
	uint64_t base;

	/* internal fields */
	uint64_t name_ofs;		/* offset in the name vector */
	int32_t name_len;

	/* splitted section link (used to get original name) */
	uint32_t base_id;

	/* section table */
	uint32_t fw_link_idx_base;
	uint32_t rv_link_idx_base;

	uint64_t reserved;
};
_static_assert(sizeof(struct ref_section_intl_s) == 48);

/**
 * @struct ref_gid_pos_s
 */
struct ref_gid_pos_s {
	uint32_t gid;
	uint32_t pos;
};

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
	/* name to id mapping */
	kvec_t(char) name;						/* name vector */
	struct hmap_s *hmap;					/* name to id hmap */

	/* sequence container */
	kvec_t(uint64_t) seq;					/* 4bit packed seq */
	uint64_t seq_rem;

	/* section info container */
	kvec_t(struct ref_section_intl_s) sec;
	kvec_t(struct ref_gid_pair_s) fw_link;
	kvec_t(struct ref_gid_pair_s) rv_link;
	uint64_t next_id;

	/* params */
	struct ref_params_s params;
};

/**
 * @struct ref_s
 */
struct ref_s {
	/* name to id mapping */
	kvec_t(char) name;						/* name vector */
	struct hmap_s *hmap;					/* name to id hmap */

	/* sequence container */
	kvec_t(uint64_t) seq;
	uint64_t seq_len;

	/* section info container */
	kvec_t(struct ref_section_intl_s) sec;
	uint64_t reserved1[2];
	uint32_t *link_idx;						/* fw and rv link array */
	uint64_t reserved2[3];
	uint64_t link_size;

	/* params */
	struct ref_params_s params;
};
_static_assert(sizeof(struct ref_prec_s) == sizeof(struct ref_s));
_static_assert_offset(struct ref_prec_s, name, struct ref_s, name, 0);
_static_assert_offset(struct ref_prec_s, hmap, struct ref_s, hmap, 0);
_static_assert_offset(struct ref_prec_s, seq, struct ref_s, seq, 0);
_static_assert_offset(struct ref_prec_s, seq_rem, struct ref_s, seq_len, 0);
_static_assert_offset(struct ref_prec_s, sec, struct ref_s, sec, 0);
_static_assert_offset(struct ref_prec_s, fw_link.a, struct ref_s, link_idx, 0);
_static_assert_offset(struct ref_prec_s, next_id, struct ref_s, link_size, 0);
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
	kv_init(prec->name);
	prec->hmap = hmap_init(1024);

	kv_init(prec->seq);
	kv_push(prec->seq, 0);
	prec->seq_rem = 64;
	prec->next_id = 0;

	kv_init(prec->sec);
	kv_init(prec->fw_link);
	kv_init(prec->rv_link);

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
		kv_destroy(prec->name);
		hmap_clean(prec->hmap);

		kv_destroy(prec->seq);

		kv_destroy(prec->sec);
		kv_destroy(prec->fw_link);
		kv_destroy(prec->rv_link);
	}
	return;
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
 * @fn ref_append_name
 */
static _force_inline
uint64_t ref_append_name(
	struct ref_prec_s *prec,
	char const *name,
	int32_t len)
{
	uint64_t ofs = kv_size(prec->name);

	kv_reserve(prec->name, kv_size(prec->name) + len + 1);
	for(int64_t i = 0; i < len; i++) {
		kv_push(prec->name, name[i]);
	}
	kv_push(prec->name, '\0');
	return(ofs);
}

/**
 * @struct ref_hmap_res_s
 */
struct ref_hmap_res_s {
	uint32_t h;
	uint32_t id;
};

/**
 * @fn ref_hash_find_name
 */
static _force_inline
struct ref_hmap_res_s ref_hash_find_name(
	struct ref_prec_s *prec,
	char const *name,
	int32_t name_len)
{
	uint32_t tmp_id = -1, base_id = -1;
	uint32_t h = hash_string(name, name_len);
	while((tmp_id = hmap_get(prec->hmap, h)) != -1) {
		/* compare string */
		if(kv_at(prec->sec, tmp_id).name_len == name_len
		&& strcmp(kv_ptr(prec->name) + kv_at(prec->sec, tmp_id).name_ofs, name) == 0) {
			/* match with existing string in the section array */
			base_id = tmp_id; break;
		}

		/* not matched, rehash */
		h = hash_rehash(h);
	}
	debug("id(%u), hash(%x), str(%s)", base_id, h, name);
	return((struct ref_hmap_res_s){
		.h = h,
		.id = base_id
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
	/* calc section number */
	uint64_t const max_sec_len = 0x80000000;

	/* put (name -> id) to hash */
	struct ref_hmap_res_s h = ref_hash_find_name(prec, name, name_len);

	/* add sequence at the tail of the seq buffer */
	struct ref_seq_interval_s iv = ref_append_sequence(prec, seq, seq_len);

	/* append the first section */ {
		uint64_t len = MIN2(iv.tail - iv.base, max_sec_len);
		if(h.id == -1) {
			/* acquire a new id */
			h.id = prec->next_id++;
			kv_push(prec->sec, ((struct ref_section_intl_s){
				.id = h.id,
				.len = len,
				.base = iv.base,
				.name_ofs = ref_append_name(prec, name, name_len),
				.name_len = name_len,
				.base_id = h.id
			}));

			debug("reserve new key(%u)", h.id);

			/* add key to hash */
			hmap_add(prec->hmap, h.h, h.id);
		} else {
			/* modify existing section */
			debug("base_id(%u), sec[base_id].id(%u)", h.id, kv_at(prec->sec, h.id).id);
			debug("name(%s, %s)", name, kv_ptr(prec->name) + kv_at(prec->sec, h.id).name_ofs);
			kv_at(prec->sec, h.id).len = len;
			kv_at(prec->sec, h.id).base = iv.base;
		}
		iv.base += len;
	}

	/* append remaining sections */
	while(iv.base < iv.tail) {
		uint32_t id = prec->next_id++;

		/* append section info */
		uint64_t len = MIN2(iv.tail - iv.base, max_sec_len);
		kv_push(prec->sec, ((struct ref_section_intl_s){
			.id = id,
			.len = len,
			.base = iv.base,
			.name_ofs = ref_append_name(prec, "", 0),
			.name_len = 0,
			.base_id = h.id
		}));
		iv.base += len;
	}
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

	/* put src info */
	struct ref_hmap_res_s s = ref_hash_find_name(prec, src, src_len);
	if(s.id == -1) {
		/* acquire a new id for src */
		s.id = prec->next_id++;
		kv_push(prec->sec, ((struct ref_section_intl_s){
			.id = s.id,
			.name_ofs = ref_append_name(prec, src, src_len),
			.name_len = src_len,
			.base_id = s.id
		}));

		/* add key to hash */
		hmap_add(prec->hmap, s.h, s.id);
	}

	/* put dst info */
	struct ref_hmap_res_s d = ref_hash_find_name(prec, dst, dst_len);
	if(d.id == -1) {
		/* acquire a new id for dst */
		d.id = prec->next_id++;
		kv_push(prec->sec, ((struct ref_section_intl_s){
			.id = d.id,
			.name_ofs = ref_append_name(prec, dst, dst_len),
			.name_len = dst_len,
			.base_id = d.id
		}));

		/* add key to hash */
		hmap_add(prec->hmap, d.h, d.id);
	}

	/* add forward link */
	kv_push(prec->fw_link, ((struct ref_gid_pair_s){
		.from = _encode_id(s.id, src_ori),
		.to = _encode_id(d.id, dst_ori)
	}));

	/* add reverse link */
	kv_push(prec->rv_link, ((struct ref_gid_pair_s){
		.from = _encode_id(d.id, _rev(dst_ori)),
		.to = _encode_id(s.id, _rev(src_ori))
	}));
	return(0);
}

/**
 * @fn ref_get_base
 */
static _force_inline
uint8_t ref_get_base(
	struct ref_prec_s *prec,
	uint32_t id,
	uint32_t dir,
	uint32_t pos)
{
	uint64_t p = kv_at(prec->sec, id).base;
	uint64_t *ptr = kv_ptr(prec->seq);

	p += (dir == 0)
		? pos
		: ((kv_at(prec->sec, id).len - 1) - pos);
	return(0x0f & (ptr[p/16]>>((p & 0x0f) * 4)));
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
			buf[j * w.curr + k] = mask &
				((buf[j * w.curr + k]<<2) | encode_2bit[c][j]);
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
	uint64_t const shift_len = 2 * (seed_len - 1);

	/* working variables */
	struct ref_pack_kmer_work_s w = {
		.curr = 1,
		.cnt_arr = 0
	};

	debug("sec.len(%u), seed_len(%llu), prefetch_len(%llu), shift_len(%llu)",
		kv_at(prec->sec, sec_id).len,
		seed_len, prefetch_len, shift_len);

	/* prefetch */
	for(int64_t i = 0; i < prefetch_len; i++) {
		w = ref_pack_kmer_sec_update(prec, w, buf,
			ref_get_base(prec, sec_id, 0, i));
	}

	/* body */
	for(int64_t i = prefetch_len; i < kv_at(prec->sec, sec_id).len; i++) {
		w = ref_pack_kmer_sec_update(prec, w, buf,
			ref_get_base(prec, sec_id, 0, i));
		ref_pack_kmer_sec_push(prec, w, buf, tuple_vec,
			sec_id, i - prefetch_len);
	}

	/* tail */
	for(int64_t j = kv_at(prec->sec, sec_id).fw_link_idx_base;
		j < kv_at(prec->sec, sec_id + 1).fw_link_idx_base;
		j++) {
		uint32_t next_sec_id = _decode_id(link_idx[j]);
		uint32_t next_sec_dir = _decode_dir(link_idx[j]);

		/* copy curr and cnt_arr */
		struct ref_pack_kmer_work_s tw = w;

		debug("link_idx(%lld), link to %u", j, next_sec_id);
		for(int64_t i = 0; i < prefetch_len; i++) {
			tw = ref_pack_kmer_sec_update(prec, tw, buf,
				ref_get_base(prec, next_sec_id, next_sec_dir, i));
			ref_pack_kmer_sec_push(prec, tw, buf, tuple_vec,
				sec_id, i + kv_at(prec->sec, sec_id).len - prefetch_len);
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
	kvec_tuple_t *tuple_vec,
	uint32_t *link_idx)
{
	/* init buffer */
	uint64_t size = (uint64_t)pow(3.0, prec->params.seed_length * 0.5);
	uint64_t *buf = (uint64_t *)malloc(size);

	for(int64_t i = 0; i < prec->next_id; i++) {
		debug("pack_kmer id(%lld), base(%llu), len(%u)",
			i, kv_at(prec->sec, i).base, kv_at(prec->sec, i).len);
		ref_pack_kmer_sec(prec, tuple_vec, buf, link_idx, i);
	}
	free(buf);
	return;
}

/**
 * @fn ref_build_kmer_array
 */
struct ref_kmer_idx_s {
	struct ref_gid_pos_s *p_arr;
	int64_t size;
};
static _force_inline
struct ref_kmer_idx_s ref_build_kmer_array(
	struct ref_prec_s *prec,
	uint32_t *link_idx)
{
	/* init kmer tuple vector */
	kvec_tuple_t v;
	kv_init(v);

	/* pack kmer */
	ref_pack_kmer(prec, &v, link_idx);

	/* sort vector by kmer */
	psort_half(
		kv_ptr(v),
		kv_size(v),
		sizeof(struct ref_hash_tuple_s),
		0);

	for(int64_t i = 0; i < kv_size(v); i++) {
		debug("kmer(%llx), id(%u, %u), pos(%u)",
			kv_at(v, i).kmer,
			_decode_id(kv_at(v, i).p.gid),
			_decode_dir(kv_at(v, i).p.gid),
			kv_at(v, i).p.pos);
	}

	return((struct ref_kmer_idx_s){
		.p_arr = (struct ref_gid_pos_s *)kv_ptr(v),
		.size = kv_size(v)
	});
}

/**
 * @fn ref_build_link_array
 */
struct ref_link_idx_s {
	uint32_t *gid_arr;
	int64_t size;
};
static _force_inline
struct ref_link_idx_s ref_build_link_array(
	struct ref_prec_s *prec)
{
	/* sort by src, build src->dst mapping */
	debug("sort src->dst mapping, size(%llu)", kv_size(prec->fw_link));
	psort_half(
		kv_ptr(prec->fw_link),
		kv_size(prec->fw_link),
		sizeof(struct ref_gid_pair_s),
		0);

	debug("forward list");
	for(int64_t i = 0; i < kv_size(prec->fw_link); i++) {
		debug("(%u-%u, %u-%u)",
			_decode_id(kv_at(prec->fw_link, i).from),
			_decode_dir(kv_at(prec->fw_link, i).from),
			_decode_id(kv_at(prec->fw_link, i).to),
			_decode_dir(kv_at(prec->fw_link, i).to));
	}

	/* sort by dst, build dst->src (reverse) mapping */
	debug("sort dst->src mapping, size(%llu)", kv_size(prec->rv_link));
	psort_half(
		kv_ptr(prec->rv_link),
		kv_size(prec->rv_link),
		sizeof(struct ref_gid_pair_s),
		0);

	debug("reverse list");
	for(int64_t i = 0; i < kv_size(prec->rv_link); i++) {
		debug("(%u-%u, %u-%u)",
			_decode_id(kv_at(prec->rv_link, i).from),
			_decode_dir(kv_at(prec->rv_link, i).from),
			_decode_id(kv_at(prec->rv_link, i).to),
			_decode_dir(kv_at(prec->rv_link, i).to));
	}

	/* build id -> index on link table mapping */
	uint32_t tail_id = prec->next_id;
	int64_t sec_arr_size = tail_id + 1;
	int64_t link_arr_size = kv_size(prec->fw_link);
	debug("build id->index table, size(%llu)", 2 * sec_arr_size);

	/* forward */ {
		uint32_t prev_id = 0;
		kv_at(prec->sec, prev_id).fw_link_idx_base = 0;
		for(int64_t i = 0; i < link_arr_size; i++) {
			uint32_t id = _decode_id(kv_at(prec->fw_link, i).from);
			debug("id(%u), prev_id(%u)", id, prev_id);

			if(prev_id == id) { continue; }

			debug("fw, index table for id(%u) ends at %lld, next_id(%u)", prev_id, i, id);

			/* sequence update detected */
			for(int64_t j = prev_id; j < id; j++) {
				debug("fill gaps j(%lld)", j);
				kv_at(prec->sec, j + 1).fw_link_idx_base = i;
			}
			prev_id = id;
		}
		for(int64_t i = prev_id; i < tail_id; i++) {
			kv_at(prec->sec, i + 1).fw_link_idx_base = link_arr_size;
		}
	}

	/* reverse */ {
		uint32_t prev_id = 0;
		kv_at(prec->sec, prev_id).rv_link_idx_base = 0;
		for(int64_t i = 0; i < link_arr_size; i++) {
			uint32_t id = _decode_id(kv_at(prec->rv_link, i).from);
			debug("id(%u), prev_id(%u)", id, prev_id);

			if(prev_id == id) { continue; }

			debug("rv, index table for id(%u) ends at %lld, next_id(%u)", prev_id, i, id);

			/* sequence update detected */
			for(int64_t j = prev_id; j < id; j++) {
				debug("fill gaps j(%lld)", j);
				kv_at(prec->sec, j + 1).rv_link_idx_base = i;
			}
			prev_id = id;
		}
		for(int64_t i = prev_id; i < tail_id; i++) {
			kv_at(prec->sec, i + 1).rv_link_idx_base = link_arr_size;
		}
	}

	for(int64_t i = 0; i < sec_arr_size; i++) {
		debug("%lld, %u, %u", i,
			kv_at(prec->sec, i).fw_link_idx_base,
			kv_at(prec->sec, i).rv_link_idx_base);
	}

	/* pack */
	uint32_t *packed_gid = (uint32_t *)kv_ptr(prec->fw_link);
	for(int64_t i = 0; i < link_arr_size; i++) {
		packed_gid[i] = kv_at(prec->fw_link, i).to;
	}
	for(int64_t i = 0; i < link_arr_size; i++) {
		packed_gid[i + link_arr_size] = kv_at(prec->rv_link, i).to;
	}

	for(int64_t i = 0; i < 2 * link_arr_size; i++) {
		debug("%lld, %u", i, packed_gid[i]);
	}

	/* cleanup rv_link */
	kv_destroy(prec->rv_link);

	return((struct ref_link_idx_s){
		.gid_arr = packed_gid,
		.size = sec_arr_size
	});
}

/**
 * @fn ref_build_index
 */
ref_t *ref_build_index(
	ref_prec_t *_prec)
{
	struct ref_prec_s *prec = (struct ref_prec_s *)_prec;

	dump(kv_ptr(prec->seq), 28 / 2);

	/* push sentinel to section array */
	kv_push(prec->sec, ((struct ref_section_intl_s){
		.id = prec->next_id,
		.base = 0,
		.len = 0,
		.name_ofs = ref_append_name(prec, "", 0),
		.name_len = 0,
		.base_id = prec->next_id
	}));

	/* build link array */
	struct ref_link_idx_s l = ref_build_link_array(prec);

	/* build kmer array */
	struct ref_kmer_idx_s k = ref_build_kmer_array(prec, l.gid_arr);

	ref_t *ref = (ref_t *)prec;
	ref->link_idx = l.gid_arr;
	ref->link_size = l.size;

	return(ref);
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
		kv_destroy(ref->name);
		hmap_clean(ref->hmap);

		kv_destroy(ref->seq);

		kv_destroy(ref->sec);
		free(ref->link_idx);
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
 * unittests
 */
#ifdef TEST

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

	ref_clean(ref);
}

#endif

/**
 * end of ref.c
 */
