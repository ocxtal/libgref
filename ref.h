
/**
 * @file ref.h
 *
 * @brief a header of ref.c
 *
 * @author Hajime Suzuki
 * @date 2015/6/23
 * @license Apache v2.
 *
 * @detail
 * Ref is a hash-based sequence indexer and exact matcher. 
 */

#ifndef _REF_H_INCLUDED
#define _REF_H_INCLUDED

#include <stdint.h>

/**
 * @enum ref_error
 * @brief error flags
 */
enum ref_error {
	/** error codes */
	REF_SUCCESS 				= 0,
	REF_ERROR 					= 1,
	REF_ERROR_INVALID_CONTEXT	= 2,
	REF_ERROR_INVALID_ARGS 		= 3,
	REF_ERROR_OVERWRITE			= 4,
	REF_ERROR_FILE_NOT_FOUND	= 5,
	REF_ERROR_BROKEN_FILE		= 6,

	/** return values */
	REF_INDEX_VALID				= 0,
	REF_INDEX_INVALID			= -1
};

/**
 * @enum ref_format
 */
enum ref_format_flags {
	REF_ASCII				= 1,
	REF_2BIT				= 2,
	REF_2BIT8PACKED			= 3
};

/**
 * @type ref_prec_t
 */
typedef struct ref_prec_s ref_prec_t;

/**
 * @type ref_t
 * @brief an alias to struct ref_s
 */
typedef struct ref_s ref_t;

/**
 * @struct ref_params_s
 */
struct ref_params_s {
	int32_t seed_length;
	uint32_t reserved;
};
typedef struct ref_params_s ref_params_t;
#define REF_PARAMS(...)			( &((struct ref_params_s const) { __VA_ARGS__ }) )

/**
 * @struct ref_pos_s
 * @brief sizeof(struct ref_pos_s) == 8
 */
struct ref_pos_s {
	int32_t id;
	int32_t pos;
};
typedef struct ref_pos_s ref_pos_t;

/**
 * @struct ref_section_s
 * @brief has equivalent fields to struct sea_section_s
 */
struct ref_section_s {
	int32_t id;
	uint32_t len;
	uint64_t base;
};
typedef struct ref_section_s ref_section_t;

/**
 * @fn ref_prec_init
 * @brief initialize mutable reference object (reference index precursor)
 */
ref_prec_t *ref_prec_init(
	ref_params_t const *params);

/**
 * @fn ref_prec_clean
 */
void ref_prec_clean(
	ref_prec_t *prec);

/**
 * @fn ref_append_segment
 * @brief append a sequence block to the context.
 */
int ref_append_segment(
	ref_prec_t *prec,
	char const *name,
	int32_t name_len,
	char const *seq,
	int64_t seq_len);

/**
 * @fn ref_append_link
 *
 * @brief append a edge on graph (not implemented yet)
 */
int ref_append_link(
	ref_prec_t *prec,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori);

/**
 * @fn ref_build_index
 * @brief build index
 */
ref_t *ref_build_index(
	ref_prec_t *prec);

/**
 * @fn ref_clean
 */
void ref_clean(
	ref_t *ref);

/**
 * @fn ref_load_index
 */
ref_t *ref_load_index(
	zf_t *fp);

/**
 * @fn ref_dump_index
 */
int ref_dump_index(
	ref_t const *ref,
	zf_t *fp);

/**
 * @fn ref_match
 */
struct ref_intv_s ref_match(
	ref_t *ref,
	void *seq);
struct ref_intv_s ref_match_2bit8packed(
	ref_t *ref,
	uint64_t seq);

/**
 * @fn ref_get_section
 */
struct ref_section_s const *ref_get_section(
	ref_t *ref,
	int32_t id);

/**
 * @fn ref_get_prefix
 */
char const *ref_get_prefix(
	ref_t const *ref);

/**
 * @fn ref_get_ptr
 */
uint8_t const *ref_get_ptr(
	ref_t const *ref);

/**
 * @fn ref_get_total_len
 */
int64_t ref_get_total_len(
	ref_t const *ref);

/**
 * @fn ref_get_ann
 *
 * @brief get annotation string at the given position.
 */
char const *ref_get_ann_name(
	ref_t const *ref,
	int64_t pos);
#if 0
/**
 * @fn ref_get_bound
 *
 * @brief get bounds of the sequence
 */
ref_interval_t ref_get_bound(
	ref_t const *ref,
	int64_t pos);
#endif
/**
 * @fn ref_is_amb
 */
int64_t ref_is_amb(
	ref_t const *ref,
	int64_t lb, int64_t ub);

/**
 * @fn ref_is_uniq
 */
int64_t ref_is_uniq(
	ref_t const *ref,
	int64_t lb, int64_t ub);

#endif /** #ifndef _REF_H_INCLUDED */
/**
 * end of ref.h
 */
