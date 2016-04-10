
/**
 * @file gref.h
 *
 * @brief a header of gref.c
 *
 * @author Hajime Suzuki
 * @date 2015/6/23
 * @license Apache v2.
 *
 * @detail
 * Ref is a hash-based sequence indexer and exact matcher. 
 */

#ifndef _GREF_H_INCLUDED
#define _GREF_H_INCLUDED

#include <stdint.h>

/**
 * @enum gref_error
 * @brief error flags
 */
enum gref_error {
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
 * @enum gref_format
 */
enum gref_format_flags {
	REF_ASCII				= 1,
	REF_2BIT				= 2,
	REF_2BIT8PACKED			= 3
};

/**
 * @type gref_prec_t
 */
typedef struct gref_prec_s gref_prec_t;

/**
 * @type gref_t
 * @brief an alias to struct gref_s
 */
typedef struct gref_s gref_t;

/**
 * @struct gref_params_s
 */
struct gref_params_s {
	int32_t seed_length;
	uint32_t reserved;
};
typedef struct gref_params_s gref_params_t;
#define REF_PARAMS(...)			( &((struct gref_params_s const) { __VA_ARGS__ }) )

/**
 * @struct gref_section_s
 * @brief has equivalent fields to struct sea_section_s
 */
struct gref_section_s {
	uint32_t id;
	uint32_t len;
	uint64_t base;
};
typedef struct gref_section_s gref_section_t;

/**
 * @struct gref_str_s
 */
struct gref_str_s {
	char const *str;
	int32_t len;
};
typedef struct gref_str_s gref_str_t;

/**
 * @struct gref_gid_pos_s
 */
struct gref_gid_pos_s {
	uint32_t gid;
	uint32_t pos;
};
typedef struct gref_gid_pos_s gref_gid_pos_t;

/* encode and decode id */
#define gref_rev(_gid)			( 0x01 ^ (_gid) )
#define gref_gid(_id, _d)		( ((_id)<<1) | (0x01 & (_d)) )
#define gref_id(_gid)			( (_gid)>>1 )
#define gref_dir(_gid)			( (_gid) & 0x01 )

/**
 * @struct gref_match_res_s
 */
struct gref_match_res_s {
	struct gref_gid_pos_s *ptr;
	int64_t len;
};
typedef struct gref_match_res_s gref_match_res_t;

/**
 * @fn gref_prec_init
 * @brief initialize mutable reference object (reference index precursor)
 */
gref_prec_t *gref_prec_init(
	gref_params_t const *params);

/**
 * @fn gref_prec_clean
 */
void gref_prec_clean(
	gref_prec_t *prec);

/**
 * @fn gref_append_segment
 * @brief append a sequence block to the context.
 */
int gref_append_segment(
	gref_prec_t *prec,
	char const *name,
	int32_t name_len,
	char const *seq,
	int64_t seq_len);

/**
 * @fn gref_append_link
 *
 * @brief append a edge on graph (not implemented yet)
 */
int gref_append_link(
	gref_prec_t *prec,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori);

/**
 * @fn gref_build_index
 * @brief build index
 */
gref_t *gref_build_index(
	gref_prec_t *prec);

/**
 * @fn gref_clean
 */
void gref_clean(
	gref_t *ref);

/**
 * @fn gref_load_index
 */
gref_t *gref_load_index(
	zf_t *fp);

/**
 * @fn gref_dump_index
 */
int gref_dump_index(
	gref_t const *ref,
	zf_t *fp);

/**
 * @fn gref_match
 */
struct gref_match_res_s gref_match(
	gref_t const *ref,
	char const *seq);
struct gref_match_res_s gref_match_2bitpacked(
	gref_t const *ref,
	uint64_t seq);

/**
 * @fn gref_get_section
 */
struct gref_section_s const *gref_get_section(
	gref_t const *ref,
	uint32_t id);

/**
 * @fn gref_get_name
 */
struct gref_str_s gref_get_name(
	gref_t const *ref,
	uint32_t id);

/**
 * @fn gref_get_ptr
 */
uint8_t const *gref_get_ptr(
	gref_t const *ref);

/**
 * @fn gref_get_total_len
 */
int64_t gref_get_total_len(
	gref_t const *ref);

/**
 * @fn gref_is_amb
 */
int64_t gref_is_amb(
	gref_t const *ref,
	int64_t lb, int64_t ub);

#endif /** #ifndef _GREF_H_INCLUDED */
/**
 * end of gref.h
 */
