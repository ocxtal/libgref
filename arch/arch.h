
/**
 * @file arch.h
 *
 * @brief architecture-dependent settings and utils
 */
#ifndef _ARCH_H_INCLUDED
#define _ARCH_H_INCLUDED

#include <stdint.h>

#ifdef __x86_64__

/* map reverse-complement sequence out of the canonical-formed address */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifdef AARCH64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifdef PPC64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifndef GREF_SEQ_LIM
#  error "No architecuture detected. Check CFLAGS."
#endif


#endif /* _ARCH_H_INCLUDED */
/**
 * end of arch.h
 */
