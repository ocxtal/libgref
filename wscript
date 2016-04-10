#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.load('compiler_c')

def configure(conf):
	conf.recurse('psort')
	conf.recurse('hmap')
	conf.recurse('zf')

	conf.load('ar')
	conf.load('compiler_c')
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.recurse('psort')
	bld.recurse('hmap')
	bld.recurse('zf')

	bld.stlib(
		source = ['ref.c'],
		target = 'ref',
		lib = bld.env.LIB_PSORT + bld.env.LIB_HMAP + bld.env.LIB_ZF,
		use = ['psort', 'hmap', 'zf'])

	bld.program(
		source = ['ref.c'],
		target = 'unittest',
		lib = bld.env.LIB_PSORT + bld.env.LIB_HMAP + bld.env.LIB_ZF,
		use = ['psort', 'hmap', 'zf'],
		defines = ['TEST'])
