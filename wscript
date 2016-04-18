#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('psort')
	opt.recurse('hmap')
	opt.recurse('zf')
	opt.load('compiler_c')

def configure(conf):
	conf.recurse('psort')
	conf.recurse('hmap')
	conf.recurse('zf')

	conf.load('ar')
	conf.load('compiler_c')
	conf.env.append_value('LIB_GREF',
		conf.env.LIB_PSORT + conf.env.LIB_HMAP + conf.env.LIB_ZF)
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.recurse('psort')
	bld.recurse('hmap')
	bld.recurse('zf')

	bld.stlib(
		source = ['gref.c'],
		target = 'gref',
		lib = bld.env.LIB_GREF,
		use = ['psort', 'hmap', 'zf'])

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		linkflags = ['-all_load'],
		use = ['gref'],
		lib = bld.env.LIB_GREF,
		defines = ['TEST'])
