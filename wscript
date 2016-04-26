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

	conf.env.append_value('CFLAGS', '-g')
	conf.env.append_value('CFLAGS', '-Wall')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')

	conf.env.append_value('LIB_GREF',
		conf.env.LIB_PSORT + conf.env.LIB_HMAP + conf.env.LIB_ZF)
	conf.env.append_value('OBJ_GREF',
		['gref.o'] + conf.env.OBJ_PSORT + conf.env.OBJ_HMAP + conf.env.OBJ_ZF)


def build(bld):
	bld.recurse('psort')
	bld.recurse('hmap')
	bld.recurse('zf')

	bld.objects(source = 'gref.c', target = 'gref.o')

	print(bld.env.LIB_GREF)
	print(bld.env.OBJ_GREF)

	bld.stlib(
		source = ['unittest.c'],
		target = 'gref',
		use = bld.env.OBJ_GREF,
		lib = bld.env.LIB_GREF)

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		use = bld.env.OBJ_GREF,
		lib = bld.env.LIB_GREF,
		defines = ['TEST'])
