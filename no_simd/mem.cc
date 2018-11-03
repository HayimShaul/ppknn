#if 0
#include <malloc.h>
#include <iostream>
//#include "mem.h"

int crash_at = -1;

static int mem_cnt = 0;
static int del_cnt = 0;

void *operator new[](size_t s) {
	if (mem_cnt == crash_at)
		*((int *)(0)) = 5;
	void *p = malloc(s);
	std::cerr << std::endl << p << " " << mem_cnt++ << " " << s << " allocated MEM\n";
	return p;
}

void operator delete[](void *p) {
	std::cerr << std::endl << p << " " << mem_cnt++ << " " << del_cnt++ << " freed MEM\n";
	free(p);
}

void *operator new(size_t s) {
	if (mem_cnt == crash_at)
		*((int *)(0)) = 5;
	void *p = malloc(s);
	std::cerr << std::endl << p << " " << mem_cnt++ << " " << s << " allocated MEM\n";
	return p;
}

void operator delete(void *p) {
	std::cerr << std::endl << p << " " << mem_cnt++ << " " << del_cnt++ << " freed MEM\n";
	free(p);
}
#endif
