#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <semaphore.h>

#define PROGRAM_NAME "ca"

#define DETERMINE_THREADCOUNT_DYNAMICALLY
#define DEFAULT_NUMBER_OF_THREADS 18

#define SPACE_PRINTING

#define fplog(     F, STR, ...) fprintf(F, STR "\n", ##__VA_ARGS__)
#define fperr(     F, STR, ...) fplog(F, "[%s:%s:%s:%d]: " STR, PROGRAM_NAME, __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__)
#define  plog(STR, ...) fplog(stderr, STR, ##__VA_ARGS__)
#define  perr(STR, ...) fperr(stderr, STR, ##__VA_ARGS__)

#define SWORD(IND) ((IND) >> 6)
#define SRANK(IND) ((IND) & 0x3Ful)
#define BIT(IND) (0x01ul << (IND))

typedef struct {
	int tid;
} thread_arg_t;

void print_usage_and_exit(void);
int parse_args(int argc, char *argv[]);
int allocate_states(void);
int init_state(int argc, char *argv[]);
void free_states(void);
//void print_uint64_t(uint64_t w);
//void print_state(uint64_t *state);
void print_results(void);
void fill_table(void);
void wrap_state(uint64_t *state);
uint64_t advance_gen_with_bounds(const uint64_t *__restrict__ is, uint64_t *__restrict__ os, const uint64_t b, const uint64_t e);
uint64_t advance_gen_last_word(const uint64_t *__restrict__ is, uint64_t *__restrict__ os);
int init_synchronization_objects(void);
int destroy_synchronization_objects(void);
void determine_num_threads(void);
int init_threads(void);
int join_threads(void);
void *thread_main(thread_arg_t *arg);

char rule = 0;
uint64_t width;
uint64_t gens;
uint64_t gen = 1;
uint64_t cnt;
uint64_t *states[2];
uint64_t statelen;
int curstate = 0;
uint8_t table[1024];
uint64_t *cnts;
pthread_barrier_t barrier;
sem_t sem;
//pthread_mutex_t mutex;
pthread_t *threads;
int nthreads;
int nthreads_left;
uint64_t minpop = 0, mingen = 0, maxpop, maxgen = 0;

int main(int argc, char *argv[]){
	if(parse_args(argc, argv)) print_usage_and_exit();
	determine_num_threads();
	if(allocate_states()) return 1;
	if(init_state(argc, argv)) print_usage_and_exit();
	fill_table();
	if(init_synchronization_objects()) return 1;
	if(init_threads()) return 1;

	uint64_t b, e;

	b = (statelen * nthreads) / (nthreads + 1);
	e = statelen - 1;

	// print_state(states[0]);

	int i;

	for(curstate = 1; gen <= gens; gen++){
		curstate = !curstate;
		nthreads_left = nthreads - 1;
		wrap_state(states[curstate]);
		pthread_barrier_wait(&barrier);

		cnt = advance_gen_with_bounds(states[curstate], states[!curstate], b, e) +
		      advance_gen_last_word  (states[curstate], states[!curstate]);

		sem_wait(&sem);
		// print_state(states[!curstate]);
		for(i = 0; i < nthreads; i++){
			cnt += cnts[i];
		}
		cnt -= states[!curstate][0] & 0x01ul;

		if(cnt > maxpop){
			maxpop = cnt;
			maxgen = gen;
		} else
		if(cnt < minpop){
			minpop = cnt;
			mingen = gen;
		}
	}
	pthread_barrier_wait(&barrier);

	free_states();
	if(join_threads()) return 1;
	if(destroy_synchronization_objects()) return 1;
	print_results();

	return 0;
}

int init_synchronization_objects(void){
	int r;
	r = pthread_barrier_init(&barrier, NULL, nthreads + 1);
	if(r == EAGAIN) perr("EAGAIN");
	if(r == EINVAL) perr("EINVAL");
	if(r == ENOMEM) perr("ENOMEM");
	if(r) return 1;

	r = sem_init(&sem, 0, 0);
	if(r == EINVAL) perr("EINVAL");
	if(r == ENOSYS) perr("ENOSYS");
	if(r) return 1;

	// r = pthread_mutex_init(&mutex, NULL);
	// if(r == EAGAIN) perr("EAGAIN");
	// if(r == ENOMEM) perr("ENOMEM");
	// if(r == EPERM ) perr("EPERM" );
	// if(r == EBUSY ) perr("EBUSY" );
	// if(r == EINVAL) perr("EINVAL");
	// if(r) return 1;

	return 0;
}

int destroy_synchronization_objects(void){
	int r;
	r = pthread_barrier_destroy(&barrier);
	if(r == EBUSY ) perr("EBUSY" );
	if(r == EINVAL) perr("EINVAL");
	if(r) return 1;

	r = sem_destroy(&sem);
	if(r == EINVAL) perr("EINVAL");
	if(r) return 1;

	// r = pthread_mutex_destroy(&mutex);
	// if(r == EBUSY ) perr("EBUSY" );
	// if(r == EINVAL) perr("EINVAL");
	// if(r) return 1;

	return 0;
}

void determine_num_threads(void){
	#if defined(DETERMINE_THREADCOUNT_DYNAMICALLY) && defined(_SC_NPROCESSORS_ONLN)
		// So, basically below a significantly large statelen, it's basically not worth multithreading.
		unsigned int ncores;
		ncores = sysconf(_SC_NPROCESSORS_ONLN);
		if(statelen > ncores * 1024){
			nthreads = (ncores * 3) / 2; // Let's really saturate those cores!
		} else
		if(statelen > ncores * 16){
			if(ncores > 3){
				nthreads = ncores - 1;
			} else {
				nthreads = 2;
			}
		} else {
			nthreads = 2;
		}
	#else
		nthreads = DEFAULT_NUMBER_OF_THREADS;
	#endif
	nthreads--; // nthreads is _really_ the count of _other_ non-main threads. The main thread does work too, though.
}

int init_threads(void){
	threads = malloc(nthreads * sizeof(pthread_t));
	if(threads == NULL){
		perr("malloc fail.");
		return 1;
	}

	cnts = malloc(nthreads * sizeof(uint64_t));
	if(cnts == NULL){
		perr("malloc fail.");
		return 1;
	}

	thread_arg_t *thread_args;
	thread_args = malloc(nthreads * sizeof(thread_arg_t));
	if(thread_args == NULL){
		perr("malloc fail.");
		return 1;
	}

	int i, r;
	for(i = 0; i < nthreads; i++){
		thread_args[i] = (thread_arg_t) {
			.tid      = i
		};
		r = pthread_create(&threads[i], NULL, (void *(*)(void *)) thread_main, (void *) &thread_args[i]);
		if(r == EAGAIN) perr("EAGAIN");
		if(r == EINVAL) perr("EINVAL");
		if(r == EPERM ) perr("EPERM" );
		if(r) return 1;
	}
	return 0;
}

int join_threads(void){
	int i, r;
	for(i = 0; i < nthreads; i++){
		r = pthread_join(threads[i], NULL);
		if(r == EDEADLK) perr("EDEADLK");
		if(r == EINVAL ) perr("EINVAL" );
		if(r == ESRCH  ) perr("ESRCH"  );
	}
	free(threads);
	threads = NULL;
	return 0;
}

void *thread_main(thread_arg_t *arg){
	int tleft;
	uint64_t b, e;

	b = (statelen * arg->tid) / (nthreads + 1);
	e = (statelen * (arg->tid + 1)) / (nthreads + 1);

	pthread_barrier_wait(&barrier);
	while(gen <= gens){

		cnts[arg->tid] = advance_gen_with_bounds(states[curstate], states[!curstate], b, e);

		do {
			tleft = nthreads_left;
		} while(!__sync_bool_compare_and_swap(&nthreads_left, tleft, tleft - 1));
		if(!tleft){
			sem_post(&sem);
		}
		pthread_barrier_wait(&barrier);
	}

	return NULL;
}

uint64_t advance_gen_with_bounds(const uint64_t *__restrict__ is, uint64_t *__restrict__ os, const uint64_t b, const uint64_t e){
	uint64_t prev, cur, next, cnt, o;
	
	if(b){
		prev = is[b - 1];
	} else {
		prev = 0;
	}
	cur = is[b];

	cnt = 0;
	
	uint64_t i;
	for(i = b; i < e; i++){
		next = is[i + 1];
		
		o = (((uint64_t) table[((cur <<  1) | (prev >> 63)) & 0x3FFul]) <<  0) +
			(((uint64_t) table[(cur >>  7) & 0x3FFul]) <<  8) +
			(((uint64_t) table[(cur >> 15) & 0x3FFul]) << 16) +
			(((uint64_t) table[(cur >> 23) & 0x3FFul]) << 24) +
			(((uint64_t) table[(cur >> 31) & 0x3FFul]) << 32) +
			(((uint64_t) table[(cur >> 39) & 0x3FFul]) << 40) +
			(((uint64_t) table[(cur >> 47) & 0x3FFul]) << 48) +
			(((uint64_t) table[((cur >> 55) & 0x1FFul) | ((next & 0x01ul) << 9)]) << 56);
		
		os[i] = o;
		
		cnt += __builtin_popcountl(o);

		prev = cur;
		cur = next;
	}

	return cnt;
}

uint64_t advance_gen_last_word(const uint64_t *__restrict__ is, uint64_t *__restrict__ os){
	uint64_t prev, cur, w, o;
	w = SWORD(width - 1);
	prev = is[w - 1];
	cur = is[w];
	o = ((((uint64_t) table[((cur <<  1) | (prev >> 63)) & 0x3FFul]) <<  0) +
		 (((uint64_t) table[(cur >>  7) & 0x3FFul]) <<  8) +
		 (((uint64_t) table[(cur >> 15) & 0x3FFul]) << 16) +
		 (((uint64_t) table[(cur >> 23) & 0x3FFul]) << 24) +
		 (((uint64_t) table[(cur >> 31) & 0x3FFul]) << 32) +
		 (((uint64_t) table[(cur >> 39) & 0x3FFul]) << 40) +
		 (((uint64_t) table[(cur >> 47) & 0x3FFul]) << 48) +
		 (((uint64_t) table[(cur >> 55) & 0x1FFul]) << 56)) & (BIT(SRANK(width - 1)) - 1);
	os[w] = o;
	return __builtin_popcountl(o);
}

void print_usage_and_exit(void){
	printf(
		"Usage: " PROGRAM_NAME " RULE_SPEC GRID_WIDTH NUM_GENERATIONS [INDICES...]\n\n"
		"\tRULE_SPEC        The rule, in ASCII-fied binary, e.g., 01110110.\n"
		"\tGRID_WIDTH       The circumference of the cylinder on which the simulation will run.\n"
		"\tNUM_GENERATIONS  The number of generations to run the simulation for.\n"
		"\t[INDICES...]     Which indices in the grid to initialize to 1.\n\n"
	);
	exit(EXIT_FAILURE);
}

int parse_args(int argc, char *argv[]){
	if(argc < 4){
		return 1;
	}

	int i;
	for(i = 0; i < 8; i++){
		switch(argv[1][i]){
			case '0':
				break;
			case '1':
				rule += BIT(i);
				break;
			default:
				return 1;
		}
	}

	errno = 0;
	width = strtoull(argv[2], NULL, 10) + 2;
	if(errno) return 1;
	gens  = strtoull(argv[3], NULL, 10);
	if(errno) return 1;

	statelen = SWORD(width - 1) + 1;

	return 0;
}

int allocate_states(void){
	uint64_t *s;

	s
	 = malloc(statelen * sizeof(uint64_t));
	if(s == NULL){
		perr("Malloc failed.");
		return 1;
	}
	memset(s, 0x00, statelen);
	states[0] = s;

	s = malloc(statelen * sizeof(uint64_t));
	if(s == NULL){
		perr("Malloc failed.");
		return 1;
	}
	memset(s, 0x00, statelen);
	states[1] = s;

	return 0;
}

int init_state(int argc, char *argv[]){
	uint64_t *s;
	uint64_t index;

	s = states[curstate];

	errno = 0;

	int i;
	for(i = 4; i < argc; i++){
		index = strtoull(argv[i], NULL, 10) + 1;
		if(errno) return 1;
		if(index >= width - 1) return 1;
		s[SWORD(index)] |= BIT(SRANK(index));
		minpop += 1;
	}

	maxpop = minpop;

	return 0;
}

void free_states(void){
	free(states[0]);
	free(states[1]);

	states[0] = NULL;
	states[1] = NULL;
}

/*
void print_uint64_t(uint64_t w){
	int j;
	char y = '1', n = '0';
	for(j = 0; j < 32; j += 8){
		printf(    "%c%c%c%c" "%c%c%c%c"
			#ifdef SPACE_PRINTING
			" "
			#endif
			,
			((w >> (j + 0)) & 1) ? y : n, ((w >> (j + 1)) & 1) ? y : n,
			((w >> (j + 2)) & 1) ? y : n, ((w >> (j + 3)) & 1) ? y : n,
			((w >> (j + 4)) & 1) ? y : n, ((w >> (j + 5)) & 1) ? y : n,
			((w >> (j + 6)) & 1) ? y : n, ((w >> (j + 7)) & 1) ? y : n
		);
	}
	for(j = 32; j < 64; j += 8){
		printf(
			#ifdef SPACE_PRINTING
			" "
			#endif
			"%c%c%c%c" "%c%c%c%c"    ,
			((w >> (j + 0)) & 1) ? y : n, ((w >> (j + 1)) & 1) ? y : n,
			((w >> (j + 2)) & 1) ? y : n, ((w >> (j + 3)) & 1) ? y : n,
			((w >> (j + 4)) & 1) ? y : n, ((w >> (j + 5)) & 1) ? y : n,
			((w >> (j + 6)) & 1) ? y : n, ((w >> (j + 7)) & 1) ? y : n
		);
	}
	return;
}
*/

/*
void print_state(uint64_t *state){
	uint64_t i, s;
	for(i = 0; i < statelen; i++){
		s = state[i];
		print_uint64_t(s);
		#ifdef SPACE_PRINTING
		printf("    ");
		#endif
	}
	printf("\n");
}
*/

void print_results(void){
	printf("Minimum popcount: %ld at step %ld\n", minpop, mingen);
	printf("Maximum popcount: %ld at step %ld\n", maxpop, maxgen);
	printf("Final popcount: %ld at step %ld\n", cnt, gen - 1);
}

void fill_table(void){
	int i, j;
	uint8_t o;
	for(i = 0; i < 1024; i++){
		o = 0;
		for(j = 0; j < 8; j++){
			o += ((rule & (1 << ((i >> j) & 7))) != 0) << j;
		}
		table[i] = o;
	}
}

void wrap_state(uint64_t *state){
	// if(state[0] & BIT(1)){
	// 	state[SWORD(width - 1)] = (state[SWORD(width - 1)] & (BIT(SRANK(width - 1)) - 1)) | BIT(SRANK(width - 1));
	// } else {
	// 	state[SWORD(width - 1)] = (state[SWORD(width - 1)] & (BIT(SRANK(width - 1)) - 1));
	// }
	// if(state[SWORD(width - 2)] & BIT(SRANK(width - 2))){
	// 	state[0] |=  BIT(0);
	// } else {
	// 	state[0] &= ~BIT(0);
	// }
	if(state[0] & BIT(1)){
		state[SWORD(width - 1)] |= BIT(SRANK(width - 1));
	} else {
		state[SWORD(width - 1)] &= BIT(SRANK(width - 1)) - 1;
	}
	if(state[SWORD(width - 2)] & BIT(SRANK(width - 2))){
		state[0] |=  BIT(0);
	} else {
		state[0] &= ~BIT(0);
	}
}
