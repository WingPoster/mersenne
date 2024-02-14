// compile : gcc-13 -fopenmp mersenne.c -o mersenne
// run : ./mersenne option.config
// console : ./mersenne console.config
// kill : ./mersenne kill.config

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <signal.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <arpa/inet.h> // socket
#include <sys/types.h> // socket
#include <sys/socket.h> // socket
#include <stdint.h> // block encrypt
#include <stddef.h> // block encrypt
#include <inttypes.h> // asymmetric encrypt
#include <sys/stat.h>
#include <complex.h>
#include <omp.h>

#ifdef _WIN32
#include <io.h>
#define F_OK 0
#define access _access
#else
#include <unistd.h>
#endif

#define MAX_SIZE            1024
#define MAX_PATH            1024
#define MAX_NOISE           1024
#define MAX_PROCESS         1024
#define MAX_BUFF_SIZE       8096
#define MAX_NUMBER         65535
#define MAX_CHAR             256
#define MAX_BIT                8
#define MAX_DIFF             256 // min 512 difference
#define MAX_PRIMES          1000
#define MAX_NOISES           100
#define MAX_DATE              26
#define MAX_INTEGER     10000000

#define MIN_CANDIDATE         37
#define MAX_CANDIDATE 1000000000
#define MIN_PRIME             37
#define MAX_PRIME     1000000000
#define MAX_CHUNK           1024 

#define MERSENNE_NUMBER        0
#define MERSENNE_CANDIDATE     1
#define KNOWN_MERSENNE        51
#define PI 3.14159265358979323846

#define EXPERIMENTAL
#ifdef EXPERIMENTAL
// AI Framework
#define INPUT_SIZE 2
#define HIDDEN_SIZE 100
#define OUTPUT_SIZE 1
#define LEARNING_RATE 0.1
#define NUM_EPOCHS 10000
#define NUM_THREADS 4

// Mersenne Twister constants
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908b0df
#define MT_UPPER_MASK 0x80000000
#define MT_LOWER_MASK 0x7fffffff

// BigInt
#define MAX_DIGITS 100000
#define SIZE_INPUT_CHUNK MAX_DIGITS/8
#define SIZE_OUTPUT_CHUNK MAX_DIGITS/12
#define BITS_PER_DIGIT 32 // Assuming each digit is a 32-bit unsigned integer
#define PRECISION 50
#define BASE 10

#endif // #ifdef EXPERIMENTAL

// Encryption
#define blockSize             16
#define KeySize               32
#define KeyExpSize           240
#define Nb                     4
#define Nk                     8
#define Nr                    14

#define E_VALUE                3 // 65535

#define get_sb(num) (sb[(num)])
#define get_sb_inv(num) (rsb[(num)])

#define multiply(x, y)                                \
      (  ((y & 1) * x) ^                              \
      ((y>>1 & 1) * xtime(x)) ^                       \
      ((y>>2 & 1) * xtime(xtime(x))) ^                \
      ((y>>3 & 1) * xtime(xtime(xtime(x)))) ^         \
      ((y>>4 & 1) * xtime(xtime(xtime(xtime(x))))))   \

// Main Process
#define RUN_ALL           0
#define RUN_TDM           1
#define RUN_LLT           2
#define RUN_SERVER        3
#define RUN_CONSOLE       4
#define RUN_TEST          5
#define RUN_BENCHMARK     6
#define RUN_KILL          7

#define FINDER_NORMAL     0
#define FINDER_NOMOD      1
#define FINDER_SKIP       2
#define FINDER_EULER      3

#define PROCESS_ID        3
#define PROCESS_ID_CLASS  0
#define PROCESS_ID_PID    1
#define PROCESS_ID_STATUS 2

#define PRIME_ID          3
#define PRIME_ID_USE      0
#define PRIME_ID_VALUE    1
#define PRIME_ID_PRIME    2

#define ALGO_KJ           0
#define ALGO_MULMOD       1

#define USE_NOFFT         0 
#define USE_FFT           1 

#define USE_NOFILE        0
#define USE_FILE          1 

#define USE_NOOMP         0 
#define USE_OMP           1 

#define USE_DIRECT        0 
#define USE_SOCKET        1 

#define LOG_NONE          0
#define LOG_ERROR         1
#define LOG_INFO          2
#define LOG_DEBUG         3
#define LOG_DUMP          4
#define LOG_FILE          5
#define LOG_NETWORK       6
#define LOG_VERBOSE       7

#define TDM_PROCESS_NAME    "TDM"
#define LLT_PROCESS_NAME    "LLT"

#define PROCESS_FILE_NAME   "process"
#define STATUS_FILE_NAME    "status"
#define LOG_FILE_NAME       "mersenne"
#define LOCK_FILE_NAME      "lock.pid"
#define DB_FILE_NAME        "db.dat"
#define CANDIDATE_FILE_NAME "candidate.dat"
#define LAST_FILE_NAME      "candidate.last"

#define STATUS_INIT         "init"
#define STATUS_READY        "ready"
#define STATUS_RUN          "run"
#define STATUS_STOP         "stop"
#define STATUS_KILL         "kill"
#define STATUS_EXIT         "exit"

#if defined(_WIN32)
    #define PLATFORM_NAME "windows" // Windows
#elif defined(_WIN64)
    #define PLATFORM_NAME "windows" // Windows
#elif defined(__CYGWIN__) && !defined(_WIN32)
    #define PLATFORM_NAME "windows" // Windows (Cygwin POSIX under Microsoft Window)
#elif defined(__ANDROID__)
    #define PLATFORM_NAME "android" // Android (implies Linux, so it must come first)
#elif defined(__linux__)
    #define PLATFORM_NAME "linux" // Debian, Ubuntu, Gentoo, Fedora, openSUSE, RedHat, Centos and other
#elif defined(__unix__) || !defined(__APPLE__) && defined(__MACH__)
    #include <sys/param.h>
    #if defined(BSD)
        #define PLATFORM_NAME "bsd" // FreeBSD, NetBSD, OpenBSD, DragonFly BSD
    #endif
#elif defined(__hpux)
    #define PLATFORM_NAME "hp-ux" // HP-UX
#elif defined(_AIX)
    #define PLATFORM_NAME "aix" // IBM AIX
#elif defined(__APPLE__) && defined(__MACH__) // Apple OSX and iOS (Darwin)
    #include <TargetConditionals.h>
    #if TARGET_IPHONE_SIMULATOR == 1
        #define PLATFORM_NAME "ios" // Apple iOS
    #elif TARGET_OS_IPHONE == 1
        #define PLATFORM_NAME "ios" // Apple iOS
    #elif TARGET_OS_MAC == 1
        #define PLATFORM_NAME "osx" // Apple OSX
    #endif
#elif defined(__sun) && defined(__SVR4)
    #define PLATFORM_NAME "solaris" // Oracle Solaris, Open Indiana
#else
    #define PLATFORM_NAME "unknown" // Unknown Platform
#endif

#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */

// Options
//#define USE_COMPRESSION
//#define USE_ENCRYPTION

#ifdef EXPERIMENTAL
#define BIGINT
#ifdef BIGINT
// Structure to represent big integers
typedef struct {
    int digits[MAX_DIGITS]; // Array to store digits
    int length;             // Number of digits
} BigInt;
#endif

typedef struct {
    uint64_t state[MT_N];
    int index;
} mt19937_state;

mt19937_state mt_state;

#endif

// Data Structure
typedef struct Node {
  int data;
  struct Node* next;
} node;

typedef struct Row {
  char class[MAX_SIZE];
  char pid[MAX_SIZE];
  char status[MAX_SIZE];
} row;

// Compression
typedef struct symbolInfo {
  unsigned char symbol;
  int frequency;
} symbolInfo;

typedef struct HFNode {
  symbolInfo data;
  struct HFNode* left;
  struct HFNode* right;
} HFNode;

typedef struct bitBuffer {
  unsigned int size;
  unsigned char* buffer;
}bitBuffer;

typedef struct HFCode {
  unsigned int size;
  unsigned char code[MAX_BIT];
} HFCode;

typedef int priorityType;

typedef struct PQNode {
  priorityType priority;
  void* data;
} PQNode;

typedef struct priorityQueue {
  PQNode* nodes;
  int capacity;
  int usedSize;
} priorityQueue;

// sequence A000043 in the OEIS
int knownMersenne[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091, 756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951, 30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933};

// Zeta Function
const long double LOWER_THRESHOLD = 1.0e-6;
const long double UPPER_BOUND = 1.0e+4;
const int MAXNUM = 100;

// Encryption
typedef uint8_t state_t[4][4];

static const uint8_t sb[256] = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

static const uint8_t rsb[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

static const uint8_t rcon[11] = {
  0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36 };

char * buf = NULL;
int len = 0;
uint8_t rkey[KeyExpSize];
uint8_t iv[blockSize];

// Database
struct Node* g_lstPrime = NULL;
struct Node* g_lstCandidate = NULL;
struct Node* g_pNextPrime = NULL;
int g_nChunk = MAX_CHUNK;
uint64_t g_nCandidate = MIN_CANDIDATE;
uint64_t p = 0;

// Memory Cache
uint64_t candidate[MAX_CANDIDATE/64 + 1];
uint64_t cache[MAX_CANDIDATE/64 + 1];

// Process
char g_pchClass[MAX_SIZE];
char g_pchPID[MAX_SIZE];
int g_nChildProcessID = 0;

// Configuration
void initMT19937();
void readCommand(int argc,char ** argv);
char g_achOptionPath[MAX_SIZE];

uint64_t g_nSystemMinCandidate = MIN_CANDIDATE;
uint64_t g_nSystemMaxCandidate = MAX_CANDIDATE;
uint64_t g_nMinCandidate = MIN_CANDIDATE;
uint64_t g_nMaxCandidate = MAX_CANDIDATE;
uint64_t g_nMinPrime = MIN_PRIME;
uint64_t g_nMaxPrime = MAX_PRIME;

int g_nRunMode = RUN_ALL;
int g_nAlgorithm = ALGO_KJ;
int g_nUseFFT = USE_FFT;
int g_nUseOMP = USE_OMP;
int g_nUseSocket = USE_SOCKET;
int g_nNoFile = USE_FILE;
int g_nFinder = FINDER_NORMAL;
int g_nTDMProcessCount = 1;
int g_nLLTProcessCount = 1;
int g_nPort = 8080;
int g_nLogFile = LOG_FILE;
int g_nLogNetwork = LOG_NETWORK;
int g_nLogVerbose = LOG_VERBOSE;
int g_nLogLevel = LOG_DEBUG;
int g_nDisplayStart = MIN_CANDIDATE;
int g_nDisplayEnd = MIN_CANDIDATE + MAX_SIZE;
int g_nDisplayCount = MAX_SIZE;
char g_achAddress[MAX_SIZE];

// Mersenne Library
uint64_t bits_max = 4294967296;
uint64_t bigA[MAX_INTEGER];
uint64_t bigB[MAX_INTEGER];
uint64_t bigC[MAX_INTEGER];
uint64_t bigS[MAX_INTEGER];
uint64_t bigL[MAX_INTEGER];
uint64_t bigBuff[MAX_INTEGER];
uint64_t bigTemp[MAX_INTEGER];
uint64_t bigTemp2[MAX_INTEGER];
double x[MAX_INTEGER];
double xi[MAX_INTEGER];
double rx[MAX_INTEGER];
double rxi[MAX_INTEGER];
double y[MAX_INTEGER];
double yi[MAX_INTEGER];
double ry[MAX_INTEGER];
double ryi[MAX_INTEGER];

void initInteger();
void add(uint64_t * src, uint64_t * dst);
void addop(uint64_t op, uint64_t * dst);
void set(uint64_t * dst);
void subtract(uint64_t * src, uint64_t * dst);
void subtractop(uint64_t op, uint64_t * dst);
void copy(uint64_t * src, uint64_t * dst);
void copyBits(uint64_t * src, uint64_t * dst, uint64_t nStart, uint64_t nEnd);
void mod(uint64_t * src, uint64_t * dst);
void mul(uint64_t * src1, uint64_t * src2, uint64_t * dst);
void mul_omp(uint64_t * src1, uint64_t * src2, uint64_t * dst);
void mul_fft(uint64_t * src1,uint64_t * src2, uint64_t * dst);
void mul_fft_internal(int * src1,int * src2, int * dst);
void digits_to_array(const int digits[], int num_digits, uint64_t result[]);
void array_to_digits(const uint64_t array[], int num_elements, int digits[]);
void square(uint64_t * src,uint64_t * dst);
void square_fft(uint64_t * src,uint64_t * dst);
void square_fft_internal(int * src,int * dst);
void mulop(uint64_t op, uint64_t * dst);
void leftshift(uint64_t * src);
void leftshiftop(uint64_t nBits, uint64_t * src);
void rightshift(uint64_t * src);
void rightshiftop(uint64_t nBits, uint64_t * src);
void PrintArray(uint64_t * src);
void PrintBits(uint64_t * src);
uint64_t getDigits(uint64_t * arr);
int getDigitsBased10(int * arr);
uint64_t getBits(uint64_t * arr);
uint64_t getBitsop(uint64_t value);
int compare(uint64_t * src, uint64_t * dst);
int compareop(uint64_t op, uint64_t * dst);
bool iszero();
void propa_carrier(uint64_t nIndex, uint64_t * dst);

void test_get_primes();
void primeDifferenceStatistics();

// Zeta Function ( need to convert big number version )
void test_riemann_zeta();
void riemann_zeta(double p,double q,double * r,double * ri,int maxNumber);
int euler_product(double p, double q, int maxNumber);
bool isPrime_kj(uint64_t n, int maxNumber);
double getPrimeCount_kj(int n);
double primeDistribution(int n);

void test_goldBach();
void goldBach();

#ifdef EXPERIMENTAL
void test_zeta();
complex double riemann_zeta_omp(complex double s);
void test_prime_gap_distribution();
void prime_gap_distribution(uint64_t start, uint64_t end);
void test_geometric_series();
complex double geometric_series(complex double first_term, complex double common_ratio, int num_terms);
void fft(complex double *x, int n);
void ifft(complex double *x, int n);
void multiply_bigint_fft(BigInt *a, BigInt *b, BigInt *result);
void square_bigint_fft(BigInt *a, BigInt *result);
void fft_omp(complex double *x, int n);
void ifft_omp(complex double *x, int n);
void multiply_bigint_fft_omp(BigInt *a, BigInt *b, BigInt *result);
void square_bigint_fft_omp(BigInt *a, BigInt *result);
void test_fft();
void mt19937_initialize(mt19937_state *state, uint32_t seed);
uint32_t mt19937_extract_number(mt19937_state *state);
uint64_t power(uint64_t base, uint64_t exp, uint64_t mod);
uint64_t power_ll(uint64_t a, uint64_t b, uint64_t mod);
bool miller_rabin(uint64_t n, int k);
void test_miller_rabin();
bool LehmannLucasPrimalityTest(uint64_t n, int num_trials);
void test_montgomery_multiply();
void montgomery_multiply_omp(BigInt *a, BigInt *b, BigInt *n, BigInt *n_inv, BigInt *r, BigInt *result);
void montgomery_multiply(BigInt *a, BigInt *b, BigInt *n, BigInt *n_inv, BigInt *r, BigInt *result);
void rsa_encrypt(BigInt *plaintext, BigInt *n, BigInt *e, BigInt *ciphertext);
void rsa_decrypt(BigInt *ciphertext, BigInt *n, BigInt *d, BigInt *plaintext);
bool test_lucaslehmer_omp(int exponent);
void test_lucaslehmer();
bool is_prime_aks(uint64_t n);
bool is_prime_aks_omp(uint64_t n);
bool is_prime_aks(uint64_t n);
void test_aks();
void test_count_primes();
uint64_t count_primes_eratosthenes(uint64_t n);
uint64_t count_primes_eratosthenes_omp(uint64_t n);
uint64_t count_primes_riemann(uint64_t n);
uint64_t count_primes_montgomery(uint64_t n);
void fermat_primes(int n);
void test_fermat_primes();
double sigmoid(double x);
double sigmoid_derivative(double x);
void forward_pass(double input[INPUT_SIZE], double hidden[HIDDEN_SIZE], double output[OUTPUT_SIZE], double weights_ih[INPUT_SIZE][HIDDEN_SIZE], double weights_ho[HIDDEN_SIZE][OUTPUT_SIZE]);
void backpropagation(double input[INPUT_SIZE], double hidden[HIDDEN_SIZE], double output[OUTPUT_SIZE], double target[OUTPUT_SIZE], double weights_ih[INPUT_SIZE][HIDDEN_SIZE], double weights_ho[HIDDEN_SIZE][OUTPUT_SIZE]);
void test_nn();
uint64_t factorial_mod_n(uint64_t n);
int is_prime_wilson(uint64_t n);
void test_wilson();
void test_monteCarloSimulation();

#ifdef BIGINT 
void test_bigint();
void free_bigint(BigInt *num);
void mod_multiply(const BigInt *a, const BigInt *b, const BigInt *n, BigInt *result);
int gcd_bigint(int a, int b);
uint64_t binomial(uint64_t n, uint64_t k);
uint64_t mod_exp(uint64_t base, uint64_t exp, uint64_t mod);

// Function to initialize a BigInt with a string representation of a number
void init_bigint_from_string(const char *num_str, BigInt *num);
void init_bigint(BigInt *num);
void init_bigint_value(BigInt *num, int value);
void assign_bigint(BigInt *dest, const BigInt *src);
void print_bigint(const BigInt *num);
void print_bigint_desc(char * desc, const BigInt *num);

// Function to compare two big integers (returns -1 if a < b, 0 if a == b, and 1 if a > b)
int compare_bigint(const BigInt *a,const BigInt *b);

// Function to add two big integers a and b and store the result in result
void add_bigint(const BigInt *a, const BigInt *b, BigInt *result);
void subtract_bigint(const BigInt *a,const BigInt *b, BigInt *result);

// Function to multiply two big integers
void multiply_bigint(const BigInt *a,const BigInt *b, BigInt *result);
void multiply_scalar_bigint(const BigInt *a,const int scalar, BigInt *result);
void divide_bigint(const BigInt *a, const BigInt *b, BigInt *quotient, BigInt *remainder);
void divide_scalar_bigint(const BigInt *a, const int scalar, BigInt *quotient);

// Function to perform modular arithmetic with big integers (a % b)
void mod_bigint(const BigInt *a, const BigInt *b, BigInt *result);
void square_bigint(const BigInt *a, BigInt *result);
void sqrt_bigint(const BigInt *a, BigInt *result);
void power_bigint(const BigInt *base, const int exponent, BigInt *result);

// Function to compute factorial of n using big integers
void factorial_bigint(const int n, BigInt *result);

// Function to perform left shift operation on a BigInt by a specified number of bits
void left_shift_bit_bigint(BigInt *num,const int shiftbit);

// Function to perform right shift operation on a BigInt by a specified number of bits
void right_shift_bit_bigint(BigInt *num,const int shiftbit);
void left_shift_bigint(BigInt *num, const int shiftdigits);
void right_shift_bigint(BigInt *num, const int shiftdigits);

// Function to compute factorial
double factorial(int n);

// Function to compute sine using Taylor series expansion
double sine(double x, int precision);

// Function to compute cosine using Taylor series expansion
double cosine(double x, int precision);

// Function to compute natural logarithm using Taylor series expansion
double logarithm(double x, int precision);

// Function to compute base 10 logarithm using Taylor series expansion
double logarithm_base10(double x, int precision);

// Function to compute pi using Machin's formula and Taylor series expansion
double pi(int precision);

void multiply_pi(int *arr, int multiplier);
void print_pi(int *arr);
void compute_pi();
#endif

#endif

// Main Process
void runAll();
void runTDM();
void runLLT();
void runServer();
void runConsole();
void runTest();
void runBenchmark();

// Mersenne Process
void launchLLT();
void launchTDM();

// LLT Process
void setMersenneNumber(uint64_t p);
bool PrimalityTesting(uint64_t p);
void LLTmethod(uint64_t * src, uint64_t *dst);
void LLTmulmod(uint64_t * src, uint64_t *dst);

// TDM Process
void findMersenne();
void findMersenneOMP();
void findChunkNormal();
void findChunkNomod();
void findChunkSkip();
void findChunkEuler();
void setLastTDM(long p);
long getLastTDM();

// Management Candidate
void initCandidate();
bool isPrime(uint64_t num);
bool isPrimeNormal(uint64_t num);
bool isPrimeFromDB(uint64_t num);
#ifdef EXPERIMENTAL
bool is_prime(uint64_t n);
#endif
bool isKnownMersenne(int number);
void setNumber(uint64_t nIndex, int bit);
bool setNumberCache(uint64_t nIndex, int bit);
bool isCandidate(uint64_t nIndex);
int getCandidateCount();
void PrintCandidateCount();
struct Node* getPrimeList();
struct Node* getCandidateList();
void freeList(struct Node* pList);
int getNextChunk(uint64_t chunk[][PRIME_ID]);

// Client Candidate Database
void readBits();
void writeBits();
void updateBits();

// Server Candidate Database
void readDB();
void writeDB();
void updateDB();

// Server Primes Database
void getPrime(uint64_t maxPrime, uint64_t * primes);
void findPrime(uint64_t maxPrime);
void writePrime(uint64_t maxPrime, uint8_t * primes);
void writeNoise(uint64_t maxPrime, uint32_t * noises);

// System Process
void updateProcess(char * pchClass, char * pchPID, char * pchStatus);
void waitTDMReady();
bool isTDMReady();
void killAll();

// Compression
#ifdef USE_COMPRESSION
void compress(char * pchFileName, uint64_t * data);
void uncompress(char * pchFileName, uint64_t * data);
void compressPrimes(uint64_t maxPrime);
void uncompressPrimes(uint64_t maxPrime, uint64_t * primes);

HFNode* createNode(symbolInfo data);
void destroyNode(HFNode* node);
void destroyTree(HFNode* tree);

HFNode* createNode(symbolInfo data);
void destroyNode(HFNode* node);
void destroyTree(HFNode* tree);

void addBit(bitBuffer* buffer, unsigned char bit);
void encode(HFNode** tree, const unsigned char *src,
                    bitBuffer* encoded, HFCode codeTable[MAX_CHAR]);
void decode(HFNode* tree, bitBuffer* encoded, unsigned char* decoded);
void buildPrefixTree(HFNode** tree, symbolInfo table[MAX_CHAR]);
void decode(HFNode* tree, bitBuffer* encoded, unsigned char* decoded);
void buildPrefixTree(HFNode** tree, symbolInfo table[MAX_CHAR]);
void buildCodeTable(HFNode* tree, HFCode table[MAX_CHAR],
                            unsigned char code[MAX_BIT], int size);
void printBinary(bitBuffer* buffer);

priorityQueue* Create(int size);
void destroy(priorityQueue* pq);
void enqueue(priorityQueue* pq, PQNode newNode);
void dequeue(priorityQueue* pq, PQNode* root);
int  getParent(int pos);
int  getLeftChild(int pos);
int  getRightChild(int pos);
void swapNodes(priorityQueue* pq, int pos1, int pos2);
int  isEmpty(priorityQueue* pq);
void printNode (PQNode* node);

void test_compress();
#endif

// Encryption
#ifdef USE_ENCRYPTION
static void key_exp(uint8_t* rkey, const uint8_t* key);
static void add_rkey(uint8_t round, state_t* state, const uint8_t* rkey);
static void subbytes(state_t* state);
static void shiftrows(state_t* state);
static uint8_t xtime(uint8_t x);
static void mixcolumns(state_t* state);
static void mixcolumns_inv(state_t* state);
static void subbytes_inv(state_t* state);
static void shiftrows_inv(state_t* state);
static void cipher(state_t* state, const uint8_t* rkey);
static void cipher_inv(state_t* state, const uint8_t* rkey);
static void xoriv(const uint8_t* iv);

void encryptFile(char * pchFileName, uint64_t * data);
void decryptFile(char * pchFileName, uint64_t * data);

// asymmetric encrypt
uint64_t findD(uint64_t e, uint64_t phi);
uint64_t gcd(uint64_t num1, uint64_t num2);
uint64_t getprime();
void setprimes(uint64_t e, uint64_t *p, uint64_t *q, uint64_t *n, uint64_t *phi);
uint64_t modpow(uint64_t base, uint64_t power, uint64_t mod);
uint64_t inverse(int a, int mod);
int test_sign ();

#endif

// Locking Protocol
void lock_enter();
void lock_leave();

bool isSingleProcess();
void getStatus(char * pchStatus);
void setStatus(char * pchStatus);

// System Time
time_t tStart = 0;

// Service Process
void server();
void agent(int sock);
void error_handling(char *message);

// Function Test Process
void test();

// Function Benchmark Process
void benchmark();

// System Console Process
void console();

// Configuation
int configInt(char *pKey);
void readConfig(char *pKey, char *pValue);
void initConfig();

// Logging
void logging(int level, const char *format, ...);
void writeLog(char * buf);
void sendLog(char *buf);

// Testing
void test_mersenne();

// Benchmark
void benchmark_multiply();

// Process Entry Point
int main(int argc, char ** argv) {
  initMT19937();

  readCommand(argc, argv);
  logging(LOG_INFO, "Platform %s\n", PLATFORM_NAME);

  if(g_nRunMode == RUN_ALL) {
    runAll();
  } else if(g_nRunMode == RUN_TDM) {
    runTDM();
  } else if(g_nRunMode == RUN_LLT) {
    runLLT();
  } else if(g_nRunMode == RUN_SERVER) {
    runServer();
  } else if(g_nRunMode == RUN_CONSOLE) {
    runConsole();
  } else if(g_nRunMode == RUN_TEST) {
    runTest();
  } else if(g_nRunMode == RUN_BENCHMARK) {
    runBenchmark();
  }else {
    logging(LOG_INFO, "Unknown Run Mode [%d].\n", g_nRunMode);
  }
}

void runAll() {
  logging(LOG_DEBUG, "This is launcher process. pid=[%ld]\n", (long)getpid()); 
  logging(LOG_INFO, "All process START!\n");
  launchTDM();
  launchLLT();
}

void runTDM(){
  logging(LOG_INFO, "TDM Process START!\n");

  char achPID[MAX_SIZE];
  memset(achPID, 0x00, MAX_SIZE);
  sprintf(achPID, "%ld", (long)getpid());
  updateProcess(TDM_PROCESS_NAME, achPID, STATUS_INIT);
  initCandidate();
  updateProcess(TDM_PROCESS_NAME, achPID, STATUS_READY);
  waitTDMReady();

  logging(LOG_INFO, "runTDM : TDM Ready!\n");
  g_lstPrime = getPrimeList();
  g_lstCandidate = getCandidateList();
  g_pNextPrime = g_lstPrime;
  if(g_nUseOMP == USE_OMP) {
    findMersenneOMP();
  } else {
    findMersenne();
  }
  freeList(g_lstPrime);
  freeList(g_lstCandidate);
}

void runLLT(){
  logging(LOG_INFO, "LLT Process START!\n");

  char achPID[MAX_SIZE];
  memset(achPID, 0x00, MAX_SIZE);
  sprintf(achPID, "%ld", (long)getpid());
  updateProcess(LLT_PROCESS_NAME, achPID, STATUS_INIT);
  waitTDMReady();
  updateProcess(LLT_PROCESS_NAME, achPID, STATUS_READY);
  logging(LOG_INFO, "runLLT : TDM Ready!\n");
  readBits();

  if(g_nRunMode == RUN_ALL || g_nRunMode == RUN_LLT) {
    while(1) {
      updateBits();
      g_nCandidate = g_nMinCandidate;
      for(int i = g_nMinCandidate; i < g_nMaxCandidate; i+=2) {
        if(isCandidate(i) && !isKnownMersenne(i)) {
          g_nCandidate = i;
          logging(LOG_DEBUG, "PrimalityTesting select=[%d]\n", g_nCandidate);
          break;
        }
      }
      if(PrimalityTesting(g_nCandidate)) {
        logging(LOG_INFO, "2^%d - 1 is prime!!\n", g_nCandidate);
      } else {
        setNumber(g_nCandidate, MERSENNE_NUMBER);
      }
    }
  } else {
    logging(LOG_INFO, "Unknown Run Mode [%d].\n", g_nRunMode);
  }
  logging(LOG_INFO, "PrimalityTesting Complete!\n");
}

void runServer() {
  logging(LOG_INFO, "Server Mode START!\n");
  server();
}

void runConsole() {
  logging(LOG_INFO, "Console Mode START!\n");
  console();
}

void runTest() {
  logging(LOG_INFO, "Test Mode START!\n");
  test();
}

void runBenchmark() {
  logging(LOG_INFO, "Benchmark Mode START!\n");
  benchmark();
}

void launchLLT() {
  pid_t pid;

  g_nSystemMinCandidate = g_nMinCandidate;
  g_nSystemMaxCandidate = g_nMaxCandidate;

  for(int i = 1; i <= g_nLLTProcessCount; i++) {
    pid = fork();
    int nMinCandidate = 0;
    int nMaxCandidate = 0;
    nMinCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nLLTProcessCount)*(i-1);
    if(i != g_nLLTProcessCount) {
      nMaxCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nLLTProcessCount)*i - 1;
    } else {
      nMaxCandidate = g_nMaxCandidate;
    }

    if(nMinCandidate < MIN_CANDIDATE) {
      nMinCandidate = MIN_CANDIDATE;
    }
    if(nMinCandidate%2==0) nMinCandidate++;

    switch(pid) {
      case -1:
        logging(LOG_INFO, "fork failed, reason=failed create child process.\n");
        return;
      case 0:
        logging(LOG_DEBUG, "This is child LLT process. pid=[%ld]\n", (long)getpid());
        g_nChildProcessID = pid;
        g_nMinCandidate = nMinCandidate;
        g_nMaxCandidate = nMaxCandidate;
        runLLT();
        logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
        break;
      default:
        break;
    }
  }
}

void launchTDM() {
  pid_t pid;

  g_nSystemMinCandidate = g_nMinCandidate;
  g_nSystemMaxCandidate = g_nMaxCandidate;

  for(int i = 1; i <= g_nTDMProcessCount; i++) {
    pid = fork();

    int nMinPrime = 0;
    int nMaxPrime = 0;
    nMinPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTDMProcessCount)*(i-1);
    if(i != g_nTDMProcessCount) {
      nMaxPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTDMProcessCount)*i - 1;
    } else {
      nMaxPrime = g_nMaxPrime;
    }
    int nMinCandidate = 0;
    int nMaxCandidate = 0;
    nMinCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nTDMProcessCount)*(i-1);
    if(i != g_nTDMProcessCount) {
      nMaxCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nTDMProcessCount)*i - 1;
    } else {
      nMaxCandidate = g_nMaxCandidate;
    }
    if(nMinPrime < MIN_PRIME) {
      nMinPrime = MIN_PRIME;
    }
    if(nMinPrime%2==0) nMinPrime++;

    if(nMinCandidate < MIN_CANDIDATE) {
      nMinCandidate = MIN_CANDIDATE;
    }
    if(nMinCandidate%2==0) nMinCandidate++;

    switch(pid) {
      case -1:
        logging(LOG_DEBUG, "fork failed, reason=failed create child process.\n");
        return;
      case 0:
        logging(LOG_DEBUG, "This is child TDM process. pid=[%ld]\n", (long)getpid());
        g_nChildProcessID = pid;
        g_nMinPrime = nMinPrime;
        g_nMaxPrime = nMaxPrime;
        g_nMinCandidate = nMinCandidate;
        g_nMaxCandidate = nMaxCandidate;
        runTDM();
        logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
        break;
      default:
        break;
    }
  }
}

void initInteger() {
  set(bigA);
  set(bigB);
  set(bigC);
  set(bigS);
  set(bigL);
  set(bigBuff);
  set(bigTemp);
  set(bigTemp2);
  logging(LOG_DEBUG, "initInteger[%ld]\n", sizeof(uint64_t)*MAX_INTEGER);
}

void PrintArray(uint64_t * src) {
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  for(uint64_t i=nIndex; i != -1; i--) {
    logging(LOG_INFO, "[%llu]", src[i]);
  }
  logging(LOG_INFO, "\n");
}

void PrintBits(uint64_t * src) {
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  for(uint64_t i=0; i <= nIndex; i++) {
    uint64_t value = src[i];
    logging(LOG_INFO, "[");
    for(uint64_t j = 0; j < 64; j++) {
      logging(LOG_INFO, "%llu", value & 0x01);
      value = value >> 1;
    }
    logging(LOG_INFO, "]\n");
  }
  logging(LOG_INFO, "\n");
}

// Lucas-Lehmer Primality Testing
// return 1 : Prime
// return 0 : NonPrime
bool PrimalityTesting(uint64_t p) {
  time_t tStepStart = 0;
  time_t tStepEnd = 0;
  set(bigS);
  addop(4, bigS);

  setMersenneNumber(p);
  for(uint64_t i = 1; i <= p -2; i++) {
    tStepStart = time(NULL);

    if(g_nAlgorithm == ALGO_MULMOD) {
      LLTmulmod(bigC, bigS);
    } else {
      LLTmethod(bigC, bigS);
    }
    tStepEnd = time(NULL);

    logging(LOG_DEBUG, "LLT [%7llu]  i=[%7llu] s=[%7llu] s0=[%20llu] t=[%3ld]\n", p, i, getBits(bigS), bigS[0], tStepEnd - tStepStart);
    logging(LOG_INFO, "LLT [%7llu]  i=[%7llu] t=[%3ld]\n", p, i, tStepEnd - tStepStart);

    if(i == p-2 && iszero(bigS)) {
      return true;
    } else if(i == p-2 && !iszero(bigS)) {
      return false;
    }
  }
  return false;
}

uint64_t getBits(uint64_t * arr) {
  uint64_t nIndex = getDigits(arr);
  uint64_t value = arr[nIndex];
  int nBits = 0;
  for(uint64_t i=0; i<64; i++) {
    if((1 & (value >> i)) == 1) {
      nBits = i;
    }
  }
  return nIndex*32 + nBits;
}

uint64_t getBitsop(uint64_t value) {
  int nBits = 0;
  for(uint64_t i=0; i<64; i++) {
    if((1 & (value >> i)) == 1) {
      nBits = i;
    }
  }
  return nBits;
}

void LLTmulmod(uint64_t * src, uint64_t * dst) {
  set(bigA);
  set(bigB);
  copy(dst, bigA);
  copy(dst, bigB);
  set(dst);
  mul(bigA, bigB, dst);
  subtractop(2, dst);
  mod(bigC,dst);
}

void LLTmethod(uint64_t * src, uint64_t * dst) {
  set(bigA);
  uint64_t nPartition = g_nCandidate + (4 - g_nCandidate%4);
  uint64_t nSubPart = nPartition/4;
  if(compare(dst, src) > 0) {
    set(dst);
    logging(LOG_ERROR,"LLTmethod error!\n");
    return;
  }

  // x < sqrt(p)
  uint64_t nIndexS = getBits(src);
  uint64_t nIndexD = getBits(dst);
  if(nIndexD*2 < nIndexS) {
    copy(dst, bigA);
    square(bigA, dst);
    if(compareop(2,dst) <= 0) {
      subtractop(2, dst);
    }
    return;
  }

  set(bigBuff);

  // A^2
  set(bigA);
  set(bigL);
  copyBits(dst, bigA, nSubPart*2,  nSubPart*4 -1);
  square(bigA, bigL);
  if(nSubPart*4 > g_nCandidate) {
    leftshiftop(nSubPart*4 - g_nCandidate, bigL);
  }
  add(bigL, bigBuff);

  // A1*B1
  set(bigA);
  set(bigB);
  set(bigL);
  copyBits(dst, bigA, nSubPart*3,  nSubPart*4 -1);
  copyBits(dst, bigB, nSubPart, nSubPart*2 -1);
  mul(bigA, bigB, bigL);
  if(nSubPart*4 + 1 > g_nCandidate) {
    leftshiftop(nSubPart*4 - g_nCandidate + 1, bigL);
  } else {
    leftshiftop(nSubPart*4 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A1*B2
  set(bigA);
  set(bigB);
  set(bigL);
  copyBits(dst, bigA, nSubPart*3,  nSubPart*4 -1);
  copyBits(dst, bigB, 0,  nSubPart -1);
  mul(bigA, bigB, bigL);
  if(nSubPart*3 + 1 > g_nCandidate) {
    leftshiftop(nSubPart*3 - g_nCandidate + 1, bigL);
  } else {
    leftshiftop(nSubPart*3 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A2*B1
  set(bigA);
  set(bigB);
  set(bigL);
  copyBits(dst, bigA, nSubPart*2,  nSubPart*3 -1);
  copyBits(dst, bigB, nSubPart, nSubPart*2 -1);
  mul(bigA, bigB, bigL);
  if(nSubPart*3 + 1 > g_nCandidate) {
    leftshiftop(nSubPart*3 - g_nCandidate + 1, bigL);
  } else {
    leftshiftop(nSubPart*3 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A2*B2
  set(bigA);
  set(bigB);
  set(bigL);
  copyBits(dst, bigA, nSubPart*2,  nSubPart*3 -1);
  copyBits(dst, bigB, 0,  nSubPart -1);
  mul(bigA, bigB, bigL);
  if(nSubPart*2 + 1 > g_nCandidate) {
    leftshiftop(nSubPart*2 - g_nCandidate + 1, bigL);
  } else {
    leftshiftop(nSubPart*2 + 1, bigL);
  }
  add(bigL, bigBuff);

  // B^2
  set(bigB);
  set(bigL);
  copyBits(dst, bigB, 0,  nSubPart*2 -1);
  square(bigB, bigL);
  add(bigL, bigBuff);
  mod(src,bigBuff);
  copy(bigBuff, dst);

  if(compareop(2, dst) <= 0) {
    subtractop(2, dst);
  }
}

void mul(uint64_t * src1, uint64_t * src2, uint64_t * dst) {
  if(g_nUseFFT == USE_FFT) {
    mul_fft(src1,src2,dst);
    return;
  }
  if(g_nUseOMP == USE_OMP) {
    mul_omp(src1,src2,dst);
    return;
  }
  uint64_t nIndexS1 = getDigits(src1);
  uint64_t nIndexS2 = getDigits(src2);
  uint64_t value = 0;
  uint64_t exponent = 0;

  set(dst);

  for(uint64_t i = 0; i <= nIndexS1; i++) {
    for(uint64_t j = 0; j <= nIndexS2; j++) {
      value = src1[i]*src2[j];
      exponent = i + j;
      dst[exponent] += ((value << 32) >> 32);
      propa_carrier(exponent, dst);
      dst[exponent+1] += value >> 32;
      propa_carrier(exponent+1, dst);
    }
  }
}

void mul_omp(uint64_t * src1, uint64_t * src2, uint64_t * dst) {
  uint64_t nIndexS1 = getDigits(src1);
  uint64_t nIndexS2 = getDigits(src2);
  uint64_t value = 0;
  uint64_t exponent = 0;

  set(dst);

  //#pragma omp parallel for num_threads(10)
  #pragma omp parallel for shared(dst)
  for(uint64_t i = 0; i <= nIndexS1; i++) {
    for(uint64_t j = 0; j <= nIndexS2; j++) {
      value = src1[i]*src2[j];
      exponent = i + j;
      dst[exponent] += ((value << 32) >> 32);
      propa_carrier(exponent, dst);
      dst[exponent+1] += value >> 32;
      propa_carrier(exponent+1, dst);
    }
  }
}

void square(uint64_t * src, uint64_t *dst) {
  if(g_nUseFFT == USE_FFT) {
    square_fft(src, dst);
    return;
  }
  if(g_nUseOMP == USE_OMP) {
    //square_omp(src, dst);
    return;
  }
  mul(src, src, dst);
}

void init_bigint_from_arr(BigInt *num, const uint64_t *arr, int length) {
  init_bigint(num); // Initialize the BigInt
  array_to_digits(arr, length, num->digits);
  for(int i = MAX_DIGITS -1 ; i >= 0; i--) {
    if(num->digits[i] != 0) {
      num->length = i + 2;
      break;
    }
  }
}

void array_from_bigint(uint64_t *arr,const BigInt *num) {
  digits_to_array(num->digits, num->length, arr);
}

void digits_to_array(const int digits[], int num_digits, uint64_t* result) {
    int nChunk = SIZE_INPUT_CHUNK;
    int index = nChunk - 1; // Start storing from the most significant end of the result array

    // Initialize the result array with zeros
    for (int i = 0; i < nChunk; i++) {
        result[i] = 0;
    }

    for (int i = 0; i < num_digits; i++) {
        // Multiply the current result by 10 and add the current digit
        for (int j = nChunk - 1; j >= index; j--) {
            result[j] *= BASE;
        }
        result[index] += digits[i];

        // Perform carry propagation if necessary
        for (int j = nChunk - 1; j > 0; j--) {
            if (result[j] >= bits_max) { // If the current element exceeds the base (2^32), propagate the carry
                result[j - 1] += (result[j] / (bits_max));
                result[j] %= (bits_max);
            }
        }
    }
    for(int i = 0; i < nChunk/2; i++) {
      if(result[i] != result[nChunk-i-1]) {
        uint64_t temp = 0;
        temp = result[i];
        result[i] = result[nChunk-i-1];
        result[nChunk-i-1] = temp;
      }
    }
}

void array_to_digits(const uint64_t array[], int num_elements, int digits[]) {
  int nChunk = SIZE_OUTPUT_CHUNK;
  uint64_t arr[nChunk];
  memset(arr, 0x00, sizeof(uint64_t)*nChunk);
  memcpy(arr, array, sizeof(uint64_t)*nChunk);
  // Start from the most significant end of the array
  int index = nChunk-1;

  // Initialize the digits array with zeros
  for (int i = 0; i < MAX_DIGITS; i++) {
    digits[i] = 0;
  }

  while (arr[index] <= 0) index--;
  if(index <= 0) {
    return;
  }
  int dpos = 0;
  int remainder = 0;
  // Loop through the array and extract digits
  for (int i = 0; i < num_elements; i++) {
    // Extract each digit by performing modulo 10
    while (arr[index] > 0) {
      if(index != 0) {
        remainder = arr[index] % BASE;
        arr[index]-=remainder;
      }
      // max int : 4294967296
      char temp[10];
      int pos = 0;
      memset(temp, 0x00, 10);
      while(arr[index] > 0) {
        temp[pos] = arr[index] % BASE;
        arr[index] /= BASE;
        pos++;
      }
      for(int j = pos -1; j >= 0; j--) {
        digits[dpos++] = temp[j];
      }
      if(index != 0) {
        arr[index-1] += remainder*bits_max;
      }
    }
    index--;
  }
  for(int i = 0; i < dpos/2; i++) {
    if(digits[i] != digits[dpos-i-1]) {
      uint64_t temp = 0;
      temp = digits[i];
      digits[i] = digits[dpos-i-1];
      digits[dpos-i-1] = temp;
    }
  }
}

void mul_fft(uint64_t * src1, uint64_t * src2, uint64_t * dst) {
  if(g_nUseOMP == USE_OMP) {
    //mul_fft_omp(src1,src2,dst);
    return;
  }
  BigInt bigA, bigB, bigResult;
  uint64_t temp[SIZE_INPUT_CHUNK];
  int nChunk = SIZE_OUTPUT_CHUNK;
  int nCountA = (getDigits(src1)+1)/nChunk + 1;
  int nCountB = (getDigits(src2)+1)/nChunk + 1; 

  for(int i = 0; i < nCountA; i++) {
    for(int j = 0; j < nCountB; j++) {
      if(i == nCountA - 1) {
        init_bigint_from_arr(&bigA, src1+i*nChunk*8, (getDigits(src1)+1)%nChunk);
      } else {
        init_bigint_from_arr(&bigA, src1+i*nChunk*8, nChunk);
      }
      if(j == nCountB - 1) {
        init_bigint_from_arr(&bigB, src2+j*nChunk*8, (getDigits(src2)+1)%nChunk);
      } else {
        init_bigint_from_arr(&bigB, src2+j*nChunk*8, nChunk);
      }
      init_bigint(&bigResult);
      mul_fft_internal(bigA.digits, bigB.digits, bigResult.digits);
      bigResult.length = getDigitsBased10(bigResult.digits);
      memset(temp, 0x00, sizeof(uint64_t)*nChunk);
      array_from_bigint(temp, &bigResult);
      for(int k = 0; k < SIZE_INPUT_CHUNK + 1; k++) {
        dst[i*nChunk+k] += temp[k];
      }
    }
  }
  for(int i = 0; i < getDigits(dst) + 1; i++) {
    if(dst[i] >= bits_max) {
      dst[i+1]+= dst[i] / bits_max;
      dst[i]= dst[i] % bits_max;
    }
  }
}

void mul_fft_internal(int * src1,int * src2, int * dst) {
  int size = getDigitsBased10(src1) + getDigitsBased10(src2) +1;

  memset(x, 0x00, sizeof(double)*MAX_DIGITS);
  memset(xi, 0x00, sizeof(double)*MAX_DIGITS);
  memset(rx, 0x00, sizeof(double)*MAX_DIGITS);
  memset(rxi, 0x00, sizeof(double)*MAX_DIGITS);
  memset(y, 0x00, sizeof(double)*MAX_DIGITS);
  memset(yi, 0x00, sizeof(double)*MAX_DIGITS);
  memset(ry, 0x00, sizeof(double)*MAX_DIGITS);
  memset(ryi, 0x00, sizeof(double)*MAX_DIGITS);

  int n = 1;
  while( n < size) {
    n <<= 1;
  }
  for (int i = 0; i < n; ++i) {
    int s = 0;
    for (int b = 1, d = n / 2; b < n; b <<= 1, d >>= 1) {
      if( b & i) {
        s += d;
      }
    }
    rx[s] = i < size ? src1[i] : 0;
    rxi[s] = 0;
  }


  for (int mte = 2; mte <= n; mte <<= 1) {
    double we = cos(2 * PI / mte);
    double wei = sin(2 * PI / mte);
    for (int i = 0; i < n; i += mte) {
      double w = 1;
      double wi = 0;
      for (int j = i; j < i + mte / 2; ++j) {
        int j2 = j + mte /2;
        double temp1R = rx[j];
        double temp1I = rxi[j];
        double temp2R = rx[j2];
        double temp2I = rxi[j2];
        rx[j] = temp1R + w * temp2R - wi * temp2I;
        rxi[j] = temp1I + wi * temp2R + w * temp2I;
        rx[j2] = temp1R - (w * temp2R - wi * temp2I);
        rxi[j2] = temp1I - (wi * temp2R + w * temp2I);
        double tmp = w * we - wi * wei;
        wi = wi * we + w * wei;
        w = tmp;
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    int s = 0;
    for (int b = 1, d = n / 2; b < n; b <<= 1, d >>= 1) {
      if( b & i) {
        s += d;
      }
    }
    ry[s] = i < size ? src2[i] : 0;
    ryi[s] = 0;
  }

  for (int mte = 2; mte <= n; mte <<= 1) {
    double we = cos(2 * PI / mte);
    double wei = sin(2 * PI / mte);
    for (int i = 0; i < n; i += mte) {
      double w = 1;
      double wi = 0;
      for (int j = i; j < i + mte / 2; ++j) {
        int j2 = j + mte /2;
        double temp1R = ry[j];
        double temp1I = ryi[j];
        double temp2R = ry[j2];
        double temp2I = ryi[j2];
        ry[j] = temp1R + w * temp2R - wi * temp2I;
        ryi[j] = temp1I + wi * temp2R + w * temp2I;
        ry[j2] = temp1R - (w * temp2R - wi * temp2I);
        ryi[j2] = temp1I - (wi * temp2R + w * temp2I);
        double tmp = w * we - wi * wei;
        wi = wi * we + w * wei;
        w = tmp;
      }
    }
  }

  memset(x, 0x00, sizeof(double)*MAX_DIGITS);
  memset(xi, 0x00, sizeof(double)*MAX_DIGITS);

  for (int i = 0; i < size; ++i) {
    x[i] = rx[i]*ry[i] - rxi[i]*ryi[i];
    xi[i] = rxi[i]*ry[i] + rx[i]*ryi[i];
  }

  for (int i = 0; i < n; ++i) {
    int s = 0;
    for (int b = 1, d = n / 2; b < n; b <<= 1, d >>= 1) {
      if( b & i) {
        s += d;
      }
    }
    rx[s] = i < size ? x[i] : 0;
    rxi[s] = i < size ? xi[i] : 0;
  }

  for (int mte = 2; mte <= n; mte <<= 1) {
    double we = cos(-2 * PI / mte);
    double wei = sin(-2 * PI / mte);
    for (int i = 0; i < n; i += mte) {
      double w = 1;
      double wi = 0;
      for (int j = i; j < i + mte / 2; ++j) {
        int j2 = j + mte /2;
        double temp1R = rx[j];
        double temp1I = rxi[j];
        double temp2R = rx[j2];
        double temp2I = rxi[j2];
        rx[j] = temp1R + w * temp2R - wi * temp2I;
        rxi[j] = temp1I + wi * temp2R + w * temp2I;
        rx[j2] = temp1R - (w * temp2R - wi * temp2I);
        rxi[j2] = temp1I - (wi * temp2R + w * temp2I);
        double tmp = w * we - wi * wei;
        wi = wi * we + w * wei;
        w = tmp;
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    rx[i] /= n;
    rxi[i] /= n;
  }

  for(int i = 0; i < size; i++) {
    int value = dst[i] + (int)llround(rx[i]);
    dst[i] = value % BASE;
    dst[i+1] = value / BASE;
  }
}

void square_fft(uint64_t * src, uint64_t * dst) {
  if(g_nUseOMP == USE_OMP) {
    //square_fft_omp(src1,src2,dst);
    return;
  }
  BigInt bigA, bigB, bigResult;
  uint64_t temp[SIZE_INPUT_CHUNK];
  int nChunk = SIZE_OUTPUT_CHUNK;
  int nCountA = (getDigits(src)+1)/nChunk + 1;
  int nCountB = (getDigits(src)+1)/nChunk + 1;

  assign_bigint(&bigA, &bigB);
  for(int i = 0; i < nCountA; i++) {
    for(int j = 0; j < nCountB; j++) {
      if(i == nCountA - 1) {
        init_bigint_from_arr(&bigA, src+i*nChunk*8, (getDigits(src)+1)%nChunk);
      } else {
        init_bigint_from_arr(&bigA, src+i*nChunk*8, nChunk);
      }
      if(j == nCountB - 1) {
        init_bigint_from_arr(&bigB, src+j*nChunk*8, (getDigits(src)+1)%nChunk);
      } else {
        init_bigint_from_arr(&bigB, src+j*nChunk*8, nChunk);
      }
      init_bigint(&bigResult);
      square_fft_internal(bigA.digits, bigResult.digits);
      bigResult.length = getDigitsBased10(bigResult.digits);
      memset(temp, 0x00, sizeof(uint64_t)*nChunk);
      array_from_bigint(temp, &bigResult);
      for(int k = 0; k < SIZE_INPUT_CHUNK + 1; k++) {
        dst[i*nChunk+k] += temp[k];
      }
    }
  }
  for(int i = 0; i < getDigits(dst) + 1; i++) {
    if(dst[i] >= bits_max) {
      dst[i+1]+= dst[i] / bits_max;
      dst[i]= dst[i] % bits_max;
    }
  }
}

void square_fft_internal(int * src, int * dst) {
  int size = getDigitsBased10(src)*2 + 1;
  memset(x, 0x00, sizeof(double)*MAX_DIGITS);
  memset(xi, 0x00, sizeof(double)*MAX_DIGITS);
  memset(rx, 0x00, sizeof(double)*MAX_DIGITS);
  memset(rxi, 0x00, sizeof(double)*MAX_DIGITS);

  int n = 1;
  while( n < size) {
    n <<= 1;
  }
  for (int i = 0; i < n; ++i) {
    int s = 0;
    for (int b = 1, d = n / 2; b < n; b <<= 1, d >>= 1) {
      if( b & i) {
        s += d;
      }
    }
    rx[s] = i < size ? src[i] : 0;
    rxi[s] = 0;
  }

  for (int mte = 2; mte <= n; mte <<= 1) {
    double we = cos(2 * PI / mte);
    double wei = sin(2 * PI / mte);
    for (int i = 0; i < n; i += mte) {
      double w = 1;
      double wi = 0;
      for (int j = i; j < i + mte / 2; ++j) {
        int j2 = j + mte /2;
        double temp1R = rx[j];
        double temp1I = rxi[j];
        double temp2R = rx[j2];
        double temp2I = rxi[j2];
        rx[j] = temp1R + w * temp2R - wi * temp2I;
        rxi[j] = temp1I + wi * temp2R + w * temp2I;
        rx[j2] = temp1R - (w * temp2R - wi * temp2I);
        rxi[j2] = temp1I - (wi * temp2R + w * temp2I);
        double tmp = w * we - wi * wei;
        wi = wi * we + w * wei;
        w = tmp;
      }
    }
  }

  memset(x, 0x00, sizeof(double)*MAX_DIGITS);
  memset(xi, 0x00, sizeof(double)*MAX_DIGITS);

  for (int i = 0; i < size; ++i) {
    x[i] = rx[i]*rx[i] - rxi[i]*rxi[i];
    xi[i] = rxi[i]*rx[i] + rx[i]*rxi[i];
  }

  for (int i = 0; i < n; ++i) {
    int s = 0;
    for (int b = 1, d = n / 2; b < n; b <<= 1, d >>= 1) {
      if( b & i) {
        s += d;
      }
    }
    rx[s] = i < size ? x[i] : 0;
    rxi[s] = i < size ? xi[i] : 0;
  }

  for (int mte = 2; mte <= n; mte <<= 1) {
    double we = cos(-2 * PI / mte);
    double wei = sin(-2 * PI / mte);
    for (int i = 0; i < n; i += mte) {
      double w = 1;
      double wi = 0;
      for (int j = i; j < i + mte / 2; ++j) {
        int j2 = j + mte /2;
        double temp1R = rx[j];
        double temp1I = rxi[j];
        double temp2R = rx[j2];
        double temp2I = rxi[j2];
        rx[j] = temp1R + w * temp2R - wi * temp2I;
        rxi[j] = temp1I + wi * temp2R + w * temp2I;
        rx[j2] = temp1R - (w * temp2R - wi * temp2I);
        rxi[j2] = temp1I - (wi * temp2R + w * temp2I);
        double tmp = w * we - wi * wei;
        wi = wi * we + w * wei;
        w = tmp;
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    rx[i] /= n;
    rxi[i] /= n;
  }

  memset(dst, 0x00, sizeof(int)*MAX_DIGITS);

  // compute carry
  for(int i = 0; i < size; i++) {
    int value = dst[i] + (int)llround(rx[i]);
    dst[i] = value % BASE;
    dst[i+1] = value / BASE;
  }
}

void mulop(uint64_t op, uint64_t * dst) {
  uint64_t nIndex = getDigits(dst);
  uint64_t value = 0;
  uint64_t exponent = 0;
  set(bigTemp);
  for(uint64_t i = 0; i <= nIndex; i++) {
    value = op*dst[i];
    exponent = i;
    bigTemp[exponent + 1] += value >> 32;
    bigTemp[exponent] += (value <<32) >> 32;
    if(bigTemp[exponent] > pow(2, 32) -1) {
      bigTemp[exponent] -= pow(2,32);
      if(bigTemp[exponent + 1] == pow(2,32) -1) {
        bigTemp[exponent + 2] += 1;
        bigTemp[exponent + 1] = 0;
      } else {
        bigTemp[exponent + 1] += 1;
      }
    }
  }
  copy(bigTemp, dst);
}

void leftshift(uint64_t * src) {
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  for(uint64_t i = nIndex+1; i != -1; i--) {
    src[i] = src[i-1];
  } 
  src[0] = 0;
}

void leftshiftop(uint64_t nBits, uint64_t * src) {
  uint64_t nIndex = 0;

  if(nBits + getBits(src) > MAX_INTEGER) {
    set(src);
    return;
  }
  set(bigTemp);
  if(nBits >= 32) {
    nIndex = getDigits(src);
    if(nIndex + nBits/32 < MAX_INTEGER) {
      for(uint64_t i = nIndex; i != -1; i--) {
        src[i + nBits/32] = src[i];
        src[i] = 0;
        if(i == 0) {
          for(uint64_t j = i + nBits/32 -1; j > 0; j--) {
            src[j] = 0;
          }
        }
      }
    } else {
      logging(LOG_ERROR,"leftshift overflow\n");
      exit(0);
    }
  }
  uint64_t shift = nBits%32;
  if(shift != 0) {
    nIndex = getDigits(src);
    if(nIndex + 1 < MAX_INTEGER) {
      for(uint64_t i = nIndex; i != -1; i--) {
        if(i!=0) {
          uint64_t highvalue = src[i];
          uint64_t lowvalue = src[i-1];
        
          if(i == nIndex) {
            bigTemp[i+1] = highvalue >> (32 - shift);
          }
          bigTemp[i] = (highvalue << (shift+32)) >> 32; 
          bigTemp[i]+= (lowvalue << shift) >> 32;
        } else {
          uint64_t lowvalue = src[0];
          if(nIndex == 0) {
            bigTemp[1] = (lowvalue << shift) >> 32;
          }
          bigTemp[0] = (lowvalue << (shift+32)) >> 32;
        }
      }
      copy(bigTemp, src);
    }
  }
}

void rightshift(uint64_t * src) {
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  for(uint64_t i = 1; i <= nIndex; i++) {
    src[i-1] = src[i];
    if(nIndex == i) {
      src[i] = 0;
    }
  } 
}

void rightshiftop(uint64_t nBits, uint64_t * src) {
  uint64_t nIndex = 0;
  if(nBits > getBits(src)) {
    set(src); 
    return; 
  }
  set(bigTemp);
  if(nBits >= 32) {
    nIndex = getDigits(src);
    if(nIndex - nBits/32 >= 0) {
      for(uint64_t i = 0; i <= nIndex -nBits/32; i++) {
        src[i] = src[i+nBits/32];
        src[i+nBits/32] = 0;
        if(i == 0) {
          for(uint64_t j = i + nBits/32 -1; j > 0; j--) {
            src[j] = 0;
          }
        }
      }
    } else {
      logging(LOG_ERROR,"rightshift underflow\n");
      exit(0);
    }
  }

  uint64_t shift = nBits%32;
  if(shift != 0) {
    nIndex = getDigits(src);
    if(nIndex + 1 < MAX_INTEGER) {
      for(uint64_t i = 0; i <= nIndex; i++) {
        uint64_t lowvalue = src[i];
        uint64_t highvalue = src[i+1];
        
        bigTemp[i]= lowvalue >> shift; 
        bigTemp[i]+= (highvalue << (64 - shift)) >> 32;
      }
      copy(bigTemp, src);
    }
  }
}

void propa_carrier(uint64_t nIndex, uint64_t * dst) {
  for(uint64_t i = nIndex; i < MAX_INTEGER; i++) {
    if(dst[i] < bits_max -1) {
      return;
    } else {
      while(dst[i] >= bits_max) {
        dst[i] = dst[i] - bits_max;
        dst[i+1]+= 1;
      }
      propa_carrier(i+1, dst);
      break;
    }
  }
}

void add(uint64_t * src, uint64_t * dst) {
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  for(uint64_t i = 0; i <= nIndex; i++) {
    uint64_t value = src[i] + dst[i];
    dst[i] = (value << 32) >> 32;
    dst[i+1] += value >> 32;
    propa_carrier(i+1, dst); 
  }
}

void addop(uint64_t op, uint64_t * dst) {
  set(bigTemp);
  bigTemp[0] = op;
  add(bigTemp, dst);
}

void set(uint64_t * dst) {
  memset(dst, 0x00, sizeof(uint64_t)*MAX_INTEGER);
}

void subtract(uint64_t * src, uint64_t * dst) {
  // The value returned by subtract must be unsigned integer.
  uint64_t nIndexS = 0;
  uint64_t nIndexD = 0;

  nIndexS = getDigits(src);
  nIndexD = getDigits(dst);

  if(compare(src, dst) <= 0) {
    for(uint64_t i = 0; i <= nIndexS; i++) {
      if(dst[i] >= src[i]) {
        dst[i] -= src[i];
      } else {
        for(uint64_t j = i+1; j <= nIndexD; j++) {
          if(dst[j] != 0) {
            dst[j] -=1;
            for(uint64_t k = j-1; k > i && k != -1; k--) {
              dst[k] = pow(2,32) -1; 
            }
            dst[i] = (dst[i] + pow(2,32)) - src[i];
            break;
          }
        }
      }
    }
  } else {
    logging(LOG_ERROR,"subtract error\n");
    exit(0);
  }
}

void subtractop(uint64_t op, uint64_t * dst) {
  set(bigBuff);
  bigBuff[0] = op;
  subtract(bigBuff, dst);
}

void mod(uint64_t * src, uint64_t * dst) {
  while(compare(src,dst) <= 0) {
    if(compare(src,dst) == 0) {
      set(dst);
      return;
    }
    uint64_t nIndexS = getBits(src);
    uint64_t nIndexD = getBits(dst);
    set(bigA);
    set(bigB);

    copyBits(dst, bigA, nIndexS+1,  nIndexD+1);
    copyBits(dst, bigB, 0, nIndexS);
    add(bigA, bigB); 
    copy(bigB, dst);
  }
}

int compareop(uint64_t op, uint64_t * dst) {
  set(bigTemp2);
  addop(2, bigTemp2);
  return compare(bigTemp2, dst);
}

int compare(uint64_t * src, uint64_t * dst) {
  uint64_t nIndexS = 0;
  uint64_t nIndexD = 0;
  nIndexS = getBits(src);
  nIndexD = getBits(dst);

  if(nIndexS>nIndexD) {
    return 1;
  } else if(nIndexS<nIndexD){
    return -1; 
  }
  for(uint64_t k = nIndexD; k != -1; k--) {
    if(src[k] > dst[k]) {
      return 1; 
    } else if(src[k] < dst[k]) {
      return -1;
    }
  }
  return 0;
}

void copy(uint64_t * src, uint64_t * dst) {
  set(dst);
  uint64_t nIndex = 0;
  nIndex = getDigits(src);
  if(nIndex >= 0 && nIndex < MAX_INTEGER) {
    for(uint64_t i=nIndex; i!= -1; i--) {
      dst[i] = src[i];
    }
  } else {
    logging(LOG_ERROR,"copy error\n");
    exit(0);
  }
}

void copyBits(uint64_t * src, uint64_t * dst, uint64_t nStart, uint64_t nEnd) {
  set(dst);
  copy(src, dst);

  rightshiftop(nStart, dst);
  uint64_t nEndBits = (nEnd - nStart + 1)%32;
  uint64_t nEndIndex = (nEnd - nStart)/32;
  uint64_t value = dst[nEndIndex];
  dst[nEndIndex] = ((value << (64 - nEndBits)) >> (64 - nEndBits));
  uint64_t nIndexS = getDigits(src);
  if(nIndexS - nEndIndex + 1 != -1) {
    memset(dst + nEndIndex + 1, 0x00, sizeof(uint64_t)* (nIndexS - nEndIndex + 1));
  }
}

uint64_t getDigits(uint64_t * arr) {
  for(uint64_t i = MAX_INTEGER-1; i != -1; i--) {
    if(arr[i] != 0) {
      return i;
    }
  }  
  return 0; 
}

int getDigitsBased10(int * arr) {
  for(int i = MAX_DIGITS-1; i != -1; i--) {
    if(arr[i] != 0) {
      return i;
    }
  }  
  return 0; 
}

void setMersenneNumber(uint64_t p) {
  set(bigC);
  int qt = (int)((p-1)/32);
  for(uint64_t i = 0; i < qt; i++) {
    bigC[i] = UINT32_MAX;
  }
  int rm = (int)((p-1)%32);
  for(uint64_t i = 0; i <= rm; i++) {
    bigC[qt] += (uint64_t)pow(2, i);
  }
}

bool iszero(uint64_t * arr) {
  for(int i = 0; i < MAX_INTEGER; i++) {
    if(arr[i] != 0) {
      return false;
    }
  }
  return true;
}

void error_handling(char *message) {
  fprintf(stderr,"ERROR message=[%s]\n", message);
}

void server() {
  int serv_sock;
  int clnt_sock;

  char message[MAX_SIZE];
  int len;

  struct sockaddr_in serv_addr;
  struct sockaddr_in clnt_addr;
  int clnt_addr_size;

  serv_sock = socket(PF_INET, SOCK_STREAM, 0);

  if(serv_sock == -1) {
    error_handling("socket() error");
  } else {
    memset(&serv_addr, 0, sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr=htonl(INADDR_ANY);
    serv_addr.sin_port = htons(g_nPort);

    if(bind(serv_sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr))==-1) {
      error_handling("bind() error");
    } else {
      if(listen(serv_sock, 5) == -1) {
        error_handling("listen() error");
      } else {
        clnt_addr_size = sizeof(clnt_addr);

        readDB();
        while(1) {
          clnt_sock = accept(serv_sock, (struct sockaddr*)&clnt_addr,(socklen_t *) &clnt_addr_size);
          if(clnt_sock == -1) {
            error_handling("accept() error");
            break;
          } else {
            agent(clnt_sock);
          }
        }
      } // if(listen(serv_sock, 5) == -1) {
    } // if(bind(serv_sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr))==-1) {
  } // if(serv_sock == -1) {
}

void agent(int sock) {
  int len = 0;
  char message[MAX_SIZE];
  pid_t pid;
  pid = fork();
  switch(pid) {
    case -1:
      logging(LOG_DEBUG, "fork failed, reason=failed create child process.\n");
      return;       
    case 0:
      logging(LOG_DEBUG, "This is child process. pid=[%ld]\n", (long)getpid());
      g_nChildProcessID = pid;
      memset(message, 0x00, MAX_SIZE);
      while((len=read(sock, message, MAX_SIZE)) != 0) {
        if(strncmp(message, "status", strlen("status")) == 0) {
          char response[MAX_SIZE];
          memset(response, 0x00, MAX_SIZE);
          sprintf(response, "{ \"candidatecount\" : \"%d\"}\n", getCandidateCount());
          write(sock, response, sizeof(response));
        }
        if(strncmp(message, "push", strlen("push")) == 0) {
          char response[MAX_SIZE];
          memset(response, 0x00, MAX_SIZE);
          sprintf(response, "{ \"result\" : \"successed\" }\n");
          write(sock, response, sizeof(response));
        }
        if(strncmp(message, "mersenne", strlen("mersenne")) == 0) {
          char response[MAX_SIZE];
          memset(response, 0x00, MAX_SIZE);
          sprintf(response, "{ \"mersenne\" : \"mersenne list\" }\n");
          write(sock, response, sizeof(response));
        }
        if(strncmp(message, "log", strlen("log")) == 0) {
          char response[MAX_SIZE];
          memset(response, 0x00, MAX_SIZE);
          sprintf(response, "{ \"log\" : \"log string\" }\n");
          write(sock, response, sizeof(response));
        }
        write(sock, "HTTP1.1 200 OK\n\n", strlen("HTTP1.1 200 OK\n\n"));
      }
      close(sock);
      logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
      break;
    default:
      break;
  }
}

void console() {
  if(g_nUseSocket == USE_DIRECT) {
    readBits();
    int nCount = 0;
    int nMinCandidate = 0;
    for(int i = g_nMinCandidate; i < g_nMaxCandidate; i+=2) {
      if(isCandidate(i) && !isKnownMersenne(i)) {
        nCount++;
        if(nMinCandidate == 0) {
          nMinCandidate = i;
        }
      }
    }
    logging(LOG_INFO, "Total Candidate Count : [%d] (M%llu ~ M%llu)\n", nCount, g_nMinCandidate, g_nMaxCandidate);
    logging(LOG_INFO, "Min Mersenne Candidate : [%d]\n", nMinCandidate);
    if(isSingleProcess()) {
      logging(LOG_INFO, "TDM Last Prime : [%ld]\n", getLastTDM());
    }

    for(int i = g_nDisplayStart; i < g_nDisplayEnd; i++) {
      if(isCandidate(i) && !isKnownMersenne(i)) {
        logging(LOG_DEBUG, "*"); 
      } else {
        logging(LOG_DEBUG, " ");
      }
      if(i%100 == 0) {
        logging(LOG_DEBUG, "\n");
      }
    }
    logging(LOG_DEBUG, "\n");
  } else if(g_nUseSocket == USE_SOCKET) {
    int sock;
    struct sockaddr_in serv_addr;
    char message[MAX_SIZE];
    int len;

    sock = socket(PF_INET, SOCK_STREAM, 0);
    if(sock == -1) {
      error_handling("socket() error");
    } else {
      memset(&serv_addr, 0, sizeof(serv_addr));
      serv_addr.sin_family = AF_INET;
      serv_addr.sin_addr.s_addr = inet_addr(g_achAddress);
      serv_addr.sin_port = htons(g_nPort);

      if(connect(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) == -1) {
        error_handling("connect() error");
      } else {
        while(1) {
          fgets(message, MAX_SIZE, stdin);
          write(sock, message, strlen(message));

          len = read(sock, message, MAX_SIZE-1);
          message[len] = 0;
          logging(LOG_DUMP, "%s\n", message);
        }
        close(sock);
      }
    }
  } else {
    logging(LOG_ERROR, "Unknown UseSocket Option [%d]\n", g_nUseSocket);
  }
}

void benchmark() {
  logging(LOG_INFO, "benchmark_multiply begin\n");
  benchmark_multiply();
  logging(LOG_INFO, "benchmark_multiply end\n");
}

void test() {

  logging(LOG_INFO, "test_mersenne begin\n");
  test_mersenne();
  logging(LOG_INFO, "test_mersenne end\n");

  logging(LOG_INFO, "test_get_primes begin\n");
  test_get_primes();
  logging(LOG_INFO, "test_get_primes end\n");

#ifdef USE_COMPRESSION
  logging(LOG_INFO, "test_compress begin\n");
  test_compress();
  logging(LOG_INFO, "test_compress end\n");
#endif

#ifdef USE_ENCRYPTION
  logging(LOG_INFO, "test_sign begin\n");
  test_sign();
  logging(LOG_INFO, "test_sign end\n");
#endif

#ifdef EXPERIMENTAL
  logging(LOG_INFO, "test_riemann_zeta begin\n");
  test_riemann_zeta();
  logging(LOG_INFO, "test_riemann_zeta end\n");

  logging(LOG_INFO, "test_zeta begin\n");
  test_zeta();
  logging(LOG_INFO, "test_zeta end\n");

  logging(LOG_INFO, "test_goldBach begin\n");
  test_goldBach();
  logging(LOG_INFO, "test_goldBach end\n");

  logging(LOG_INFO, "test_prime_gap_distribution begin\n");
  test_prime_gap_distribution();
  logging(LOG_INFO, "test_prime_gap_distribution end\n");

  logging(LOG_INFO, "test_geometric_series begin\n");
  test_geometric_series();
  logging(LOG_INFO, "test_geometric_series end\n");

  logging(LOG_INFO, "test_fft begin\n");
  test_fft();
  logging(LOG_INFO, "test_fft end\n");

  logging(LOG_INFO, "test_miller_rabin begin\n");
  test_miller_rabin();
  logging(LOG_INFO, "test_miller_rabin end\n");

  logging(LOG_INFO, "test_lucaslehmer begin\n");
  test_lucaslehmer();
  logging(LOG_INFO, "test_lucaslehmer end\n");

  logging(LOG_INFO, "test_aks begin\n");
  test_aks();
  logging(LOG_INFO, "test_aks end\n");

  logging(LOG_INFO, "test_montgomery_multiply begin\n");
  test_montgomery_multiply();
  logging(LOG_INFO, "test_montgomery_multiply end\n");

  logging(LOG_INFO, "test_count_primes begin\n");
  test_count_primes();
  logging(LOG_INFO, "test_count_primes end\n");

  logging(LOG_INFO, "test_fermat_primes begin\n");
  test_fermat_primes();
  logging(LOG_INFO, "test_fermat_primes end\n");

  logging(LOG_INFO, "test_nn begin\n");
  test_nn();
  logging(LOG_INFO, "test_nn end\n");

  logging(LOG_INFO, "test_wilson begin\n");
  test_wilson();
  logging(LOG_INFO, "test_winson end\n"); 

  logging(LOG_INFO, "test_monteCarloSimulation begin\n");
  test_monteCarloSimulation();
  logging(LOG_INFO, "test_montenCarloSimulation end\n"); 

#ifdef BIGINT
  logging(LOG_INFO, "test_bigint begin\n");
  test_bigint();
  logging(LOG_INFO, "test_bigint end\n"); 
#endif // #ifdef BIGINT

#endif
} 

void benchmark_multiply() {
  BigInt a, b, result;
  int n = MAX_DIGITS/2;
  int un = MAX_INTEGER/(2*32);
  time_t start = 0;
  initInteger();

  logging(LOG_DEBUG, "multiply benchmark test digits[%d] bits[%d]\n", n, un);

  set(bigA);
  set(bigB);
  set(bigS);
  for(int i = 0; i < un + 1; i++) {
    bigA[i] = rand() % bits_max;
    bigB[i] = rand() % bits_max;
  }

  g_nUseOMP = USE_NOOMP;
  g_nUseFFT = USE_NOFFT;
  start = time(NULL);
  mul(bigA,bigB,bigS);
  logging(LOG_DEBUG, "Result : ");
  printf("multiply time [%ld]\n", time(NULL) - start);

  set(bigA);
  set(bigB);
  set(bigS);

  logging(LOG_DEBUG, "multiply with OMP\n");
  g_nUseOMP = USE_OMP;
  g_nUseFFT = USE_NOFFT;
  start = time(NULL);
  mul(bigA,bigB,bigS);
  logging(LOG_DEBUG, "Result : ");
  printf("multiply OMP time [%ld]\n", time(NULL) - start);

  set(bigA);
  set(bigB);
  set(bigS);

  logging(LOG_DEBUG, "multiply with FFT\n");
  g_nUseOMP = USE_NOOMP;
  g_nUseFFT = USE_FFT;
  start = time(NULL);
  mul(bigA,bigB,bigS);
  logging(LOG_DEBUG, "Result : ");
  printf("multiply FFT time [%ld]\n", time(NULL) - start);

  set(bigA);
  set(bigB);
  set(bigS);

  logging(LOG_DEBUG, "multiply with FFT + OMP\n");
  g_nUseOMP = USE_OMP;
  g_nUseFFT = USE_FFT;
  start = time(NULL);
  mul(bigA,bigB,bigS);
  logging(LOG_DEBUG, "Result : ");
  printf("multiply FFT+OMP time [%ld]\n", time(NULL) - start);
}

void test_mersenne() {
  initInteger();

  logging(LOG_DEBUG, "addop test + 2\n");
  addop(2, bigS);
  PrintArray(bigS);

  logging(LOG_DEBUG, "subtractop test - 1\n");
  subtractop(1, bigS);
  PrintArray(bigS);

  logging(LOG_DEBUG, "add test [1][1][1]\n");
  bigBuff[2] = 1;
  bigBuff[1] = 1;
  bigBuff[0] = 1; 
  add(bigBuff, bigS);
  PrintArray(bigS);

  logging(LOG_DEBUG, "subtract test [1][0][0]\n");
  bigBuff[2] = 1;
  bigBuff[1] = 0;
  bigBuff[0] = 0; 
  subtract(bigBuff, bigS);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigS);

  logging(LOG_DEBUG, "mul test [1][1][1]*[1][1][1] -> [1][2][3][2][1]\n");
  set(bigBuff);
  set(bigS);
  bigA[2] = 1;
  bigA[1] = 1;
  bigA[0] = 1; 
  bigB[2] = 1;
  bigB[1] = 1;
  bigB[0] = 1; 
  mul(bigA,bigB,bigS);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigS);

  logging(LOG_DEBUG, "mod test [0][0][10] mod [0][0][7] -> [0][0][4]\n");
  set(bigC);
  set(bigS);
  bigS[2] = 0;
  bigS[1] = 0;
  bigS[0] = 10; 
  bigC[2] = 0; 
  bigC[1] = 0; 
  bigC[0] = 7; 
  mod(bigC, bigS);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigS);

  logging(LOG_DEBUG, "square test [1][1][1]^2 -> [1][2][3][2][1]\n");
  set(bigC);
  set(bigS);
  bigC[2] = 1; 
  bigC[1] = 1; 
  bigC[0] = 1; 
  square(bigC, bigS);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigS);

  logging(LOG_DEBUG, "leftshift test [1][1][1] -> [1][1][1][0]\n");
  set(bigC);
  bigC[2] = 1;
  bigC[1] = 1;
  bigC[0] = 1;
  leftshift(bigC);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigC);

  logging(LOG_DEBUG, "leftshiftop test 16 bit [1][1][1] -> [65536][65536][65536]\n");
  set(bigC);
  bigC[2] = 1;
  bigC[1] = 1;
  bigC[0] = 1;
  leftshiftop(16, bigC);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigC);

  logging(LOG_DEBUG, "rightshift test [1][1][1] -> [1][1]\n");
  set(bigC);
  bigC[2] = 1;
  bigC[1] = 1;
  bigC[0] = 1;
  rightshift(bigC);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigC);

  logging(LOG_DEBUG, "rightshiftop test 16 bit [1][1][1] -> [65536][65536]\n");
  set(bigC);
  bigC[2] = 1;
  bigC[1] = 1;
  bigC[0] = 1;
  rightshiftop(16, bigC);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigC);

  logging(LOG_DEBUG, "copyBits test [1][1][31] [3...4] -> [3]\n");
  set(bigC);
  bigS[2]=0;
  bigS[1]=0;
  bigS[0]=0;
  bigC[2]=1;
  bigC[1]=1;
  bigC[0]=31;
  copyBits(bigC, bigS, 3,  4);
  logging(LOG_DEBUG, "Result : ");
  PrintArray(bigS);
  logging(LOG_DEBUG, "-----------------------rightshiftop test start\n");
  set(bigA);
  bigA[0]=1;
  rightshiftop(1, bigA);
  bigA[2]=1;
  bigA[1]=1;
  bigA[0]=1;
  rightshiftop(61, bigA);
  logging(LOG_DEBUG, "-----------------------rightshiftop test end\n");
  logging(LOG_DEBUG, "----------------------mod test start\n");
  set(bigC);
  set(bigS);
  bigS[2]=3;
  bigS[1]=4294967295;
  bigS[0]=4294967295;
  bigC[1]=4294967295;
  bigC[0]=4294967295;
  PrintArray(bigS); 
  PrintArray(bigC); 
  mod(bigC, bigS);
  PrintArray(bigS);
  PrintArray(bigC);
  logging(LOG_DEBUG, "-----------------------mod test end\n");
  logging(LOG_DEBUG, "----------------------mul test start\n");
  g_nCandidate = 107;
  set(bigA);
  set(bigB);
  set(bigS);
  bigS[3]=0;
  bigS[2]=0;
  bigS[1]=0;
  bigS[0]=0;
  bigA[3]=0;
  bigA[2]=1;
  bigA[1]=1;
  bigA[0]=1;
  square(bigA, bigS);
  PrintArray(bigS);
  logging(LOG_DEBUG, "-----------------------mul test end\n");
}

void logging(int level, const char *format, ...) {
  va_list ap;
  char buf[MAX_BUFF_SIZE];

  va_start(ap, format);
  vsprintf(buf, format, ap);
  va_end(ap);

  if(g_nLogLevel >= level) {
    if(g_nLogFile == LOG_FILE) {
      writeLog(buf);
    }
    if(g_nLogNetwork == LOG_NETWORK) {
      sendLog(buf);
    }
    if(g_nLogVerbose == LOG_VERBOSE) {
      printf("%s", buf);
    }
  }
}

void initCandidate() {
  int nCount = 0;
  logging(LOG_DEBUG, "initCandidate Begin\n");
  if(g_nNoFile == USE_FILE && access(CANDIDATE_FILE_NAME, F_OK) == 0) {
    logging(LOG_DEBUG, "initCandidate From File Min[%llu]~Max[%llu]\n", g_nMinCandidate, g_nMaxCandidate);
    readBits();
  } else {
    logging(LOG_DEBUG, "initCandidate From Process Min[%llu]~Max[%llu]\n", g_nMinCandidate, g_nMaxCandidate);
    for(uint64_t i = g_nMinCandidate; i < g_nMaxCandidate; i+=2) {
      if(isPrime(i)) {
        setNumberCache(i, MERSENNE_CANDIDATE);
        nCount++;
      }
    }
    updateBits();
  }
  logging(LOG_DEBUG, "initCandidate End pid=[%ld] [%d]\n",getpid(),  nCount);
}

void findMersenne() {
  logging(LOG_INFO, "findMersenne [%d]~[%d]\n", g_nMinPrime, g_nMaxPrime);
  int nChunkCount = (g_nMaxPrime - g_nMinPrime)/g_nChunk + 1;
  for(int i = 0; i < nChunkCount; i++) {
    switch(g_nFinder) {
      case FINDER_NORMAL:
        findChunkNormal();
        break;
      case FINDER_NOMOD:
        findChunkNomod();
        break;
      case FINDER_SKIP:
        findChunkSkip();
        break;
      case FINDER_EULER:
        findChunkEuler();
        break;
      default:
        break;
    }
    updateBits();
    PrintCandidateCount();
  }
}

void findMersenneOMP() {
  logging(LOG_INFO, "findMersenne [%d]~[%d]\n", g_nMinPrime, g_nMaxPrime);
  int nChunkCount = (g_nMaxPrime - g_nMinPrime)/g_nChunk + 1;
  #pragma omp parallel for num_threads(10)
  for(int i = 0; i < nChunkCount; i++) {
    switch(g_nFinder) {
      case FINDER_NORMAL:
        findChunkNormal();
        break;
      case FINDER_NOMOD:
        findChunkNomod();
        break;
      case FINDER_SKIP:
        findChunkSkip();
        break;
      default: 
        break;
    }
    PrintCandidateCount();
  } 
}

void findChunkNormal() {
  uint64_t chunk[g_nChunk][PRIME_ID];
  memset(chunk, 0x00, PRIME_ID*g_nChunk*sizeof(uint64_t));
  int count = getNextChunk(chunk);
  int nCount = 0;
  struct Node* pNextCandidate = g_lstCandidate;
  while(pNextCandidate) {
    int exponent = pNextCandidate->data;
    int oldexponent = 0; 
    for(int i =0; i < exponent - oldexponent; i++) {
      for(int j=0; j<count; j++) {
        if(chunk[j][PRIME_ID_USE] == 0) continue;
        chunk[j][PRIME_ID_VALUE] <<= 1;
        if(chunk[j][PRIME_ID_VALUE] < chunk[j][PRIME_ID_PRIME]) continue;
        chunk[j][PRIME_ID_VALUE] -= chunk[j][PRIME_ID_PRIME];
        if(chunk[j][PRIME_ID_VALUE] == 1) setNumberCache(pNextCandidate->data, MERSENNE_NUMBER);
        if(chunk[j][PRIME_ID_VALUE] == 0 || chunk[j][PRIME_ID_VALUE] == 1) {
          chunk[j][PRIME_ID_USE] = 0;
        }
      }
    }
    nCount++;
    if(nCount % 10000 == 0) {
      int prime = 0;
      for(int i = 0; i < count; i++) {
        if(chunk[i][PRIME_ID_USE] != 0) {
          prime++;
        }
      }
      if(prime == 0) {
        //logging(LOG_DEBUG, "findChunkSkip no live prime End\n");
        return;
      }
      //logging(LOG_DEBUG, "findChunkSkip [%llu] Candidate Count [%d] live[%d]\n", pNextCandidate->data, nCount, prime);
    }
    oldexponent = exponent;
    pNextCandidate=pNextCandidate->next;
  }
}

void findChunkNomod() {
  uint64_t chunk[g_nChunk][PRIME_ID];
  memset(chunk, 0x00, PRIME_ID*g_nChunk*sizeof(uint64_t));
  int count = getNextChunk(chunk);
  int nCount = 0;
  struct Node* pNextCandidate = g_lstCandidate;
  while(pNextCandidate) {
    int exponent = pNextCandidate->data;
    int oldexponent = 0;
    for(int i =0; i < exponent - oldexponent; i++) {
      for(int j=0; j<count; j++) {
        if(chunk[j][PRIME_ID_USE] == 0) continue;
        chunk[j][PRIME_ID_VALUE] <<= 1;
        chunk[j][PRIME_ID_VALUE] %= chunk[j][PRIME_ID_PRIME];
        if(chunk[j][PRIME_ID_VALUE] == 1) setNumberCache(pNextCandidate->data, MERSENNE_NUMBER);
        if(chunk[j][PRIME_ID_VALUE] == 0 || chunk[j][PRIME_ID_VALUE] == 1) {
          chunk[j][PRIME_ID_USE] = 0;
        }
      }
    }
    nCount++;
    if(nCount % 10000 == 0) {
      int prime = 0;
      for(int i = 0; i < count; i++) {
        if(chunk[i][PRIME_ID_USE] != 0) {
          prime++;
        }
      }
      if(prime == 0) {
        //logging(LOG_DEBUG, "findChunkSkip no live prime End\n");
        return;
      }
      //logging(LOG_DEBUG, "findChunkSkip [%llu] Candidate Count [%d] live[%d]\n", pNextCandidate->data, nCount, prime);
    }
    oldexponent = exponent;
    pNextCandidate=pNextCandidate->next;
  }
}

void findChunkSkip() {
  uint64_t chunk[g_nChunk][PRIME_ID];
  memset(chunk, 0x00, PRIME_ID*g_nChunk*sizeof(uint64_t));
  logging(LOG_DEBUG, "findChunkSikp Start\n");
  int count = getNextChunk(chunk);
  //logging(LOG_DEBUG, "findChunkSikp ChunkCount[%d]\n", count);
  int nCount = 0;
  struct Node* pNextCandidate = g_lstCandidate;
  while(pNextCandidate) {
    int exponent = pNextCandidate->data;
    int oldexponent = 0;
    int jump = (exponent - oldexponent);
    int unit = 64 - getBitsop(chunk[count-1][PRIME_ID_PRIME])-1;

    for(int i=0; i < jump/unit + 1; i++) {
      for(int j=0; j<count; j++) {
        if(chunk[j][PRIME_ID_USE] == 0) continue;
        int shift = 0;
        if(i == jump/unit) {
          shift = (jump%unit);
        } else {
          shift = unit;
        }
        chunk[j][PRIME_ID_VALUE] <<= shift;
        chunk[j][PRIME_ID_VALUE] %= chunk[j][PRIME_ID_PRIME];
        uint64_t comparer = 1;
        for(int k=0; k < shift; k++) {
          if(chunk[j][PRIME_ID_VALUE] == comparer) {
            setNumberCache(pNextCandidate->data, MERSENNE_NUMBER);
            //logging(LOG_DEBUG, "findChunkSikp setNumberCache[%llu] prime[%llu]\n", pNextCandidate->data, chunk[j][PRIME_ID_PRIME]);
            chunk[j][PRIME_ID_USE] = 0;
          }
          comparer <<= 1;
        }
      }
    }
    nCount++;
    if(nCount % 10000 == 0) {
      int prime = 0;
      for(int i = 0; i < count; i++) {
        if(chunk[i][PRIME_ID_USE] != 0) {
          prime++; 
        }
      }
      if(prime == 0) {
        //logging(LOG_DEBUG, "findChunkSkip no live prime End\n");
        return;
      }
      //logging(LOG_DEBUG, "findChunkSkip [%llu] Candidate Count [%d] live[%d]\n", pNextCandidate->data, nCount, prime);
    }
    oldexponent=pNextCandidate->data;
    pNextCandidate=pNextCandidate->next;
  }
  logging(LOG_DEBUG, "findChunkSkip End\n");
}

void findChunkEuler() {
  uint64_t chunk[g_nChunk][PRIME_ID];
  memset(chunk, 0x00, PRIME_ID*g_nChunk*sizeof(uint64_t));
  int count = getNextChunk(chunk);
  struct Node* pNextCandidate = g_lstCandidate;
  while(pNextCandidate) {
    int exponent = pNextCandidate->data;
    int oldexponent = 0;
    int jump = (exponent - oldexponent);
    int unit = 64 - getBitsop(chunk[count-1][PRIME_ID_PRIME])-1;

    for(int i=0; i < jump/unit + 1; i++) {
      for(int j=0; j<count; j++) {
        uint64_t p = chunk[j][PRIME_ID_PRIME];
        // only 2p + 1
        for(int k=2; k<=2; k+=2) {
          // p > 3
          if((p-1)/k <= 3) break;
          if((p-1)/k > g_nMaxCandidate) continue;
          // k < 100 and k is not a multiple of 3.
          if( k%3 ==0 ) continue;

          // if p ≡ 3 (mod 4) and p > 3, then the prime 2p+1 divides the Mersenne number M(p).
          if(k==2 && ((p-1)/k)%4 == 3) {
            setNumberCache((p -1)/k, MERSENNE_NUMBER);
            continue;
          }
        }
      }
    }
    oldexponent = pNextCandidate->data;
    pNextCandidate=pNextCandidate->next;
  }
}

struct Node* getPrimeList() {
  g_nMinPrime = getLastTDM();

  struct Node* prime = (struct Node*)malloc(sizeof(struct Node*));
  prime->data = g_nMinPrime;
  prime->next = NULL;
  struct Node* pNextPrime = prime;
  for(int i = g_nMinPrime + 2; i < g_nMaxPrime; i+=2) {
    int remain = i%8;
    if(remain == 1 || remain == 7) {
      if(isPrime(i)) {
        struct Node * next = (struct Node*)malloc(sizeof(struct Node *));
        next->data = i;
        next->next = NULL;
        pNextPrime->next = (struct Node*)next;
        pNextPrime = next;
      }
    }
  }

  int nCount = 0;
  pNextPrime = prime;
  while(pNextPrime) {
    nCount++;
    pNextPrime=pNextPrime->next;
  }
  logging(LOG_DEBUG, "Prime Count [%d]\n", nCount);

  pNextPrime = prime;
  return pNextPrime;
}

struct Node* getCandidateList() {

  readBits();
  struct Node* candidate = (struct Node*)malloc(sizeof(struct Node *));
  candidate->data = g_nSystemMinCandidate;
  candidate->next = NULL;

  struct Node* pNextCandidate = candidate;
  for(int i = g_nSystemMinCandidate + 2; i < g_nSystemMaxCandidate; i+=2) {
    if(isCandidate(i) && !isKnownMersenne(i)) {
      struct Node * next = (struct Node*)malloc(sizeof(struct Node *));
      next->data = i;
      next->next = NULL;
      pNextCandidate->next = (struct Node*)next;
      pNextCandidate = next;
    }
  }

  int nCandidateCount = 0;
  pNextCandidate = candidate;
  while(pNextCandidate) {
    nCandidateCount++;
    pNextCandidate=pNextCandidate->next;
  }
  logging(LOG_DEBUG, "getCandidateList [%d]\n", nCandidateCount);

  pNextCandidate = candidate;
  return pNextCandidate;
}

int getNextChunk(uint64_t chunk[][PRIME_ID]) {
  lock_enter();
  int count = 0;
  struct Node* ppos = g_pNextPrime;
  for(int i = 0; i < g_nChunk && ppos->next; i++) {
    chunk[i][PRIME_ID_USE] = 1;
    chunk[i][PRIME_ID_VALUE] = 1;
    chunk[i][PRIME_ID_PRIME] = ppos->data;
    ppos=ppos->next;
    count++;
  }
  g_pNextPrime=ppos->next;
  lock_leave();
  return count;
}

void freeList(struct Node* pList) {
  struct Node* pNextPrime = pList;
  while(pNextPrime) {
    struct Node * pos = pNextPrime;
    pNextPrime= pos->next;
    free(pos);
  }
  pList = NULL;
}

void PrintCandidateCount() {
  static int oldcount = MAX_CANDIDATE;
  int candidateCount = getCandidateCount();
  if(oldcount > candidateCount) {
    logging(LOG_DEBUG, "Mersenne Candidate : [%d], pid=[%ld]\n", candidateCount, (long)getpid());
    oldcount = candidateCount;
  }
}

bool isPrime(uint64_t num) {
#ifdef EXPERIMENTAL
    is_prime(num);
#else
  isPrimeNormal(num);
#endif
}

bool isPrimeNormal(uint64_t num) {
  if(num == 2) return true;
  if(num % 2 == 0) return false;
  for(uint64_t i = 3; i*i<=num; i+=2) {
    if(num % i == 0) return false;
  }
  return true;
}

#ifdef EXPERIMENTAL
bool is_prime(uint64_t n) {
    if (n <= 1) {
        return false; // 1 and numbers less than 1 are not prime
    }
    else if (n <= 3) {
        return true; // 2 and 3 are prime
    }
    else if (n % 2 == 0 || n % 3 == 0) {
        return false; // numbers divisible by 2 or 3 are not prime
    }
    uint64_t i = 5;
    while (i * i <= n) { // check divisibility up to square root of n
        if (n % i == 0 || n % (i + 2) == 0) { // check divisibility by numbers of the form 6k +/- 1
            return false;
        }
        i += 6;
    }
    return true;
}
#endif

bool isPrimeFromDB(uint64_t num) {
  
  return true;
}

void test_get_primes() {
  uint64_t ppos = 2;
  uint64_t maxDiff = 1;
  uint64_t primeCount = 1;
  uint64_t arr[MAX_DIFF];
  int nThreads = 10;

  if(g_nUseOMP == USE_OMP) {
#pragma omp parallel for num_threads(10)
    for(int i = 1; i <= nThreads; i++) {
      findPrime(i*100000000);
    }
  } else {
    for(int i = 1; i <= nThreads; i++) {
      findPrime(i*100000000);
    }
  }
}

void primeDifferenceStatistics() {
  uint64_t ppos = 2;
  uint64_t maxDiff = 1;
  uint64_t primeCount = 1;
  uint64_t arr[MAX_DIFF];

  memset(arr, 0x00, MAX_DIFF*sizeof(uint64_t));
  for(uint64_t i = 3; i < 1000000000; i+=2) {
    if(isPrime(i)) {
      int diff = i - ppos;
      arr[diff/2]++;
      if(maxDiff < diff) {
        maxDiff = diff;
      }
      primeCount++;
      if(primeCount % 1000000 == 0) {
        logging(LOG_DEBUG, "[i=%llu] maxDiff : [%llu] primeCount : [%llu]\n", i, maxDiff, primeCount);
        for(int j = 0; j < MAX_DIFF; j++) {
          logging(LOG_DEBUG, "[%llu]", arr[j]);
        }
        logging(LOG_DEBUG, "\n");
      }
      ppos = i;
    }
  }
  logging(LOG_DEBUG, "[i=1000000000] maxDiff : [%llu] primeCount : [%llu]\n", maxDiff, primeCount);
  for(int i = 1; i < MAX_DIFF; i++) {
    logging(LOG_DEBUG, "[%llu]", arr[i]);
  }
}

#ifdef USE_COMPRESSION
void test_compress() {
  char* source = "This is compression sample";
  char* decoded = "";

  HFNode* tree = NULL;
  bitBuffer encoded = {0, NULL};
  HFCode codeTable[MAX_CHAR];

  priorityQueue* pq = Create(3);
  PQNode result;

  PQNode nodes[6] = {
      {34, (void*)"Node1"},
      {12, (void*)"Node2"},
      {87, (void*)"Node3"},
      {45, (void*)"Node4"},
      {35, (void*)"Node5"},
      {66, (void*)"Node6"}
  };

  enqueue(pq, nodes[0]);
  enqueue(pq, nodes[1]);
  enqueue(pq, nodes[2]);
  enqueue(pq, nodes[3]);
  enqueue(pq, nodes[4]);
  enqueue(pq, nodes[5]);

  logging(LOG_DEBUG, "Count of jobs on the queue : %d\n", pq->usedSize);
  while(!isEmpty(pq)) {
    dequeue(pq, &result);
    printNode(&result);
  }

  memset(&codeTable, 0, sizeof(HFCode)*MAX_CHAR);
  encode(&tree, (unsigned char*)source, &encoded, codeTable);

  logging(LOG_DEBUG, "Original Size: %d-bit, encoded Size: %d-bit\n",
    (int)((strlen(source) + 1) * sizeof(char) * 8), encoded.size);

  decoded = (char*)malloc(sizeof(char) * (strlen(source) + 1));
  decode(tree, &encoded, (unsigned char*)decoded);

  logging(LOG_DEBUG, "Original: %s\n", source);
  logging(LOG_DEBUG, "encoded: ");
  printBinary(&encoded);
  logging(LOG_DEBUG, "\ndecoded: %s\n", decoded);

  free(decoded);
  destroyTree(tree);

  return 0;
}
#endif

void readCommand(int argc, char ** argv) {

  memset(g_achOptionPath, 0x00, MAX_SIZE);
  
  if(argv[1] == NULL) {
    g_nMinCandidate = MIN_CANDIDATE;
    g_nMaxCandidate = MAX_CANDIDATE;
    g_nMinPrime    = MIN_PRIME;
    g_nMaxPrime    = MAX_PRIME;
    g_nRunMode = RUN_ALL;
    g_nAlgorithm = ALGO_KJ;
    g_nUseFFT = USE_NOFFT;
    g_nUseOMP = USE_NOOMP;
    g_nUseSocket = USE_DIRECT;
    g_nNoFile = USE_FILE;
    g_nFinder = FINDER_NORMAL;
    g_nTDMProcessCount = 1;
    g_nLLTProcessCount = 1;
    g_nLogFile = LOG_FILE;
    g_nLogNetwork = LOG_NONE;
    g_nLogVerbose = LOG_VERBOSE;
    g_nLogLevel = LOG_INFO;
    g_nDisplayCount = 1024;
    g_nDisplayStart = MIN_CANDIDATE;
    g_nDisplayEnd = MIN_CANDIDATE + g_nDisplayCount;
  } else {
    strcpy(g_achOptionPath, argv[1]);

    g_nMinCandidate = configInt("MinCandidate");
    g_nMaxCandidate = configInt("MaxCandidate");
    g_nMinPrime = configInt("MinPrime");
    g_nMaxPrime = configInt("MaxPrime");
    g_nRunMode = configInt("RunMode");
    g_nAlgorithm = configInt("Algorithm");
    g_nUseFFT = configInt("UseFFT");
    g_nUseOMP = configInt("UseOMP");
    g_nUseSocket = configInt("UseSocket");
    g_nNoFile = configInt("NoFile");
    g_nFinder = configInt("Finder");
    g_nTDMProcessCount = configInt("TDMProcess");
    g_nLLTProcessCount = configInt("LLTProcess");
    g_nLogFile = configInt("LogFile");
    g_nLogNetwork = configInt("LogNetwork");
    g_nLogVerbose = configInt("LogVerbose");
    g_nLogLevel = configInt("LogLevel");
    g_nDisplayCount = configInt("DisplayCount");
    g_nDisplayStart = configInt("DisplayStart");
    g_nDisplayEnd = configInt("DisplayEnd");

    logging(LOG_DEBUG, "COMMAND [%s]\n", argv[1]);
  }

  if(g_nRunMode == RUN_KILL) {
    logging(LOG_DEBUG, "killAll received!\n");
    killAll();
    initConfig();
    exit(0);
  }
 
  if(g_nRunMode == RUN_CONSOLE) {
  } else {
    // if app was crashed, we clear a previous lock.
    initConfig();
  }
  memset(candidate, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nSystemMaxCandidate/64 + 1));
}

void initMT19937() {
  mt19937_initialize(&mt_state, omp_get_thread_num());
}

int configInt(char *pKey) {
  char pTemp[MAX_SIZE];
  memset(pTemp, 0x00, MAX_SIZE);
  readConfig(pKey, pTemp);
  return atoi(pTemp);
}

void readConfig(char * pKey, char * pValue) {
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;
  if((fp =fopen(g_achOptionPath, "r")) != NULL ){
    while(fgets(buff, MAX_SIZE, fp)){
      buff[strlen(buff) -1] = '\0';
      if(strncmp(buff, pKey, strlen(pKey)) == 0) {
        memcpy(pValue, buff + strlen(pKey) + 1, strlen(buff) - strlen(pKey) -1);
        return;
      }
    }
    fclose(fp);
  }
}

void readBits() {
  lock_enter();
#ifdef USE_COMPRESSION
  uncompress(CANDIDATE_FILE_NAME, candidate);
#else
#ifdef USE_ENCRYPTION
  decryptFile(CANDIDATE_FILE_NAME, candidate);
#else
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;

  memset(candidate, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nSystemMaxCandidate/64 + 1));

  if((fp =fopen(CANDIDATE_FILE_NAME, "r")) != NULL ){
    while(!feof(fp)){
      memset(buff, 0x00, MAX_SIZE);
      size = fread(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      memcpy(candidate+pos,buff, size);
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
#endif
#endif
  lock_leave();
}

void writeBits() {
  lock_enter();
#ifdef USE_COMPRESSION
  compress(CANDIDATE_FILE_NAME, candidate);
#else
#ifdef USE_ENCRYPTION
  encryptFile(CANDIDATE_FILE_NAME, candidate);
#else
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;
  if((fp = fopen(CANDIDATE_FILE_NAME, "w")) != NULL ) {
    for(int i = 0; i < (g_nSystemMaxCandidate/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
      memset(buff, 0x00, MAX_SIZE);
      memcpy(buff, candidate+pos, MAX_SIZE);
      size = fwrite(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
#endif
#endif
  lock_leave();
}

void updateBits() {
  memcpy(cache, candidate, sizeof(uint64_t)*(g_nSystemMaxCandidate/64 + 1));
  readBits();
  int shift = g_nMinCandidate%64;
  candidate[g_nMinCandidate/64] = ((candidate[g_nMinCandidate/64] >> shift) << shift) + ((cache[g_nMinCandidate/64] << shift) >> shift);
  for(int i = g_nMinCandidate/64 + 1; i < g_nMaxCandidate/64; i++) {
    candidate[i] = cache[i];
  }
  candidate[g_nMaxCandidate/64 + 1] = ((cache[g_nMaxCandidate/64 + 1] >> shift) << shift) + ((candidate[g_nMaxCandidate/64 + 1] << shift) >> shift);
  writeBits();
}

void readDB() {
  lock_enter();
#ifdef USE_COMPRESSION
  uncompress(DB_FILE_NAME, candidate);
#else
#ifdef USE_ENCRYPTION
  decryptFile(DB_FILE_NAME, candidate);
#else
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;

  memset(candidate, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nSystemMaxCandidate/64 + 1));

  if((fp =fopen(DB_FILE_NAME, "r")) != NULL ){
    while(!feof(fp)){
      memset(buff, 0x00, MAX_SIZE);
      size = fread(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      memcpy(candidate+pos,buff, size);
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
#endif
#endif
  lock_leave();
}

void writeDB() {
  lock_enter();
#ifdef USE_COMPRESSION
  compress(DB_FILE_NAME, candidate);
#else
#ifdef USE_ENCRYPTION
  encryptFile(DB_FILE_NAME, candidate);
#else
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;
  if((fp = fopen(DB_FILE_NAME, "w")) != NULL ) {
    for(int i = 0; i < (g_nSystemMaxCandidate/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
      memset(buff, 0x00, MAX_SIZE);
      memcpy(buff, candidate+pos, MAX_SIZE);
      size = fwrite(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
#endif
#endif
  lock_leave();
}

void updateDB() {
  memcpy(cache, candidate, sizeof(uint64_t)*(g_nSystemMaxCandidate/64 + 1));
  readDB();
  int shift = g_nMinCandidate%64;
  candidate[g_nMinCandidate/64] = ((candidate[g_nMinCandidate/64] >> shift) << shift) + ((cache[g_nMinCandidate/64] << shift) >> shift);
  for(int i = g_nMinCandidate/64 + 1; i < g_nMaxCandidate/64; i++) {
    candidate[i] = cache[i];
  }
  candidate[g_nMaxCandidate/64 + 1] = ((cache[g_nMaxCandidate/64 + 1] >> shift) << shift) + ((candidate[g_nMaxCandidate/64 + 1] << shift) >> shift);
  writeDB();
}

int getCompressSize(uint64_t maxPrime) {
  if(maxPrime < 100000000) {
    return 32;
  } else if(maxPrime >= 100000000 && maxPrime < 1000000000) {
    return 64;
  } else if(maxPrime >= 1000000000 && maxPrime < 100000000000) {
    return 128;
  } else if(maxPrime >= 10000000000 && maxPrime < 1000000000000) {
    return 256;
  }
  return 256;
}

#ifdef USE_COMPRESSION
void getPrime(uint64_t maxPrime, uint64_t * primes) {
  uncompressPrimes(maxPrime, primes);
}
#endif

void findPrime(uint64_t maxPrime) {
  int nPrimeIndex = 0;
  int nNoiseIndex = 0;
  int nCompressSize = getCompressSize(maxPrime);
  uint8_t primes[MAX_PRIMES];
  uint32_t noises[MAX_NOISES];
  uint64_t ppos = maxPrime - 100000000 -1;
  memset(primes, 0x00, sizeof(uint8_t)*MAX_PRIMES);
  memset(noises, 0x00, sizeof(uint8_t)*MAX_NOISES);
  for(uint64_t i = maxPrime - 100000000; i < maxPrime; i++) {
    if(isPrime(i)) {
      int diff = (i - ppos)/2;
      if(diff > nCompressSize) {
        diff -= nCompressSize;
        noises[nNoiseIndex++] = i - maxPrime + 100000000;
      }
      primes[nPrimeIndex++] = diff;

      ppos = i;
      if(nPrimeIndex % 1000 == 0 || i == maxPrime - 1) {
        // TODO : Create Index File
        writePrime(maxPrime, primes);
        writeNoise(maxPrime, noises);
        memset(primes, 0x00, sizeof(uint8_t)*MAX_PRIMES);
        memset(noises, 0x00, sizeof(uint8_t)*MAX_NOISES);
        nPrimeIndex = 0;
      }
    }
  }
#ifdef USE_COMPRESSION
  compressPrimes(maxPrime);
#endif
}

void writePrime(uint64_t maxPrime, uint8_t * primes) {
  FILE *fp = NULL;
  char filename[MAX_PATH];
  mkdir("primes", 0755);
  sprintf(filename, "primes/%llu.primes", maxPrime);
  if((fp = fopen(filename, "a+")) != NULL ) {
    fseek(fp, 0, SEEK_END);
    fwrite(&primes, 1, MAX_PRIMES, fp);
    fclose(fp);
  }
}

void writeNoise(uint64_t maxPrime, uint32_t * noises) {
  FILE *fp = NULL;
  char filename[MAX_PATH];
  mkdir("primes", 0755);
  sprintf(filename, "primes/%llu.noises", maxPrime);
  if((fp = fopen(filename, "a+")) != NULL ) {
    fseek(fp, 0, SEEK_END);
    fwrite(&noises, 1, MAX_NOISES, fp);
    fclose(fp);
  }
}

#ifdef USE_COMPRESSION

HFNode* createNode(symbolInfo data) {
  HFNode* node = (HFNode*)calloc(1, sizeof(HFNode));
  node->data = data;
  return node;
}

void destroyNode(HFNode* node) {
  free(node);
}

void destroyTree(HFNode* tree) {
  if(tree == NULL) return;
  destroyTree(tree->left);
  destroyTree(tree->right);
  destroyNode(tree);
}

void addBit(bitBuffer* buffer, unsigned char bit) {
  unsigned char mask = 0x80;
  if (buffer->size % 8 == 0) {
    buffer->buffer = realloc(buffer->buffer,
                             sizeof(unsigned char)*((buffer->size / 8) + 1));
    buffer->buffer[buffer->size / 8] = 0x00; // null character
  }
  mask >>= buffer->size % 8;

  if (bit == 1) buffer->buffer[buffer->size / 8] |= mask;
  else buffer->buffer[buffer->size / 8] &= ~mask;

  buffer->size++;
}

void encode(HFNode** tree, const unsigned char* src, bitBuffer* encoded, HFCode codeTable[MAX_CHAR]) {
  symbolInfo table[MAX_CHAR];
  unsigned char temp[MAX_BIT];
  int i = 0;

  for(; i < MAX_CHAR; ++i) {
    table[i].symbol = i;
    table[i].frequency = 0;
  }

  i = 0;
  while (src[i] != '\0') table[src[i++]].frequency++;

  buildPrefixTree(tree, table);
  buildCodeTable(*tree, codeTable, temp, 0);

  i = 0;
  while(src[i] != '\0') {
    for(int j = 0; j < codeTable[src[i]].size; ++j)
      addBit(encoded, codeTable[src[i]].code[j]);
    i++;
  }
}

void decode(HFNode* tree, bitBuffer* encoded, unsigned char* decoded) {
  int index = 0;
  HFNode* cur = tree;

  for(int i = 0; i <= encoded->size; ++i) {
    unsigned mask = 0x80;
    if(cur->left == NULL && cur->right == NULL) {
      decoded[index++] = cur->data.symbol;
      cur = tree;
    }
    mask >>= i % 8;
    cur = ((encoded->buffer[i / 8] & mask) != mask) ? cur->left: cur->right;
  }
  decoded[index] = '\0';
}

void buildPrefixTree(HFNode** tree, symbolInfo table[MAX_CHAR]) {
  PQNode result;
  priorityQueue* pq = Create(0);

  for( int i = 0; i < MAX_CHAR; ++i) {
    if(table[i].frequency > 0) {
      PQNode newNode;
      newNode.priority = table[i].frequency;
      newNode.data = createNode(table[i]);
      enqueue(pq, newNode);
    }
  }
  while(pq->usedSize > 1) {
    symbolInfo newData = {0, 0};
    HFNode* bitNode = createNode(newData);
    HFNode* left, *right;
    PQNode leftNode, rightNode, newNode;

    dequeue(pq, &leftNode);
    dequeue(pq, &rightNode);

    left = (HFNode*)leftNode.data;
    right = (HFNode*)rightNode.data;

    bitNode->data.symbol = 0;
    bitNode->data.frequency = left->data.frequency + right->data.frequency;
    bitNode->left = left;
    bitNode->right = right;

    newNode.priority = bitNode->data.frequency;
    newNode.data = bitNode;

    enqueue(pq, newNode);
  }
  dequeue(pq, &result);
  *tree = (HFNode*)result.data;
}

void buildCodeTable(HFNode* tree, HFCode codeTable[MAX_CHAR],
                            unsigned char code[MAX_BIT], int size) {
  if(tree == NULL) return;
  if(tree->left) {
    code[size] = 0;
    buildCodeTable(tree->left, codeTable, code, size + 1);
  }
  if(tree->right) {
    code[size] = 1;
    buildCodeTable(tree->right, codeTable, code, size + 1);
  }
  if(tree->left == NULL && tree->right == NULL) {
    for(int i = 0; i < size; ++i)
      codeTable[tree->data.symbol].code[i] = code[i];
    codeTable[tree->data.symbol].size = size;
  }
}

void printBinary(bitBuffer* buffer) {
  for(int i = 0; i < buffer->size; ++i) {
    unsigned char mask = 0x80 >> (i % 8);
    printf("%d", (buffer->buffer[i / 8] & mask) == mask);
  }
}

priorityQueue* Create(int size) {
  priorityQueue* pq = (priorityQueue*)calloc(1, sizeof(priorityQueue));
  pq->nodes = (PQNode*)calloc(size, sizeof(PQNode));
  pq->capacity = size;
  return pq;
}

void destroy(priorityQueue* pq) {
  free(pq->nodes);
  free(pq);
}

void enqueue(priorityQueue* pq, PQNode newNode) {
  int currentPos = pq->usedSize;
  int parentPos = getParent(currentPos);

  if(currentPos == pq->capacity) {
    if(pq->capacity == 0) pq->capacity = 1;
    pq->capacity <<=1;
    pq->nodes = (PQNode*)realloc(pq->nodes, sizeof(PQNode) * pq->capacity);
  }
  pq->nodes[currentPos] = newNode;

  while(currentPos > 0 && pq->nodes[currentPos].priority < pq->nodes[parentPos].priority) {
    swapNodes(pq, currentPos, parentPos);
    currentPos = parentPos;
    parentPos = getParent(currentPos);
  }
  pq->usedSize++;
}

void dequeue(priorityQueue* pq, PQNode* root) {
  int selectedChild = 0;
  int parentPos = 0, leftPos = 0, rightPos = 0;

  memcpy(root, &pq->nodes[0], sizeof(PQNode));
  memset(&pq->nodes[0], 0, sizeof(PQNode));
  swapNodes(pq, 0, --(pq->usedSize));

  leftPos = getLeftChild(parentPos);
  rightPos = getRightChild(parentPos);

  while(1) {
    if (leftPos >= pq->usedSize) break;
    if (rightPos >= pq->usedSize) selectedChild = leftPos;
    else selectedChild = (pq->nodes[leftPos].priority > pq->nodes[rightPos].priority) ? rightPos : leftPos;

    if (pq->nodes[selectedChild].priority >= pq->nodes[parentPos].priority) break;
    swapNodes(pq, parentPos, selectedChild);
    parentPos = selectedChild;
    leftPos = getLeftChild(parentPos);
    rightPos = getRightChild(parentPos);
  }
  if (pq->usedSize < pq->capacity / 2) {
    pq->capacity >>= 1;
    pq->nodes = (PQNode*)realloc(pq->nodes, sizeof(PQNode) * pq->capacity);
  }
}

int getParent(int pos) {
  return (int)((pos - 1) / 2);
}

int getLeftChild(int pos) {
  return (2 * pos + 1);
}

int getRightChild(int pos) {
  return (2 * pos + 2);
}

void swapNodes(priorityQueue* pq, int pos1, int pos2) {
  PQNode* temp = (PQNode*)malloc(sizeof(PQNode));
  memcpy(temp,             &pq->nodes[pos1], sizeof(PQNode));
  memcpy(&pq->nodes[pos1], &pq->nodes[pos2], sizeof(PQNode));
  memcpy(&pq->nodes[pos2], temp, sizeof(PQNode));
  free(temp);
}

int isEmpty(priorityQueue* pq) {
  return (pq->usedSize == 0);
}

void printNode (PQNode* node) {
  printf("Task name: %s, Priority: %d\n", (char*)node->data, node->priority);
}

void compress(char * pchFileName, uint64_t * data) {
  uint64_t * compress_data = NULL;
  
#ifdef USE_ENCRYPTION
  encryptFile(pchFileName, compress_data);
#else
   
#endif
}

void uncompress(char * pchFileName, uint64_t * data) {
  uint64_t * compress_data = NULL;
#ifdef USE_ENCRYPTION
  decryptFile(pchFileName, compress_data);
#else

#endif
}

void compressPrimes(uint64_t maxPrime) {

}

void uncompressPrimes(uint64_t maxPrime, uint64_t * primes){

}

#endif

#ifdef USE_ENCRYPTION
// block encrypt
static void key_exp(uint8_t* rkey, const uint8_t* key) {
  unsigned i, j, k;
  uint8_t tempa[4];
  for (i = 0; i < Nk; ++i) {
    rkey[(i * 4) + 0] = key[(i * 4) + 0];
    rkey[(i * 4) + 1] = key[(i * 4) + 1];
    rkey[(i * 4) + 2] = key[(i * 4) + 2];
    rkey[(i * 4) + 3] = key[(i * 4) + 3];
  }
  for (i = Nk; i < Nb * (Nr + 1); ++i) {
    {
      k = (i - 1) * 4;
      tempa[0]=rkey[k + 0];
      tempa[1]=rkey[k + 1];
      tempa[2]=rkey[k + 2];
      tempa[3]=rkey[k + 3];

    }
    if (i % Nk == 0) {
      {
        const uint8_t u8tmp = tempa[0];
        tempa[0] = tempa[1];
        tempa[1] = tempa[2];
        tempa[2] = tempa[3];
        tempa[3] = u8tmp;
      }
      {
        tempa[0] = get_sb(tempa[0]);
        tempa[1] = get_sb(tempa[1]);
        tempa[2] = get_sb(tempa[2]);
        tempa[3] = get_sb(tempa[3]);
      }

      tempa[0] = tempa[0] ^ rcon[i/Nk];
    }
    if (i % Nk == 4) {
      {
        tempa[0] = get_sb(tempa[0]);
        tempa[1] = get_sb(tempa[1]);
        tempa[2] = get_sb(tempa[2]);
        tempa[3] = get_sb(tempa[3]);
      }
    }
    j = i * 4; k=(i - Nk) * 4;
    rkey[j + 0] = rkey[k + 0] ^ tempa[0];
    rkey[j + 1] = rkey[k + 1] ^ tempa[1];
    rkey[j + 2] = rkey[k + 2] ^ tempa[2];
    rkey[j + 3] = rkey[k + 3] ^ tempa[3];
  }
}
static void add_rkey(uint8_t round, state_t* state, const uint8_t* rkey) {
  uint8_t i,j;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      (*state)[i][j] ^= rkey[(round * Nb * 4) + (i * Nb) + j];
    }
  }
}
static void subbytes(state_t* state) {
  uint8_t i, j;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      (*state)[j][i] = get_sb((*state)[j][i]);
    }
  }
}
static void shiftrows(state_t* state) {
  uint8_t temp;
  temp           = (*state)[0][1];
  (*state)[0][1] = (*state)[1][1];
  (*state)[1][1] = (*state)[2][1];
  (*state)[2][1] = (*state)[3][1];
  (*state)[3][1] = temp;
  temp           = (*state)[0][2];
  (*state)[0][2] = (*state)[2][2];
  (*state)[2][2] = temp;
  temp           = (*state)[1][2];
  (*state)[1][2] = (*state)[3][2];
  (*state)[3][2] = temp;
  temp           = (*state)[0][3];
  (*state)[0][3] = (*state)[3][3];
  (*state)[3][3] = (*state)[2][3];
  (*state)[2][3] = (*state)[1][3];
  (*state)[1][3] = temp;
}
static uint8_t xtime(uint8_t x) {
  return ((x<<1) ^ (((x>>7) & 1) * 0x1b));
}
static void mixcolumns(state_t* state) {
  uint8_t i;
  uint8_t Tmp, Tm, t;
  for (i = 0; i < 4; ++i) {
    t   = (*state)[i][0];
    Tmp = (*state)[i][0] ^ (*state)[i][1] ^ (*state)[i][2] ^ (*state)[i][3] ;
    Tm  = (*state)[i][0] ^ (*state)[i][1] ; Tm = xtime(Tm);  (*state)[i][0] ^= Tm ^ Tmp ;
    Tm  = (*state)[i][1] ^ (*state)[i][2] ; Tm = xtime(Tm);  (*state)[i][1] ^= Tm ^ Tmp ;
    Tm  = (*state)[i][2] ^ (*state)[i][3] ; Tm = xtime(Tm);  (*state)[i][2] ^= Tm ^ Tmp ;
    Tm  = (*state)[i][3] ^ t ;              Tm = xtime(Tm);  (*state)[i][3] ^= Tm ^ Tmp ;
  }
}
static void mixcolumns_inv(state_t* state) {
  int i;
  uint8_t a, b, c, d;
  for (i = 0; i < 4; ++i) {
    a = (*state)[i][0];
    b = (*state)[i][1];
    c = (*state)[i][2];
    d = (*state)[i][3];

    (*state)[i][0] = multiply(a, 0x0e) ^ multiply(b, 0x0b) ^ multiply(c, 0x0d) ^ multiply(d, 0x09);
    (*state)[i][1] = multiply(a, 0x09) ^ multiply(b, 0x0e) ^ multiply(c, 0x0b) ^ multiply(d, 0x0d);
    (*state)[i][2] = multiply(a, 0x0d) ^ multiply(b, 0x09) ^ multiply(c, 0x0e) ^ multiply(d, 0x0b);
    (*state)[i][3] = multiply(a, 0x0b) ^ multiply(b, 0x0d) ^ multiply(c, 0x09) ^ multiply(d, 0x0e);
  }
}
static void subbytes_inv(state_t* state) {
  uint8_t i, j;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      (*state)[j][i] = get_sb_inv((*state)[j][i]);
    }
  }
}
static void shiftrows_inv(state_t* state) {
  uint8_t temp;
  temp = (*state)[3][1];
  (*state)[3][1] = (*state)[2][1];
  (*state)[2][1] = (*state)[1][1];
  (*state)[1][1] = (*state)[0][1];
  (*state)[0][1] = temp;
  temp = (*state)[0][2];
  (*state)[0][2] = (*state)[2][2];
  (*state)[2][2] = temp;
  temp = (*state)[1][2];
  (*state)[1][2] = (*state)[3][2];
  (*state)[3][2] = temp;
  temp = (*state)[0][3];
  (*state)[0][3] = (*state)[1][3];
  (*state)[1][3] = (*state)[2][3];
  (*state)[2][3] = (*state)[3][3];
  (*state)[3][3] = temp;
}
static void cipher(state_t* state, const uint8_t* rkey) {
  uint8_t round = 0;
  add_rkey(0, state, rkey);
  for (round = 1; ; ++round) {
    subbytes(state);
    shiftrows(state);
    if (round == Nr) {
      break;
    }
    mixcolumns(state);
    add_rkey(round, state, rkey);
  }
  add_rkey(Nr, state, rkey);
}
static void cipher_inv(state_t* state, const uint8_t* rkey) {
  uint8_t round = 0;
  add_rkey(Nr, state, rkey);
  for (round = (Nr - 1); ; --round) {
    shiftrows_inv(state);
    subbytes_inv(state);
    add_rkey(round, state, rkey);
    if (round == 0) {
      break;
    }
    mixcolumns_inv(state);
  }
}
static void xoriv(const uint8_t* iv) {
  uint8_t i;
  for (i = 0; i < blockSize; ++i) {
    buf[i] ^= iv[i];
  }
}
void encryptFile(char * pchFileName, uint64_t * data) {
  size_t i;
  uint8_t *liv = iv;
  for (i = 0; i < len; i += blockSize) {
    xoriv(liv);
    cipher((state_t*)buf, rkey);
    liv = (uint8_t*)buf;
    buf += blockSize;
  }
  memcpy(iv, liv, blockSize);
}
void decryptFile(char * pchFileName, uint64_t * data) {
  size_t i;
  uint8_t liv[blockSize];
  for (i = 0; i < len; i += blockSize) {
    memcpy(liv, buf, blockSize);
    cipher_inv((state_t*)buf, rkey);
    xoriv(iv);
    memcpy(iv, liv, blockSize);
    buf += blockSize;
  }
}

// asymmetric encrypt
uint64_t findD(uint64_t e, uint64_t phi) {
  uint64_t eprev, dprev, d = 1, etemp, dtemp;
  eprev = phi, dprev = phi;
  while (e != 1) {
    etemp = e;
    dtemp = d;
    e = eprev - eprev / etemp * e;
    d = dprev - eprev / etemp * d;
    eprev = etemp;
    dprev = dtemp;
    while (d < 0)
       d += phi;
  }
  return d;
}
uint64_t gcd(uint64_t num1, uint64_t num2) {
  uint64_t i, temp;
  if (num1 > num2) {
    temp = num1;
    num1 = num2;
    num2 = temp;
  }
  for (i = num1; i > 0; i--) {
    if (num1 % i == 0 && num2 % i == 0)
      return i;
  }
  return 0;
}
uint64_t getprime() {
  uint64_t n;
  do {
    srand(time(NULL));
    n = rand() % MAX_NUMBER + 5;
  } while(!isPrime(n));
  return n;
}
void setprimes(uint64_t e, uint64_t *p, uint64_t *q, uint64_t *n, uint64_t *phi) {
  do {
    *p = getprime();
    do
      *q = getprime();
  while(*p == *q);
    *n = *p * *q;
    *phi = *n - *p - *q + 1;
   } while (gcd(e,*phi) != 1);
}
uint64_t modpow(uint64_t base, uint64_t power, uint64_t mod) {
  uint64_t result = 1;
  for (int i = 0; i < power; i++) {
    result = (result * base) % mod;
  }
  return result;
}
// find inverse of number a in % mod
// inverse of a = i
uint64_t inverse(int a, int mod) {
  uint64_t aprev, iprev, i = 1, atemp, itemp;
  aprev = mod, iprev = mod;
  while (a != 1) {
    atemp = a;
    itemp = i;
    a = aprev - aprev / atemp * a;
    i = iprev - aprev / atemp * i;
    aprev = atemp;
    iprev = itemp;
    while (i < 0)
      i += mod;
  }
  return i;
}

void test_sign () {
    uint64_t e = E_VALUE;
    uint64_t phi = 0;

    uint64_t d, n, p, q, h, m, qInv, m1m2;
    uint64_t c, dP, dQ, m1, m2;

    logging(LOG_DEBUG, "p and q must be prime numbers, e must be coprime to (p - 1)*(q - 1)\n");
    setprimes(e, &p, &q, &n, &phi);

    logging(LOG_DEBUG, "p: %llu q: %llu n: %llu phi: %llu\n", p, q, n, phi);

    d = findD(e,phi);
    logging(LOG_DEBUG, "find d: %llu\n", d);

    logging(LOG_DEBUG, "Public Key:  (n,e) = (%llu, %llu)\n", n, e);
    logging(LOG_DEBUG, "Private Key: (n,d) = (%llu, %llu)\n", n, d);
    logging(LOG_DEBUG, ("p = %llu, q = %llu, phi = %llu\n", p, q, phi);
    logging(LOG_DEBUG, ("Key generation success!\n\n");

    char message[MAX_SIZE];
    memset(message, 0x00, MAX_SIZE);
    strcpy(message, "This is a plain Text!");
    logging(LOG_DEBUG, "original message [%s]\n", message);
 
    uint64_t cipher[MAX_SIZE];
    memset(cipher, 0x00, sizeof(uint64_t)*MAX_SIZE);
    for (int i = 0; i < strlen(message) ; i++) {
      cipher[i] = modpow(message[i],e,n);
    }
    
    char decrypted[MAX_SIZE];
    memset(decrypted, 0x00, MAX_SIZE);
   
    for (int i = 0; i < strlen(message) ; i++) {
      uint64_t c = cipher[i];
      dP = d % (p - 1);
      dQ = d % (q - 1);
      qInv = inverse(q,p);
      m1 = modpow(c,dP,p);
      m2 = modpow(c,dQ,q);
      m1m2 = m1 - m2;
      if (m1m2 < 0)
        m1m2 += p;
      h = (qInv * m1m2) % p;
      m = m2 + h * q;
        decrypted[i] = m;
    }
    logging(LOG_DEBUG, "decrypted [%s]\n", decrypted);
    return 0;
}
#endif

bool setNumberCache(uint64_t nIndex, int newBit) {
  uint64_t value = candidate[nIndex/64];
  uint64_t residue = nIndex % 64;
  bool flag = false;
  int bit = (value >> residue) & 0x01;
  if(newBit == 0 && bit == 1) {
    candidate[nIndex/64]-=(uint64_t)pow(2, residue);
    flag = true;
  }
  if(newBit == 1 && bit == 0) {
    candidate[nIndex/64]+=(uint64_t)pow(2, residue);
    flag = true;
  }
  return flag;
}

void setNumber(uint64_t nIndex, int bit) {
  if(setNumberCache(nIndex, bit)) {
    writeBits();
  }
}

bool isCandidate(uint64_t nIndex) {
  uint64_t value = candidate[nIndex/64];
  uint64_t residue = nIndex % 64;
  return (value >> residue) & 0x01;
}

bool isKnownMersenne(int number) {
  for(int i = 0; i < KNOWN_MERSENNE; i++) {
    if(knownMersenne[i] == number) {
      return true;
    } 
  }
  return false;
}

int getCandidateCount() {
  int nCount = 0;
  for(int i = g_nSystemMinCandidate; i < g_nSystemMaxCandidate; i+=2) {
    if(isCandidate(i) && !isKnownMersenne(i)) {
      nCount++; 
    }
  }
  return nCount;
}

void lock_enter() {
  // You must not enter another lock between lock and unlock.
  while(access(LOCK_FILE_NAME, F_OK) == 0);
  FILE *fp = NULL;
  if((fp = fopen(LOCK_FILE_NAME, "w")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    sprintf(buff, "%ld", (long)getpid());
    fwrite(&buff, 1, MAX_SIZE, fp);
    fclose(fp);
  }
}

void lock_leave() {
  unlink(LOCK_FILE_NAME);
}

void setStatus(char * pchStatus) {
  FILE *fp = NULL;
  if((fp = fopen(STATUS_FILE_NAME, "w")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    sprintf(buff, "%s", pchStatus);
    fwrite(&buff, 1, MAX_SIZE, fp);
    fclose(fp);
  }
}

void getStatus(char * pchStatus) {
  FILE *fp = NULL;
  if((fp = fopen(STATUS_FILE_NAME, "r")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    fread(&buff, 1, MAX_SIZE, fp);
    strcpy(pchStatus, buff);
    fclose(fp);
  } else {
    setStatus(STATUS_EXIT);
  }
}

bool isSingleProcess() {
  return (g_nTDMProcessCount == 1);
}

void updateProcess(char * pchClass, char * pchPID, char * pchStatus) {
  lock_enter();
  struct Row processStatus[MAX_PROCESS];

  memset(g_pchClass, 0x00, MAX_SIZE);
  memset(g_pchPID, 0x00, MAX_SIZE);
  strcpy(g_pchClass, pchClass);
  strcpy(g_pchPID, pchPID);

  memset(processStatus, 0x00, MAX_PROCESS*sizeof(struct Row));
  int nIndex = 0;
  FILE *fp = NULL;
  if((fp = fopen(PROCESS_FILE_NAME, "r")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    while(fgets(buff, MAX_SIZE, fp)){
      buff[strlen(buff) -1] = ',';
      char *ptr = strtok(buff, ",");
      strcpy(processStatus[nIndex].class, ptr);
      ptr = strtok(NULL, ",");
      strcpy(processStatus[nIndex].pid, ptr);
      ptr = strtok(NULL, ",");
      strcpy(processStatus[nIndex].status, ptr);

      nIndex++;
    }
    fclose(fp);
  }
  // If existing process information is available, update the status.
  bool flag = false;
  for(int i=0; i < nIndex; i++) {
    if(strcmp(pchClass, processStatus[i].class) == 0 && strcmp(pchPID, processStatus[i].pid) == 0) {
      strcpy(processStatus[i].status, pchStatus);
      flag = true;
      break;
    } 
  }
  // If it doesn't find it, it adds a new process. 
  if(!flag) {
    strcpy(processStatus[nIndex].class, pchClass); 
    strcpy(processStatus[nIndex].pid, pchPID); 
    strcpy(processStatus[nIndex].status, pchStatus);
  }
  if((fp = fopen(PROCESS_FILE_NAME, "w")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    for(int i=0; i <= nIndex; i++) {
      if(strlen(processStatus[i].pid) > 0) {
        fprintf(fp, "%s,%s,%s\n", processStatus[i].class, processStatus[i].pid, processStatus[i].status);
      }
    }
    fclose(fp);
  }
  lock_leave();
}

void waitTDMReady() {
  do {
#ifdef _WIN32
    Sleep(1000);
#else
    sleep(1);
#endif
  } while(!isTDMReady());
}

bool isTDMReady() {
  lock_enter();
  struct Row processStatus[MAX_PROCESS];
          
  memset(processStatus, 0x00, MAX_PROCESS*sizeof(struct Row));
  int nIndex = 0;
  FILE *fp = NULL;
  if((fp = fopen(PROCESS_FILE_NAME, "r")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    while(fgets(buff, MAX_SIZE, fp)){
      buff[strlen(buff) -1] = ','; 
      char *ptr = strtok(buff, ",");
      strcpy(processStatus[nIndex].class, ptr);
      ptr = strtok(NULL, ","); 
      strcpy(processStatus[nIndex].pid, ptr);
      ptr = strtok(NULL, ",");
      strcpy(processStatus[nIndex].status, ptr);
      nIndex++;
    } 
    fclose(fp);
  }

  int count = 0;
  for(int i=0; i < nIndex; i++) {
    if(strcmp(TDM_PROCESS_NAME, processStatus[i].class) == 0 && strcmp(STATUS_READY, processStatus[i].status) != 0) {
      lock_leave();
      return false;
    }
    if(strcmp(TDM_PROCESS_NAME, processStatus[i].class) == 0 && strcmp(STATUS_READY, processStatus[i].status) == 0) {
      count++;
    } 
  }
  if(g_nTDMProcessCount == count) {
    lock_leave();
    return true;
  }
  lock_leave();
  return false;
}

void killAll() {
  lock_enter();
  struct Row processStatus[MAX_PROCESS];

  memset(processStatus, 0x00, MAX_PROCESS*sizeof(struct Row));
  int nIndex = 0;
  FILE *fp = NULL;
  if((fp = fopen(PROCESS_FILE_NAME, "r")) != NULL ) {
    char buff[MAX_SIZE];
    memset(buff, 0x00, MAX_SIZE);
    while(fgets(buff, MAX_SIZE, fp)){
      buff[strlen(buff) -1] = ',';
      char *ptr = strtok(buff, ",");
      strcpy(processStatus[nIndex].class, ptr);
      ptr = strtok(NULL, ",");
      strcpy(processStatus[nIndex].pid, ptr);
      ptr = strtok(NULL, ",");
      strcpy(processStatus[nIndex].status, ptr);
      nIndex++;
    }
    fclose(fp);
  }
  lock_leave();

  for(int i=0; i < nIndex; i++) {
    kill(atoi(processStatus[i].pid), SIGKILL);
    logging(LOG_INFO, "kill [%s,%s,%s]\n", processStatus[i].class, processStatus[i].pid, processStatus[i].status);
  }
  // TODO : get system process id for kill

  char line[MAX_SIZE];
  pid_t killpid;
  memset(line, 0x00, MAX_SIZE);
  fp = popen("ps -ef | grep mersenne | grep -v grep | grep -v kill","r");
  while (fgets(line,MAX_SIZE,fp)) {
    char temp[MAX_SIZE];
    char temp2[MAX_SIZE];
    char pid[MAX_SIZE];
    memset(pid, 0x00, MAX_SIZE);
    sscanf(line, "%s %s %s", temp,pid,temp2);  
    logging(LOG_INFO, "kill pid=%s\n",pid);
    kill(atoi(pid),SIGKILL);
  }
  pclose(fp);
}

void setLastTDM(long p) {
  if(isSingleProcess()) {
    FILE *fp = NULL;
    if((fp = fopen(LAST_FILE_NAME, "w")) != NULL ) {
      char buff[MAX_SIZE];
      memset(buff, 0x00, MAX_SIZE);
      sprintf(buff, "%ld", p);
      fwrite(&buff, 1, MAX_SIZE, fp);
      fclose(fp);
    }
  }
}

long getLastTDM() {
  long lLastPrime = 0;
  if(isSingleProcess()) {
    FILE *fp = NULL;
    if((fp = fopen(LAST_FILE_NAME, "r")) != NULL ) {
      char buff[MAX_SIZE];
      memset(buff, 0x00, MAX_SIZE);
      fread(&buff, 1, MAX_SIZE, fp);
      lLastPrime = atol(buff);
      fclose(fp);
    } else {
      lLastPrime = g_nMinPrime;
    }
  } else {
    lLastPrime = g_nMinPrime;
  }
  return lLastPrime;
}

void writeLog(char * buf) {
  FILE *fp = NULL;
  time_t timer;
  char date[MAX_DATE];
  struct tm* tm_info;

  timer = time(NULL);
  tm_info = localtime(&timer);

  strftime(date, MAX_DATE, "%Y%m%d%H", tm_info);
  char filename[MAX_PATH];
  memset(filename, 0x00, MAX_PATH);
  sprintf(filename, "./%s.%s.%ld.log", LOG_FILE_NAME, date, (long)getpid());
  if((fp = fopen(filename, "a+")) != NULL ) {
    char buff[MAX_BUFF_SIZE];
    memset(buff, 0x00, MAX_BUFF_SIZE);
    fseek(fp, 0, SEEK_END);
    fprintf(fp, "%s", buf);
    fclose(fp);
  }
}

void sendLog(char * buf) {
  int sock;
  struct sockaddr_in serv_addr;
  int len;

  sock = socket(PF_INET, SOCK_STREAM, 0);
  if(sock == -1) {
    error_handling("socket() error");
  } else {
    memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr(g_achAddress);
    serv_addr.sin_port = htons(g_nPort);

    if(connect(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) == -1) {
      error_handling("connect() error");
    } else {
      char achBuff[MAX_BUFF_SIZE];
      memset(achBuff, 0x00, MAX_BUFF_SIZE);
      sprintf(achBuff, "{ \"log\" : \"%s\" }", buf);
      write(sock, achBuff, strlen(buf));

      len = read(sock, achBuff, MAX_BUFF_SIZE-1);
      close(sock);
    }
  }
}

void initConfig() {
  unlink(LOCK_FILE_NAME);
  unlink(PROCESS_FILE_NAME);
  unlink(STATUS_FILE_NAME);
  unlink(LAST_FILE_NAME);
}

void test_riemann_zeta() {
  double a = 2;
  double ai = 0.0;
  double r = 0.0;
  double ri = 0.0;

  riemann_zeta(a, ai, &r, &ri, 1000000);
  logging(LOG_DEBUG, "riemann_zeta(%.10f + %.10fi) = %.10f + %.10fi\n", a, ai, r, ri);
  logging(LOG_DEBUG, "test_riemann_zeta\n");
}

// Riemman Zeta Function
// Riemann Zeta function can be defined in the complex plane.
//unsigned int maxn=1000;
//double s=0.50;
//double q=14.13472514173470;
//double q=21.02203963877156;
//double q=25.01085758014569;
//double q=30.42487612585951;
//double q=32.93506158773919;
//double q=37.58617815882568;
//double q=40.91871901214750;
//double q=43.32707328091500;
//double q=48.00515088116716;
//double q=49.77383247767230;

// ζ(s) = ζ(p + qi) = 0 means p + qi is non-trivial solution.
void riemann_zeta(double p,double q,double * r,double * ri,int maxNumber) {
  logging(LOG_DEBUG, "initialize riemann_zeta(%.10f + %.10fi) = %.10f + %.10fi\n", p, q, *r, *ri);
  for(int n = 1; n < maxNumber; n++) {
    *r += cos(q*log(n*1.0))/pow(n, p);
    *ri += -sin(q*log(n*1.0))/pow(n, p);
    logging(LOG_DEBUG, "sigma(k=1,k=%d) riemann_zeta(%.10f + %.10fi) = %.10f + %.10fi\n", n + 1, p, q, *r, *ri);
  }
  logging(LOG_DEBUG, "result riemann_zeta(%.10f + %.10fi) = %.10f + %.10fi\n", p, q, *r, *ri);
}

// The Euler product attached to the Riemann zeta function ζ(s), also using the sum of the geometric series.
// If infinite series of i is zero, p + qi is non-trivial solution of riemman_zeta function and i is prime.
int euler_product(double p, double q, int maxNumber) {
  double s = 1.0;
  double si = 0.0;
  for(int i = 3; i < g_nMaxPrime; i += 2) {
    for(int n = 1; n < maxNumber; n++) {
      s += cos(q*log(n*1.0))/pow(n, p);
      si += -sin(q*log(n*1.0))/pow(n, p);
      logging(LOG_DEBUG, "sigma(k=1,k=%d) euler_product(%.2f + %.2fi) = %.2f + %.2fi\n", n + 1, p, q, s, si);
    }
    if(s == 0 && si == 0) {
      return i;
    }
  }
  logging(LOG_DEBUG, "result[%d] riemann_zeta(%.2f + %.2fi) = %.2f + %.2fi\n", p, q, s, si);
}

#ifdef EXPERIMENTAL
void test_zeta() {

    // Define the complex number s
    complex double s = 2.0 + 0.0 * I;

    // Calculate the Riemann zeta function at s
    complex double zeta_s = riemann_zeta_omp(s);

    // Print the result
    logging(LOG_DEBUG, "Zeta(%lf + %lfi) = %lf + %lfi\n", creal(s), cimag(s), creal(zeta_s), cimag(zeta_s));

}

// Calculate Riemann zeta function for a given complex number
complex double riemann_zeta_omp(complex double s) {
    int N = 100; // Number of terms in the series approximation

    complex double result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int n = 1; n <= N; n++) {
        result += cpow(n, -s);
    }
    return result;
}

void test_geometric_series(){
    // Define the first term and common ratio of the geometric series
    complex double first_term = 1.0;
    complex double common_ratio = 0.5;

    // Define the number of terms in the geometric series
    int num_terms = 100;

    // Calculate the sum of the geometric series
    complex double sum = geometric_series(first_term, common_ratio, num_terms);

    // Print the result
    logging(LOG_DEBUG, "Sum of the geometric series: %lf + %lfi\n", creal(sum), cimag(sum));
}

// Compute the sum of a geometric series with a given first term, common ratio, and number of terms
complex double geometric_series(complex double first_term, complex double common_ratio, int num_terms) {
    complex double sum = 0.0;
    complex double term = first_term;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_terms; i++) {
        sum += term;
        term *= common_ratio;
    }
    return sum;
}

// Perform FFT using Cooley-Tukey algorithm
void fft(complex double *x, int n) {
    if (n <= 1) return;
    complex double *even = (complex double *)malloc(n / 2 * sizeof(complex double));
    complex double *odd = (complex double *)malloc(n / 2 * sizeof(complex double));

    for (int i = 0; i < n / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    fft(even, n / 2);
    fft(odd, n / 2);

    for (int i = 0; i < n / 2; i++) {
        complex double t = cexp(-2.0 * I * M_PI * i / n) * odd[i];
        x[i] = even[i] + t;
        x[i + n / 2] = even[i] - t;
    }

    free(even);
    free(odd);
}

// Perform iFFT using Cooley-Tukey algorithm
void ifft(complex double *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = conj(x[i]);
    }

    fft(x, n);

    for (int i = 0; i < n; i++) {
        x[i] = conj(x[i]) / n;
    }
}

// Multiply two BigIntegers using FFT
void multiply_bigint_fft(BigInt *a, BigInt *b, BigInt *result) {
    int n = 2 * max(a->length, b->length);
    complex double *x = (complex double *)malloc(n * sizeof(complex double));
    complex double *y = (complex double *)malloc(n * sizeof(complex double));
    complex double *z = (complex double *)malloc(n * sizeof(complex double));

    // Initialize arrays x and y with values from BigIntegers a and b
    for (int i = 0; i < a->length; i++) {
        x[i] = a->digits[i];
    }
    for (int i = a->length; i < n; i++) {
        x[i] = 0;
    }
    for (int i = 0; i < b->length; i++) {
        y[i] = b->digits[i];
    }
    for (int i = b->length; i < n; i++) {
        y[i] = 0;
    }

    // Perform FFT on arrays x and y
    fft(x, n);
    fft(y, n);

    // Element-wise multiplication of the FFT results
    for (int i = 0; i < n; i++) {
        z[i] = x[i] * y[i];
    }

    // Perform iFFT on array z
    ifft(z, n);

    // Extract the result from array z and update the BigInt result
    int carry = 0;
    for (int i = 0; i < n; i++) {
        int digit = round(creal(z[i])) + carry;
        result->digits[i] = digit % 10;
        carry = digit / 10;
    }

    // Find the length of the result
    result->length = n;
    while (result->digits[result->length - 1] == 0 && result->length > 1) {
        result->length--;
    }

    free(x);
    free(y);
    free(z);
}

// Perform FFT using Cooley-Tukey algorithm
void fft_omp(complex double *x, int n) {
    if (n <= 1) return;
    complex double *even = (complex double *)malloc(n / 2 * sizeof(complex double));
    complex double *odd = (complex double *)malloc(n / 2 * sizeof(complex double));

    #pragma omp parallel for
    for (int i = 0; i < n / 2; i++) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    fft_omp(even, n / 2);
    fft_omp(odd, n / 2);

    #pragma omp parallel for
    for (int i = 0; i < n / 2; i++) {
        complex double t = cexp(-2.0 * I * PI * i / n) * odd[i];
        x[i] = even[i] + t;
        x[i + n / 2] = even[i] - t;
    }

    free(even);
    free(odd);
}

// Perform iFFT using Cooley-Tukey algorithm
void ifft_omp(complex double *x, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        x[i] = conj(x[i]);
    }

    fft_omp(x, n);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        x[i] = conj(x[i]) / n;
    }
}

// Multiply two BigIntegers using FFT
void multiply_bigint_fft_omp(BigInt *a, BigInt *b, BigInt *result) {
    int n = 2 * max(a->length, b->length);
    complex double *x = (complex double *)malloc(n * sizeof(complex double));
    complex double *y = (complex double *)malloc(n * sizeof(complex double));
    complex double *z = (complex double *)malloc(n * sizeof(complex double));

    // Initialize arrays x and y with values from BigIntegers a and b
    for (int i = 0; i < a->length; i++) {
        x[i] = a->digits[i];
    }
    for (int i = a->length; i < n; i++) {
        x[i] = 0;
    }
    for (int i = 0; i < b->length; i++) {
        y[i] = b->digits[i];
    }
    for (int i = b->length; i < n; i++) {
        y[i] = 0;
    }

    // Perform FFT on arrays x and y
    fft_omp(x, n);
    fft_omp(y, n);

    // Element-wise multiplication of the FFT results
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        z[i] = x[i] * y[i];
    }

    // Perform iFFT on array z
    ifft_omp(z, n);

    // Extract the result from array z and update the BigInt result
    int carry = 0;
    for (int i = 0; i < n; i++) {
        int digit = round(creal(z[i])) + carry;
        result->digits[i] = digit % 10;
        carry = digit / 10;
    }

    // Find the length of the result
    result->length = n;
    while (result->digits[result->length - 1] == 0 && result->length > 1) {
        result->length--;
    }

    free(x);
    free(y);
    free(z);
}

// Square a BigInt using FFT
void square_bigint_fft(BigInt *a, BigInt *result) {
    int n = 2 * a->length;
    complex double *x = (complex double *)malloc(n * sizeof(complex double));

    // Initialize array x with values from BigInt a
    for (int i = 0; i < a->length; i++) {
        x[i] = a->digits[i];
    }
    for (int i = a->length; i < n; i++) {
        x[i] = 0;
    }

    // Perform FFT on array x
    fft(x, n);

    // Element-wise squaring of the FFT results
    for (int i = 0; i < n; i++) {
        x[i] *= x[i];
    }

    // Perform iFFT on array x
    ifft(x, n);

    // Extract the result from array x and update the BigInt result
    int carry = 0;
    for (int i = 0; i < n; i++) {
        int digit = round(creal(x[i])) + carry;
        result->digits[i] = digit % 10;
        carry = digit / 10;
    }

    // Find the length of the result
    result->length = n;
    while (result->digits[result->length - 1] == 0 && result->length > 1) {
        result->length--;
    }

    free(x);
}

// Square a BigInt using FFT with OpenMP parallelization
void square_bigint_fft_omp(BigInt *a, BigInt *result) {
    int n = 2 * a->length;
    complex double *x = (complex double *)malloc(n * sizeof(complex double));

    // Initialize array x with values from BigInt a
    for (int i = 0; i < a->length; i++) {
        x[i] = a->digits[i];
    }
    for (int i = a->length; i < n; i++) {
        x[i] = 0;
    }

    // Perform FFT on array x
    fft_omp(x, n);

    // Element-wise squaring of the FFT results
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        x[i] *= x[i];
    }

    // Perform iFFT on array x
    ifft_omp(x, n);

    // Extract the result from array x and update the BigInt result
    int carry = 0;
    for (int i = 0; i < n; i++) {
        int digit = round(creal(x[i])) + carry;
        result->digits[i] = digit % 10;
        carry = digit / 10;
    }

    // Find the length of the result
    result->length = n;
    while (result->digits[result->length - 1] == 0 && result->length > 1) {
        result->length--;
    }

    free(x);
}

void test_fft() {
    // Initialize BigIntegers a and b with values
    BigInt a, b, result;
    init_bigint_value(&a, 123);
    init_bigint_value(&b, 456);

    // Multiply BigIntegers a and b using FFT
    if(g_nUseOMP == USE_OMP) {
      multiply_bigint_fft_omp(&a, &b, &result);
    } else {
      multiply_bigint_fft(&a, &b, &result);
    }

    // Print the result
    logging(LOG_DEBUG, "Result of multiplication: ");
    print_bigint(&result);

    if(g_nUseOMP == USE_OMP) {
      square_bigint_fft_omp(&a, &result);
    } else {
      // Square BigInt a using FFT
      square_bigint_fft(&a, &result);
    }

    logging(LOG_DEBUG, "Result of squaring: ");
    print_bigint(&result);

}

void mt19937_initialize(mt19937_state *state, uint32_t seed) {
    state->state[0] = seed;
    for (int i = 1; i < MT_N; ++i) {
        state->state[i] = (1812433253UL * (state->state[i-1] ^ (state->state[i-1] >> 30)) + i);
    }
    state->index = MT_N;
}

uint32_t mt19937_extract_number(mt19937_state *state) {
    if (state->index >= MT_N) {
        int i;
        for (i = 0; i < MT_N - MT_M; ++i) {
            uint32_t y = (state->state[i] & MT_UPPER_MASK) | (state->state[i+1] & MT_LOWER_MASK);
            state->state[i] = state->state[i + MT_M] ^ (y >> 1) ^ (-(int32_t)(y & 1) & MT_MATRIX_A);
        }
        for (; i < MT_N - 1; ++i) {
            uint32_t y = (state->state[i] & MT_UPPER_MASK) | (state->state[i+1] & MT_LOWER_MASK);
            state->state[i] = state->state[i + (MT_M - MT_N)] ^ (y >> 1) ^ (-(int32_t)(y & 1) & MT_MATRIX_A);
        }
        uint32_t y = (state->state[MT_N - 1] & MT_UPPER_MASK) | (state->state[0] & MT_LOWER_MASK);
        state->state[MT_N - 1] = state->state[MT_M - 1] ^ (y >> 1) ^ (-(int32_t)(y & 1) & MT_MATRIX_A);
        state->index = 0;
    }

    uint32_t y = state->state[state->index++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);
    return y;
}
uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

// Function to calculate (a^b) % mod
uint64_t power_ll(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result = 1;
    a = a % mod;
    while (b > 0) {
        if (b & 1)
            result = (result * a) % mod;
        b = b >> 1;
        a = (a * a) % mod;
    }
    return result;
}

bool miller_rabin(uint64_t n, int k) {
    if (n <= 1 || n == 4) {
        return false;
    }
    if (n <= 3) {
        return true;
    }

    uint64_t d = n - 1;
    while (d % 2 == 0) {
        d /= 2;
    }

    for (int i = 0; i < k; i++) {
        uint64_t a = mt19937_extract_number(&mt_state) % (n - 4) + 2; // Generate a random base between [2, n - 2]
        uint64_t x = power(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }
        bool prime = false;
        for (uint64_t r = d; r != n - 1; r *= 2) {
            x = (x * x) % n;
            if (x == 1) {
                return false;
            }
            if (x == n - 1) {
                prime = true;
                break;
            }
        }
        if (!prime) {
            return false;
        }
    }
    return true;
}

void test_miller_rabin() {
    uint64_t num = 2147483647ULL; // Number to be tested for primality
    int iterations = 5; // Number of iterations for the Miller-Rabin test

    mt19937_initialize(&mt_state, omp_get_thread_num()); // Initialize Mersenne Twister with thread-specific seed

    if (miller_rabin(num, iterations)) {
        logging(LOG_DEBUG, "%llu is prime.\n", num);
    } else {
        logging(LOG_DEBUG, "%llu is composite.\n", num);
    }
}

void lucas_lehmer_test(const int p) {
    BigInt s, m;
    init_bigint(&s);
    init_bigint(&m);

    // Initialize s = 4
    s.digits[0] = 4;
    s.length = 1;

    // Calculate Mersenne number M = 2^p - 1
    m.digits[0] = 1;
    m.length = 1;
    left_shift_bigint(&m, p);
    subtract_bigint(&m, &s, &m);

    // Perform p - 2 iterations of the Lucas-Lehmer test
    for (int i = 0; i < p - 2; ++i) {
        BigInt square;
        multiply_bigint(&s, &s, &square);
        right_shift_bigint(&square, p);
        subtract_bigint(&square, &m, &s);
    }

    // If s is zero, M is prime
    if (s.length == 0) {
        logging(LOG_DEBUG, "2^%d - 1 is prime.\n", p);
    } else {
        logging(LOG_DEBUG, "2^%d - 1 is not prime.\n", p);
    }
}

void lucas_lehmer_test_omp(const int p) {
    BigInt m, s, result;
    init_bigint(&m);
    init_bigint(&s);
    init_bigint(&result);

    // Mersenne number M_p = 2^p - 1
    m.digits[0] = 1;
    m.length = 1;
    left_shift_bigint(&m, p);
    subtract_bigint(&m, &(BigInt){{1}, 1}, &m);

    // Initialize s = 4
    s.digits[0] = 4;
    s.length = 1;

    // Perform p - 2 iterations of squaring and subtracting 2
    #pragma omp parallel for
    for (int i = 0; i < p - 2; i++) {
        square_bigint(&s, &result);  // Square s
        subtract_bigint(&result, &(BigInt){{2}, 1}, &s);  // Subtract 2
    }

    // Check if s is divisible by M_p
    mod_bigint(&s, &m, &result);
    if (result.length == 1 && result.digits[0] == 0) {
        logging(LOG_DEBUG, "Mersenne prime M_%d = 2^%d - 1 is prime\n", p, p);
    } else {
        logging(LOG_DEBUG, "Mersenne number M_%d = 2^%d - 1 is not prime\n", p, p);
    }
}

void test_lucaslehmer() {
    int exponent = 9; // Exponent of the Mersenne prime to test: 2^exponent - 1

    if(g_nUseOMP == USE_OMP) {
      lucas_lehmer_test_omp(exponent);

    } else {
      lucas_lehmer_test(exponent);
    }
}

#ifdef BIGINT
void test_bigint() {
    BigInt a, b, qt, rm, result;
    init_bigint(&a);
    init_bigint(&b);
    init_bigint(&qt);
    init_bigint(&rm);
    init_bigint(&result);

    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});

    logging(LOG_DEBUG, "a = ");
    print_bigint(&a);
    logging(LOG_DEBUG,"b = ");
    print_bigint(&b);

    logging(LOG_DEBUG, "a + b = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    add_bigint(&a, &b, &result);
    print_bigint(&result);

    logging(LOG_DEBUG, "a - b = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    subtract_bigint(&a, &b, &result);
    print_bigint(&result);

    logging(LOG_DEBUG, "a / b = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    divide_bigint(&a, &b, &qt, &rm);
    print_bigint(&qt);
    print_bigint(&rm);

    logging(LOG_DEBUG,"a * b = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    multiply_bigint(&a, &b, &result);
    print_bigint(&result);

    logging(LOG_DEBUG, "a mod b = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    mod_bigint(&a, &b, &result);
    print_bigint(&result);

    logging(LOG_DEBUG, "a^2 = ");
    assign_bigint(&a, &(BigInt){{ 1, 1, 1}, 3});
    assign_bigint(&b, &(BigInt){{ 1, 1}, 2});
    square_bigint(&a, &result);
    print_bigint(&a);
    print_bigint(&result);

    logging(LOG_DEBUG, "sqrt(a=12) = ");
    assign_bigint(&a, &(BigInt){{ 4, 4, 1}, 3});
    sqrt_bigint(&a, &result);
    print_bigint(&result);

    logging(LOG_DEBUG,"(a=12)^(n=2) = ");
    assign_bigint(&a, &(BigInt){{ 2, 1}, 2});
    power_bigint(&a, 2, &result);
    print_bigint(&result);

    logging(LOG_DEBUG, "(n=5)! = ");
    factorial_bigint(5, &result);
    print_bigint(&result);

    BigInt num;
    num.length = 3;
    num.digits[0] = 0xFFFFFFFF; // Maximum value for 32-bit unsigned integer
    num.digits[1] = 0x00000001; // Least significant digit
    num.digits[2] = 0x00000000; // Additional digit

    logging(LOG_DEBUG, "Original BigInt: ");
    for (int i = num.length - 1; i >= 0; i--) {
        logging(LOG_DEBUG, "%08X", num.digits[i]);
    }
    logging(LOG_DEBUG,"\n");

    // Left shift by 36 bits
    left_shift_bit_bigint(&num, 36);
    logging(LOG_DEBUG, "Left Shifted Bits By BigInt: ");
    for (int i = num.length - 1; i >= 0; i--) {
        logging(LOG_DEBUG, "%08X", num.digits[i]);
    }
    logging(LOG_DEBUG, "\n");

    num.length = 3;
    num.digits[0] = 0xFFFFFFFF; // Maximum value for 32-bit unsigned integer
    num.digits[1] = 0x00000001; // Least significant digit
    num.digits[2] = 0x00000000; // Additional digit

    logging(LOG_DEBUG, "Original BigInt: ");
    for (int i = num.length - 1; i >= 0; i--) {
        logging(LOG_DEBUG, "%08X", num.digits[i]);
    }
    logging(LOG_DEBUG, "\n");

    // Right shift by 36 bits
    right_shift_bit_bigint(&num, 36);
    logging(LOG_DEBUG,"Right Shifted Bits By BigInt: ");
    for (int i = num.length - 1; i >= 0; i--) {
        logging(LOG_DEBUG, "%08X", num.digits[i]);
    }
    logging(LOG_DEBUG, "\n");

    //BigInt num;
    init_bigint(&num);
    num.digits[0] = 1;
    num.digits[1] = 2;
    num.digits[2] = 3;
    num.length = 3;

    logging(LOG_DEBUG, "Before left shift: ");
    print_bigint(&num);

    left_shift_bigint(&num, 1);
    logging(LOG_DEBUG, "After left shift: ");
    print_bigint(&num);

    logging(LOG_DEBUG, "Before right shift: ");
    print_bigint(&num);

    right_shift_bigint(&num, 1);
    logging(LOG_DEBUG, "After right shift: ");
    print_bigint(&num);

    double angle_degrees = 45.0;  // Angle in degrees
    double angle_radians = angle_degrees * M_PI / 180.0; // Convert angle to radians
    int precision = PRECISION;  // Desired precision

    // Compute and print sine
    double sin_value = sine(angle_radians, precision);
    logging(LOG_DEBUG, "Sine of %.2f degrees with precision %d is %.50f\n", angle_degrees, precision, sin_value);

    // Compute and print cosine
    double cos_value = cosine(angle_radians, precision);
    logging(LOG_DEBUG, "Cosine of %.2f degrees with precision %d is %.50f\n", angle_degrees, precision, cos_value);

    // Compute and print natural logarithm
    double log_value = logarithm(2.0, precision);
    logging(LOG_DEBUG, "Natural logarithm of 2.0 with precision %d is %.50f\n", precision, log_value);

    // Compute and print base 10 logarithm
    double log10_value = logarithm_base10(2.0, precision);
    logging(LOG_DEBUG, "Base 10 logarithm of 2.0 with precision %d is %.50f\n", precision, log10_value);

    // Compute and print pi
    double pi_value = pi(precision);
    logging(LOG_DEBUG, "Value of pi with precision %d is %.50f\n", precision, pi_value);

    compute_pi();
}

// Function to initialize a BigInt with a string representation of a number
void init_bigint_from_string(const char *num_str, BigInt *num) {
    num->length = strlen(num_str);
    for (int i = 0; i < num->length; i++) {
        num->digits[i] = num_str[num->length - i - 1] - '0';
    }
}

void init_bigint(BigInt *num) {
    memset(num->digits, 0, MAX_DIGITS * sizeof(int));
    num->length = 0;
}

void assign_bigint(BigInt *dest, const BigInt *src) {
    memcpy(dest->digits, src->digits, MAX_DIGITS * sizeof(int));
    dest->length = src->length;
}

void print_bigint(const BigInt *num) {
    if(num->length -1 < 0) {
      printf("[]\n");
    } else {
    for (int i = num->length - 1; i >= 0; i--) {
        printf("[%d]", num->digits[i]);
    }
      printf("\n");
    }
}

void print_bigint_desc(char * desc, const BigInt *num) {
    printf("[%s]:", desc);
    if(num->length -1 < 0) {
      printf("[]\n");
    } else {
      for (int i = num->length - 1; i >= 0; i--) {
        printf("[%d]", num->digits[i]);
      }
      printf("\n");
    }
}

// Function to compare two big integers (returns -1 if a < b, 0 if a == b, and 1 if a > b)
int compare_bigint(const BigInt *a,const BigInt *b) {
    // Compare the lengths of the numbers
    if (a->length < b->length) {
        return -1;
    } else if (a->length > b->length) {
        return 1;
    } else {
        // Compare digits starting from the most significant digit
        for (int i = a->length - 1; i >= 0; i--) {
            if (a->digits[i] < b->digits[i]) {
                return -1;
            } else if (a->digits[i] > b->digits[i]) {
                return 1;
            }
        }
        // If all digits are equal, return 0
        return 0;
    }
}

// Function to add two big integers a and b and store the result in result
void add_bigint(const BigInt *a,const BigInt *b, BigInt *result) {
    int carry = 0;
    int max_length = (a->length > b->length) ? a->length : b->length;

    for (int i = 0; i < max_length; i++) {
        int sum = carry;
        if (i < a->length) {
            sum += a->digits[i];
        }
        if (i < b->length) {
            sum += b->digits[i];
        }
        result->digits[i] = sum % 10;
        carry = sum / 10;
    }

    if (carry > 0) {
        result->digits[max_length] = carry;
        result->length = max_length + 1;
    } else {
        result->length = max_length;
    }
}

void subtract_bigint(const BigInt *a,const BigInt *b, BigInt *result) {
    int borrow = 0;
    int i;
    for (i = 0; i < a->length; i++) {
        int diff = a->digits[i] - borrow;
        if (i < b->length) diff -= b->digits[i];
        if (diff < 0) {
            diff += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        result->digits[i] = diff;
    }
    result->length = a->length;
    while (result->digits[result->length - 1] == 0 && result->length > 1) {
        result->length--;
    }
}

// Function to multiply two big integers
void multiply_bigint(const BigInt *a,const BigInt *b, BigInt *result) {
    // Initialize result to zero
    memset(result->digits, 0, sizeof(result->digits));
    result->length = 0;

    // Perform multiplication
    for (int i = 0; i < a->length; i++) {
        int carry = 0;
        for (int j = 0; j < b->length || carry; j++) {
            int tmp = result->digits[i + j] + a->digits[i] * (j < b->length ? b->digits[j] : 0) + carry;
            result->digits[i + j] = tmp % BASE;
            carry = tmp / BASE;
            if (i + j + 1 > result->length && (result->digits[i + j] != 0 || carry)) {
                result->length = i + j + 1;
            }
        }
    }
}

void multiply_scalar_bigint(const BigInt *a,const int scalar, BigInt *result) {
    int carry = 0;

    // Initialize the result to zero
    init_bigint(result);

    // Iterate over each digit of BigInt 'a'
    for (int i = 0; i < a->length || carry; ++i) {
        int product = carry;
        if (i < a->length) {
            product += a->digits[i] * scalar;
        }
        result->digits[i] = product % 10; // Store the least significant digit
        carry = product / 10; // Update carry for the next iteration
    }

    // Set the length of the result BigInt
    result->length = a->length + 1; // The length might increase by one due to carry

    // Remove leading zeros from the result BigInt
    while (result->length > 1 && result->digits[result->length - 1] == 0) {
        result->length--;
    }
}

void divide_bigint(const BigInt *a, const BigInt *b, BigInt *quotient, BigInt *remainder) {
    BigInt temp; // Temporary BigInt for storing intermediate results
    init_bigint(&temp);
    // Initialize quotient and remainder BigInts
    init_bigint(quotient);
    init_bigint(remainder);

    // Loop through each digit of the dividend 'a' starting from the most significant digit
    for (int i = a->length - 1; i >= 0; --i) {
        // Shift the remainder to the left by one digit and add the current digit of the dividend
        init_bigint(&temp);
        multiply_scalar_bigint(remainder, 10, &temp);
        assign_bigint(remainder, &temp);
        remainder->digits[0] = a->digits[i];
        if(remainder->length <= 0) {
          remainder->length = 1;
        }

        // Initialize the quotient digit to zero
        int qDigit = 0;

        // Binary search for the quotient digit that makes the current remainder less than or equal to the divisor 'b'
        int low = 0, high = 9;
        while (low <= high) {
            int mid = (low + high) / 2;

            // Multiply the divisor 'b' by the current quotient digit 'mid'
            multiply_scalar_bigint(b, mid, &temp);

            // Compare the result with the current remainder
            if (compare_bigint(&temp, remainder) <= 0) {
                qDigit = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }

        // Subtract the product of the divisor 'b' and the quotient digit 'qDigit' from the remainder
        multiply_scalar_bigint(b, qDigit, &temp);
        subtract_bigint(remainder, &temp, remainder);

        // Set the current digit of the quotient
        quotient->digits[i] = qDigit;
        if(quotient->length < i + 1) {
          quotient->length = i+1;
        }
    }

    // Remove leading zeros from the quotient
    quotient->length = a->length;
    while (quotient->length > 1 && quotient->digits[quotient->length - 1] == 0) {
        quotient->length--;
    }

    // Remove leading zeros from the remainder
    while (remainder->length > 1 && remainder->digits[remainder->length - 1] == 0) {
        remainder->length--;
    }
}

void divide_scalar_bigint(const BigInt *a, const int scalar, BigInt *quotient) {
    // Initialize quotient
    init_bigint(quotient);

    int remainder = 0;
    for (int i = a->length - 1; i >= 0; --i) {
        int dividend = remainder * 10 + a->digits[i];
        quotient->digits[i] = dividend / scalar;
        remainder = dividend % scalar;
    }

    // Adjust the length of the quotient
    quotient->length = a->length;
    while (quotient->length > 1 && quotient->digits[quotient->length - 1] == 0) {
        --quotient->length;
    }
}

// Function to perform modular arithmetic with big integers (a % b)
void mod_bigint(const BigInt *a, const BigInt *b, BigInt *result) {
    BigInt temp;
    BigInt temp2;
    init_bigint(&temp);
    init_bigint(&temp2);
    BigInt quotient;
    if (compare_bigint(a, b) < 0) {
        // If a < b, result is a itself
        *result = *a;
    } else {
        // Otherwise, perform long division to compute a % b
        temp = *a;
        quotient.length = 0;
        while (compare_bigint(&temp, b) >= 0) {
            int num_shifts = temp.length - b->length;
            BigInt shifted_b = *b;
            for (int i = shifted_b.length - 1; i >= 0; i--) {
                shifted_b.digits[i + num_shifts] = shifted_b.digits[i];
            }
            for (int i = 0; i < num_shifts; i++) {
                shifted_b.digits[i] = 0;
            }
            shifted_b.length += num_shifts;
            while (compare_bigint(&temp, &shifted_b) >= 0) {
                init_bigint(&temp2);
                subtract_bigint(&temp, &shifted_b, &temp2);
                assign_bigint(&temp, &temp2);
                quotient.digits[quotient.length++] = 1;
            }
        }
        *result = temp;
    }
}

void square_bigint(const BigInt *a, BigInt *result) {
    multiply_bigint(a, a, result);
}

void sqrt_bigint(const BigInt *a, BigInt *result) {
    BigInt low, high, mid, mid_squared;
    BigInt one, two, temp;

    init_bigint(&low);
    init_bigint(&high);
    init_bigint(&mid);
    init_bigint(&mid_squared);
    init_bigint(&one);
    init_bigint(&two);
    init_bigint(&temp);

    // Initialize constants
    one.digits[0] = 1;
    one.length = 1;
    two.digits[0] = 2;
    two.length = 1;

    // Set initial range for binary search
    assign_bigint(&low, &one);
    assign_bigint(&high, a);

    // Binary search for the square root
    while (compare_bigint(&low, &high) <= 0) {
        // Compute mid = (low + high) / 2
        add_bigint(&low, &high, &temp);
        divide_scalar_bigint(&temp, 2, &mid);

        // Compute mid_squared = mid * mid
        square_bigint(&mid, &mid_squared);

        // Compare mid_squared with a
        int cmp = compare_bigint(&mid_squared, a);
        if (cmp == 0) {
            // Found exact square root
            assign_bigint(result, &mid);
            return;
        } else if (cmp < 0) {
            // Mid is too low, adjust the range
            assign_bigint(&low, &mid);
            init_bigint(&temp);
            add_bigint(&low, &one, &temp);
            assign_bigint(&low, &temp);
        } else {
            // Mid is too high, adjust the range
            assign_bigint(&high, &mid);
            init_bigint(&temp);
            subtract_bigint(&high, &one, &temp);
            assign_bigint(&high, &temp);
        }
    }

    // The result is the lower bound of the range
    assign_bigint(result, &low);
}

void power_bigint(const BigInt *base, const int exponent, BigInt *result) {
    BigInt temp;
    init_bigint(&temp);
    init_bigint(result);
    result->digits[0] = 1;
    if(result->length <= 0) {
      result->length = 1;
    }
    for (int i = 0; i < exponent; i++) {
        multiply_bigint(result, base, &temp);
        assign_bigint(result, &temp);
    }
}

// Function to compute factorial of n using big integers
void factorial_bigint(const int n, BigInt *result) {
    BigInt temp;
    BigInt temp2;
    // Initialize result to 1
    init_bigint(result);
    result->digits[0] = 1;
    if(result->length < 1) {
      result->length = 1;
    }

    // Compute factorial
    for (int i = 2; i <= n; i++) {
        init_bigint(&temp);
        init_bigint(&temp2);
        temp.digits[0] = i;
        temp.length++;
        assign_bigint(&temp2, result);
        multiply_bigint(&temp, &temp2, result);
    }
}

// Function to perform left shift operation on a BigInt by a specified number of bits
void left_shift_bit_bigint(BigInt *num,const int shiftbit) {
    if (shiftbit <= 0) return;

    // Calculate the number of digit shifts and bit shifts
    int digit_shift = shiftbit / BITS_PER_DIGIT;
    int bit_shift = shiftbit % BITS_PER_DIGIT;

    // Shift digits to the left by digit_shift
    for (int i = num->length - 1; i >= 0; i--) {
        num->digits[i + digit_shift] = num->digits[i];
    }

    // Fill the shifted digits with zeros
    for (int i = 0; i < digit_shift; i++) {
        num->digits[i] = 0;
    }

    // If there are remaining bit shifts, perform bitwise left shift
    if (bit_shift > 0) {
        uint32_t carry = 0;
        for (int i = 0; i < num->length; i++) {
            uint32_t next_carry = (num->digits[i] >> (BITS_PER_DIGIT - bit_shift));
            num->digits[i] <<= bit_shift;
            num->digits[i] |= carry;
            carry = next_carry;
        }
    }

    // Update length
    num->length += digit_shift;
    if (num->digits[num->length - 1] == 0 && num->length > 1) {
        num->length--;
    }
}

// Function to perform right shift operation on a BigInt by a specified number of bits
void right_shift_bit_bigint(BigInt *num,const int shiftbit) {
    if (shiftbit <= 0) return;

    // Calculate the number of digit shifts and bit shifts
    int digit_shift = shiftbit / BITS_PER_DIGIT;
    int bit_shift = shiftbit % BITS_PER_DIGIT;

    if(digit_shift*BITS_PER_DIGIT + bit_shift <= shiftbit) {
      init_bigint(num);
      num->length = 1;
      return;
    }
    // Shift digits to the right by digit_shift
    for (int i = 0; i < num->length; i++) {
        num->digits[i - digit_shift] = num->digits[i];
    }

    // Fill the shifted digits with zeros
    for (int i = num->length - 1; i >= num->length - digit_shift; i--) {
        num->digits[i] = 0;
    }

    // If there are remaining bit shifts, perform bitwise right shift
    if (bit_shift > 0) {
        uint32_t carry = 0;
        for (int i = num->length - 1; i >= 0; i--) {
            uint32_t next_carry = (num->digits[i] << (BITS_PER_DIGIT - bit_shift));
            num->digits[i] >>= bit_shift;
            num->digits[i] |= carry;
            carry = next_carry;
        }
    }

    // Update length
    num->length -= digit_shift;
    printf("num->length [%d], num->digits[num->length - 1] [%d]\n", num->length, num->digits[num->length -1]);
    if (num->digits[num->length - 1] == 0 && num->length > 1) {
        num->length--;
    }
}


void left_shift_bigint(BigInt *num, const int shiftdigits) {
    // Shift each digit to the left by shiftdigits positions
    for (int i = num->length - 1; i >= 0; i--) {
        num->digits[i + shiftdigits] = num->digits[i];
    }

    // Fill the vacated positions with zeros
    for (int i = 0; i < shiftdigits; i++) {
        num->digits[i] = 0;
    }

    // Update the length of the BigInt
    num->length += shiftdigits;
}

void right_shift_bigint(BigInt *num, const int shiftdigits) {
    // Shift each digit to the right by shiftdigits positions
    for (int i = 0; i < num->length - shiftdigits; i++) {
        num->digits[i] = num->digits[i + shiftdigits];
    }

    // Fill the vacated positions with zeros
    for (int i = num->length - shiftdigits; i < num->length; i++) {
        num->digits[i] = 0;
    }

    // Update the length of the BigInt
    num->length -= shiftdigits;
}

// Function to compute factorial
double factorial(int n) {
    if (n <= 1) return 1;
    else return n * factorial(n - 1);
}

// Function to compute sine using Taylor series expansion
double sine(double x, int precision) {
    double result = 0;
    int sign = 1;
    for (int n = 0; n < precision; n++) {
        result += sign * pow(x, 2 * n + 1) / factorial(2 * n + 1);
        sign *= -1;
    }
    return result;
}

// Function to compute cosine using Taylor series expansion
double cosine(double x, int precision) {
    double result = 0;
    int sign = 1;
    for (int n = 0; n < precision; n++) {
        result += sign * pow(x, 2 * n) / factorial(2 * n);
        sign *= -1;
    }
    return result;
}

// Function to compute natural logarithm using Taylor series expansion
double logarithm(double x, int precision) {
    if (x <= 0) return NAN; // Handle non-positive values
    double result = 0;
    for (int n = 1; n <= precision; n++) {
        result += pow((x - 1) / x, n) / n;
    }
    return result;
}


// Function to compute base 10 logarithm using Taylor series expansion
double logarithm_base10(double x, int precision) {
    return logarithm(x, precision) / log(10);
}

// Function to compute pi using Machin's formula and Taylor series expansion
double pi(int precision) {
    double result = 0;
    for (int n = 0; n < precision; n++) {
        result += pow(-1, n) / (2 * n + 1);
    }
    return 4 * result;
}


void multiply_pi(int *arr, int multiplier) {
    int carry = 0;
    for (int i = PRECISION - 1; i >= 0; i--) {
        int product = arr[i] * multiplier + carry;
        arr[i] = product % 10;
        carry = product / 10;
    }
}

void print_pi(int *arr) {
    printf("π = 3.");
    for (int i = 1; i <= PRECISION; i++) {
        printf("%d", arr[i - 1]);
    }
}

void compute_pi() {
    int pi[PRECISION];
    for (int i = 0; i < PRECISION; i++) {
        pi[i] = 2;
    }

    int factor = 10;
    int temporary[PRECISION];
    for (int i = 0; i < PRECISION; i++) {
        temporary[i] = 0;
    }

    for (int i = 0; i < PRECISION * 10; i++) {
        multiply_pi(temporary, factor);
        temporary[PRECISION - 1] = 1;
        for (int j = 0; j < PRECISION; j++) {
            pi[j] += temporary[j];
        }
        factor += 2;
    }

    print_pi(pi);
}

void mod_multiply(const BigInt *a, const BigInt *b, const BigInt *n, BigInt *result) {
    BigInt temp;
    init_bigint(&temp);

    for (int i = 0; i < b->length; i++) {
        for (int j = 0; j < b->digits[i]; j++) {
            BigInt partial;
            init_bigint(&partial);
            partial = *a;  // Copy a
            for (int k = 0; k < i; k++) {
                left_shift_bigint(&partial, 1);  // Multiply by 10 for each digit of b
            }
            add_bigint(&temp, &partial, &temp);  // Add partial result to temp
        }
    }

    mod_bigint(&temp, n, result);  // Take the modulus of temp with n
}

// Function to calculate the greatest common divisor (gcd) of two integers
int gcd_bigint(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}   

    
// Function to calculate binomial coefficient
uint64_t binomial(uint64_t n, uint64_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
        
    uint64_t res = 1;
    for (uint64_t i = 0; i < k; i++) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res; 
}       
    
// Function to calculate modular exponentiation
uint64_t mod_exp(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

// Function to initialize a big integer with a value
void init_bigint_value(BigInt *num, int value) {
    num->length = 0;
    while (value > 0) {
        num->digits[num->length++] = value % 10;
        value /= 10;
    }
}

#endif // #ifdef BIGINT

void test_montgomery_multiply(){
    BigInt a, b, n, n_inv, r, result;
    init_bigint(&a);
    init_bigint(&b);
    init_bigint(&n);
    init_bigint(&n_inv);
    init_bigint(&r);
    init_bigint(&result);

    // Initialize values for demonstration
    a.digits[0] = 12345;
    a.length = 1;
    b.digits[0] = 67890;
    b.length = 1;
    n.digits[0] = 10000;
    n.length = 1;
    n_inv.digits[0] = 9000; // This should be the modular inverse of n
    n_inv.length = 1;
    r.digits[0] = 100000;
    r.length = 1;

    montgomery_multiply_omp(&a, &b, &n, &n_inv, &r, &result);
    print_bigint(&result);
}

void montgomery_multiply(BigInt *a, BigInt *b, BigInt *n, BigInt *n_inv, BigInt *r, BigInt *result) {
    BigInt temp;
    init_bigint(&temp);
    init_bigint(result);

    for (int i = 0; i < a->length; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < b->length; j++) {
            uint64_t product = a->digits[i] * b->digits[j] + temp.digits[i + j] + carry;
            temp.digits[i + j] = product % 10;
            carry = product / 10;
        }
        temp.digits[i + b->length] = carry;
    }

    for (int i = 0; i < a->length + b->length; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < n->length; j++) {
            uint64_t product = result->digits[i] + temp.digits[i] * n_inv->digits[j] + carry;
            result->digits[i] = product % 10;
            carry = product / 10;
        }
        result->digits[i + n->length] += carry;
    }

    if (compare_bigint(result, n) >= 0) {
        init_bigint(&temp);
        subtract_bigint(result, n, &temp);
        // Update result
        *result = temp;
    }
}

void montgomery_multiply_omp(BigInt *a, BigInt *b, BigInt *n, BigInt *n_inv, BigInt *r, BigInt *result) {
    BigInt temp;
    init_bigint(&temp);
    init_bigint(result);

    // Perform Montgomery multiplication in parallel
    #pragma omp parallel for shared(a, b, n, n_inv, r, temp)
    for (int i = 0; i < a->length; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < b->length; j++) {
            uint64_t product = a->digits[i] * b->digits[j] + temp.digits[i + j] + carry;
            temp.digits[i + j] = product % 10;
            carry = product / 10;
        }
        temp.digits[i + b->length] = carry;
    }

    // Combine results from parallel threads
    #pragma omp barrier
    // Reduction step can be performed here to merge partial results
    // Assign result
    for (int i = 0; i < a->length + b->length; i++) {
        result->digits[i] = temp.digits[i];
    }
    result->length = a->length + b->length;
}

void rsa_encrypt(BigInt *plaintext, BigInt *n, BigInt *e, BigInt *ciphertext) {
    montgomery_multiply_omp(plaintext, e, n, NULL, NULL, ciphertext);
}

void rsa_decrypt(BigInt *ciphertext, BigInt *n, BigInt *d, BigInt *plaintext) {
    montgomery_multiply_omp(ciphertext, d, n, NULL, NULL, plaintext);
}

void test_montgomery_rsa(){
    BigInt plaintext, ciphertext, n, e, d;
    init_bigint(&plaintext);
    init_bigint(&ciphertext);
    init_bigint(&n);
    init_bigint(&e);
    init_bigint(&d);

    // Encrypt plaintext
    rsa_encrypt(&plaintext, &n, &e, &ciphertext);
    printf("Ciphertext: ");
    // Print ciphertext

    // Decrypt ciphertext
    rsa_decrypt(&ciphertext, &n, &d, &plaintext);
    printf("Plaintext: ");
    // Print plaintext
}

// Function to calculate the Montgomery zeta function
double montgomery_zeta(double x) {
    // Implement the Montgomery zeta function
    double result = 0.0;
    for (double n = 1.0; n <= x; ++n) {
        result += 1.0 / pow(n, 0.5);
    }
    return result;
}

// Function to calculate the prime count using Montgomery's prime counting function
uint64_t count_primes_montgomery(uint64_t n) {
    // Calculate the integral term
    double integral_term = 0.0;
    for (double t = n; t <= 2 * n; t += 0.01) {
        integral_term += 1.0 / (t * (t * t - 1) * log(t));
    }

    // Calculate the sum over non-trivial zeros of the Riemann zeta function
    double sum_over_zeros = 0.0;
    for (double rho = 0.5; rho <= 100; rho += 0.5) {
        sum_over_zeros += montgomery_zeta(pow(n, rho));
    }

    // Calculate the prime count using Montgomery's prime counting function
    uint64_t prime_count = round(montgomery_zeta(n) - sum_over_zeros - log(2) / 2 + integral_term);
    return prime_count;
}

// Function to check if a number is prime using AKS primality test
bool is_prime_aks(uint64_t n) {
    // Check for some base cases
    if (n <= 1) return false;
    if (n <= 3) return true;

    // Find the smallest r such that ord_p(r) > log2(n)
    uint64_t r = 2;
    while (r < n) {
        if (gcd_bigint(r, n) != 1) return false;
        uint64_t k = 2;
        uint64_t ord = 1;
        while (ord < n) {
            if (mod_exp(r, ord, n) == 1) break;
            ord = k * (n - 1);
            k++;
        }
        if (ord > log2(n)) break;
        r++;
    }

    // Check if (X - a)^n ≡ (X^n - a) (mod X^r - 1, n)
    for (uint64_t a = 1; a <= sqrt(n); a++) {
        for (uint64_t b = 1; b <= sqrt(n); b++) {
            uint64_t lhs = binomial(n, a) * mod_exp(b, n, n);
            uint64_t rhs = binomial(n, a) * mod_exp(b, n, n) + mod_exp(n, a, n);
            if (lhs % (n * n) != rhs % (n * n)) return false;
        }
    }

    // Check if n is a prime power
    for (uint64_t a = 2; a <= sqrt(n); a++) {
        if ((n % a == 0) && is_prime_aks(a)) return false;
    }

    return true;
}

void test_aks() {
    uint64_t num = 127;
    if(g_nUseOMP == USE_OMP) {
      // Set number of threads for OpenMP
      omp_set_num_threads(10);
      // Check if n is prime using AKS primality test
      bool result = is_prime_aks_omp(num);

      // Output the result
      if (result)
        printf("%d is prime.\n", num);
      else
        printf("%d is not prime.\n", num);
    } else {
      if (is_prime_aks(num)) {
        printf("%llu is prime.\n", num);
      } else {
        printf("%llu is composite.\n", num);
      }
    }
}

// AKS primality test function
bool is_prime_aks_omp(uint64_t n) {
    // Base cases
    if (n <= 1) {
        return false;
    }
    if (n <= 3) {
        return true;
    }
    if (n % 2 == 0 || n % 3 == 0) {
        return false;
    }

    // Initialize BigInt for comparison
    BigInt num, factor;
    init_bigint_value(&num, n);

    // AKS primality test algorithm
    //#pragma omp parallel for private(factor)
    for (uint64_t r = 2; r * r <= n; r++) {
        init_bigint_value(&factor, r);
        if (compare_bigint(&factor, &num) == 0) {
            continue;
        }
        BigInt exp;
        exp = factor;
        while (compare_bigint(&exp, &num) < 0) {
            int a = 1;
            for (uint64_t m = 1; m <= n; m++) {
                a = (a * r) % (int)m;
                if (compare_bigint(&factor, &num) != 0 && compare_bigint(&exp, &num) == 0 && a != 1) {
                    return false;
                }
            }
            init_bigint_value(&exp, exp.length + 1);
        }
    }

    return true;
}

// Function to count prime numbers up to n using the Sieve of Eratosthenes algorithm
uint64_t count_primes_eratosthenes_omp(uint64_t n) {
    if (n <= 1) return 0;

    // Allocate memory to store flags indicating prime or composite
    bool *is_prime = (bool *)malloc((n + 1) * sizeof(bool));
    if (is_prime == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize all numbers as prime
    #pragma omp parallel for
    for (uint64_t i = 0; i <= n; ++i) {
        is_prime[i] = true;
    }

    // Apply Sieve of Eratosthenes algorithm
    //#pragma omp parallel for schedule(dynamic)
    for (uint64_t p = 2; p * p <= n; ++p) {
        if (is_prime[p]) {
            // Mark multiples of p as composite
            for (uint64_t i = p * p; i <= n; i += p) {
                is_prime[i] = false;
            }
        }
    }

    // Count prime numbers
    uint64_t count = 0;
    #pragma omp parallel for reduction(+:count)
    for (uint64_t p = 2; p <= n; ++p) {
        if (is_prime[p]) {
            ++count;
        }
    }

    free(is_prime);
    return count;
}

// Function to count prime numbers up to n
uint64_t count_primes_eratosthenes(uint64_t n) {
    if (n <= 1) return 0;

    // Allocate memory to store flags indicating prime or composite
    bool *is_prime = (bool *)malloc((n + 1) * sizeof(bool));
    if (is_prime == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize all numbers as prime
    for (uint64_t i = 0; i <= n; ++i) {
        is_prime[i] = true;
    }

    // Apply Sieve of Eratosthenes algorithm
    for (uint64_t p = 2; p * p <= n; ++p) {
        if (is_prime[p]) {
            // Mark multiples of p as composite
            for (uint64_t i = p * p; i <= n; i += p) {
                is_prime[i] = false;
            }
        }
    }

    // Count prime numbers
    uint64_t count = 0;
    for (uint64_t p = 2; p <= n; ++p) {
        if (is_prime[p]) {
            ++count;
        }
    }

    free(is_prime);
    return count;
}

uint64_t count_primes_riemann(uint64_t n) {
    if (n <= 1) return 0;

    double x = n;
    double log_x = log(x);
    double pi_x = x / log_x;

    return (uint64_t)pi_x;
}

void test_fermat_primes() {
    int count = 10; // Number of Fermat primes to generate
    fermat_primes(count);
}

void fermat_primes(int n) {
    uint64_t num;
    for (int i = 0; i < n; ++i) {
        num = pow(2, pow(2, i)) + 1;
        if (is_prime(num)) {
            printf("Fermat prime #%d: %llu\n", i+1, num);
        }
    }
}

void test_count_primes() {
     uint64_t n;
    printf("Enter a number: ");
    scanf("%llu", &n);

    if(g_nUseOMP == USE_OMP) {
      uint64_t primes_count = count_primes_eratosthenes_omp(n);

      printf("Number of prime numbers less than or equal to %llu: %llu\n", n, primes_count);
    } else {
      uint64_t primes_count = count_primes_eratosthenes(n);

      printf("Number of prime numbers less than or equal to %llu: %llu\n", n, primes_count);
    }

    uint64_t primes_count = count_primes_riemann(n);

    printf("Number of prime numbers less than %llu (estimated using Riemann Hypothesis): %llu\n", n, primes_count);


    uint64_t x = 1000; // Calculate the prime count up to x
    uint64_t count = count_primes_montgomery(x);
    printf("The number of primes less than or equal to %llu is %llu\n", x, count);
}
// Function to calculate factorial modulo n
uint64_t factorial_mod_n(uint64_t n) {
    if (n <= 1) return 1;

    uint64_t result = 1;
    for (uint64_t i = 2; i < n; ++i) {
        result = ((result * i) % n);
    }
    return result;
}

// Function to check if a number satisfies Wilson's theorem
int is_prime_wilson(uint64_t n) {
    if (n <= 1) return 0; // 0 and 1 are not prime
    // Calculate (n - 1)! modulo n
    uint64_t factorial_mod = factorial_mod_n(n);

    // Check if (n - 1)! ≡ -1 (mod n)
    return (factorial_mod == n - 1);
}

void test_wilson() {
    uint64_t number;
    printf("Enter a number to check if it's prime using Wilson's theorem: ");
    scanf("%llu", &number);

    if (is_prime_wilson(number)) {
        printf("%llu is a prime number according to Wilson's theorem.\n", number);
    } else {
        printf("%llu is not a prime number according to Wilson's theorem.\n", number);
    }
}

void test_monteCarloSimulation() {
  int num_points = 1000000;

  int inside_circle = 0;

  #pragma omp parallel
  {

     #pragma omp for reduction(+:inside_circle)
     for (int i = 0; i < num_points; ++i) {
       double x = (double)mt19937_extract_number(&mt_state) / (double)UINT32_MAX;
       double y = (double)mt19937_extract_number(&mt_state) / (double)UINT32_MAX;
       double distance = x * x + y * y;
       if (distance <= 1) {
         inside_circle++;
      }
    }
  }

  double pi_estimate = (double)inside_circle / num_points * 4;
  printf("Estimate of pi using Monte Carlo simulation: %.10f\n", pi_estimate);
}
#endif

bool isPrimeby_kj(uint64_t n, int maxNumber) {
  double s = 1.0;
  double si = 0.0;
  // The Riemann hypothesis is the conjecture that the Riemann zeta function has its zeros only at the negative even integers and complex numbers with real part 1/2. Many consider it to be the most important unsolved problem in pure mathematics.It is of great interest in number theory because it implies results about the distribution of prime numbers. (from en.wikipedia.org)
  // The P is 1/2 if riemanne hypothesis is true.
  double p = 0.5;
  // Montgomery's pair correlation conjecture is a conjecture made by Hugh Montgomery (1973) that the pair correlation between pairs of zeros of the Riemann zeta function. (from en.wikipedia.org)
  // 1 - (sin(pi*u)/(pi*u))^2 + delta(u).
  double q = 0.0;
  double u = 1.0;
  while(1) {
    for(int n = 1; n < maxNumber; n++) {
      q = 1 - pow(2, (sin(PI*u)/(PI*u)));
      s += cos(q*log(n*1.0))/pow(n, p);
      si += -sin(q*log(n*1.0))/pow(n, p);
      logging(LOG_DEBUG, "sigma(k=1,k=%d) euler_product(%.2f + %.2fi) = %.2f + %.2fi\n", n + 1, p, q, s, si);
    }
    if(s == 0 && si == 0) {
      return true; 
    }
  }
}

// Riemann Prime Counting Function
// Riemann's explicit formula for the number of primes less than a given number states that, in terms of a sum over the zeros of the Riemann zeta function, the magnitude of the oscillations of primes around their expected position is controlled by the real parts of the zeros of the zeta function. (from en.wikipedia.org)
// π(x) -li(x) = O(x^beta * logx)
// It was conjectured in the end of the 18th century by Gauss and by Legendre to be approximately.(from en.wikipedia.org)
double getPrimeCountRiemann(int n) {
  return n / log(n);
}

double getPrimeCount_kj(int n) {
  return (sqrt(pow(n,2) + pow((n/PI), 2))) / (log(n)*primeDistribution(n));
}

// The difference between primes is not a standard of gauss distribution.
// So, we can prime counting by riemann prime counting function multiply with probability distribution.
// The difference in prime is not Gaussian distribution.
// So we can calculate prime numbers count by adding the Riemann prime calculation function and the difference of primes probability distribution.
double primeDistribution(int n) {
 // C is adjust parameter.(This is test variable.
 int C = 1.1;
 return ((1 - pow(2, ((sin(log(n)*PI*C))) / log(n)*PI*C)))*(1 - pow(2, ((cos(log(n)*PI*C)) / log(n)*PI*C))) / 2;
}

void test_goldBach() {
 goldBach();
 test_prime_gap_distribution();
  
}
// Goldbach's conjecture is one of the oldest and best-known unsolved problems in number theory and all of mathematics. It states that every even natural number greater than 2 is the sum of two prime numbers. (from en.wikipedia.org)
void goldBach() {
  uint64_t nPrimePair = 0;
  uint64_t nPrimeA = 1;// 2 add
  uint64_t nPrimeB = 0;
  bool bPrimeNumber1 = false;
  bool bPrimeNumber2 = false;
  uint64_t nPrimePairCount = 0;
  uint64_t dbPrimePair = 0;
  uint64_t nPrimeCount = 0;
  uint64_t dbPrime = 0;
  uint64_t nPrimeTotal = 0;
  uint64_t nPrimePairTotal = 0;
  uint64_t nTwinPrimeCount = 0;
  if(g_nMinPrime%2 == 0) {
    g_nMinPrime++;
  }
  for(int i = 3; i < 10000; i+=2) {
    for(int j = 1; j < i; j+=2) {
      bPrimeNumber1 = isPrime(j);
      bPrimeNumber2 = isPrime(i*2-j);
      if(bPrimeNumber1 && bPrimeNumber2) {
        nPrimePair++;
        nPrimePairTotal++;
        if(j == i -1) {
          //TODO : Insert into Twin Primes database.
          nTwinPrimeCount++;
        }
      }
      if(bPrimeNumber1) {
        nPrimeA++;
      }
      if(bPrimeNumber2) {
        nPrimeB++;
      }
    }
  }
  nPrimePairCount+=nPrimePair;
  dbPrimePair+=(uint64_t)(nPrimePair*100.0/((nPrimeA + nPrimeB)*1.0));
  nPrimeCount+=(nPrimeA + nPrimeB);
  dbPrime+=(uint64_t)((nPrimeA + nPrimeB)*100.0/(((g_nMaxPrime - g_nMinPrime)*2.0)*1.0));
  nPrimeTotal++;
  logging(LOG_DEBUG, "goldBach [%llu][%llu][%llu][%llu][%llu][%llu][%llu]\n", nPrimePairCount, dbPrimePair, nPrimeCount, dbPrime, nPrimeTotal, nPrimePairTotal, nTwinPrimeCount);

  logging(LOG_INFO, "%6llu~%llu\t\t%6.2f %3.2f %6.2f %3.2f %llu\n", g_nMinPrime,
                                        g_nMaxPrime,
                                        nPrimeCount*1.0/nPrimeTotal*1.0,
                                        dbPrime*1.0/nPrimeTotal*1.0,
                                        nPrimePairCount*1.0/nPrimePairTotal*1.0,
                                        dbPrimePair*1.0/nPrimePairTotal*1.0,
                                        nTwinPrimeCount);
}

#ifdef EXPERIMENTAL
void test_prime_gap_distribution() {
    uint64_t start=1, end=10000000;
    //printf("Enter the range to compute prime gap distribution (start end): ");
    //scanf("%d %d", &start, &end);
    prime_gap_distribution(start, end);
}
// Function to count primes within a given range and calculate the gaps
void prime_gap_distribution(uint64_t start, uint64_t end) {
    uint64_t prev_prime = 0;
    uint64_t max_gap = 0;
    uint64_t *gap_count = (uint64_t *)calloc(1, sizeof(uint64_t)); // Dynamically allocate memory for gap counts

    for (uint64_t i = start; i <= end; i++) {
        if (is_prime(i)) {
            if (prev_prime != 0) {
                uint64_t gap = i - prev_prime;
                if (gap > max_gap) {
                    gap_count = (uint64_t *)realloc(gap_count, (gap + 1) * sizeof(uint64_t)); // Resize array if needed
                    for (uint64_t j = max_gap + 1; j <= gap; j++) {
                        gap_count[j] = 0;
                    }
                    max_gap = gap;
                }
                gap_count[gap]++;
            }
            prev_prime = i;
        }
    }

    printf("Prime gap distribution from %llu to %llu:\n", start, end);
    printf("Gap\tCount\n");
    for (uint64_t i = 1; i <= max_gap; i++) {
        printf("%llu\t%llu\n", i, gap_count[i]);
    }
    free(gap_count); // Free dynamically allocated memory
}

// Function to count primes within a given range and calculate the gaps
void prime_gap_distribution_omp(uint64_t start, uint64_t end) {
    uint64_t prev_prime = 0;
    uint64_t max_gap = 0;
    uint64_t *gap_count = (uint64_t *)calloc(1, sizeof(uint64_t)); // Dynamically allocate memory for gap counts

    #pragma omp parallel for
    for (uint64_t i = start; i <= end; i++) {
        if (isPrime(i)) {
            uint64_t local_max_gap = 0;
            uint64_t local_prev_prime = 0;
            if (prev_prime != 0) {
                local_max_gap = i - prev_prime;
                local_prev_prime = prev_prime;
            }
            prev_prime = i;

            #pragma omp critical
            {
                if (local_max_gap > max_gap) {
                    gap_count = (uint64_t *)realloc(gap_count, (local_max_gap + 1) * sizeof(uint64_t)); // Resize array if needed
                    for (uint64_t j = max_gap + 1; j <= local_max_gap; j++) {
                        gap_count[j] = 0;
                    }
                    max_gap = local_max_gap;
                }
                gap_count[local_max_gap]++;
            }
        }
    }

    printf("Prime gap distribution from %llu to %llu:\n", start, end);
    printf("Gap\tCount\n");
    for (uint64_t i = 1; i <= max_gap; i++) {
        printf("%llu\t%llu\n", i, gap_count[i]);
    }
    free(gap_count); // Free dynamically allocated memory
}


// Sigmoid activation function
double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}

// Derivative of the sigmoid function
double sigmoid_derivative(double x) {
    return x * (1.0 - x);
}

// Forward pass through the neural network
void forward_pass(double input[INPUT_SIZE], double hidden[HIDDEN_SIZE], double output[OUTPUT_SIZE], double weights_ih[INPUT_SIZE][HIDDEN_SIZE], double weights_ho[HIDDEN_SIZE][OUTPUT_SIZE]) {
    // Compute hidden layer values
    for (int i = 0; i < HIDDEN_SIZE; i++) {
        hidden[i] = 0;
        for (int j = 0; j < INPUT_SIZE; j++) {
            hidden[i] += input[j] * weights_ih[j][i];
        }
        hidden[i] = sigmoid(hidden[i]);
    }

    // Compute output layer values
    for (int i = 0; i < OUTPUT_SIZE; i++) {
        output[i] = 0;
        for (int j = 0; j < HIDDEN_SIZE; j++) {
            output[i] += hidden[j] * weights_ho[j][i];
        }
        output[i] = sigmoid(output[i]);
    }
}


// Backpropagation algorithm to update weights
void backpropagation(double input[INPUT_SIZE], double hidden[HIDDEN_SIZE], double output[OUTPUT_SIZE], double target[OUTPUT_SIZE], double weights_ih[INPUT_SIZE][HIDDEN_SIZE], double weights_ho[HIDDEN_SIZE][OUTPUT_SIZE]) {
    // Compute output layer errors
    double output_errors[OUTPUT_SIZE];
    for (int i = 0; i < OUTPUT_SIZE; i++) {
        output_errors[i] = target[i] - output[i];
    }

    // Compute hidden layer errors
    double hidden_errors[HIDDEN_SIZE];
    for (int i = 0; i < HIDDEN_SIZE; i++) {
        hidden_errors[i] = 0;
        for (int j = 0; j < OUTPUT_SIZE; j++) {
            hidden_errors[i] += output_errors[j] * weights_ho[i][j];
        }
    }

    // Update weights between hidden and output layers
    for (int i = 0; i < HIDDEN_SIZE; i++) {
        for (int j = 0; j < OUTPUT_SIZE; j++) {
            double delta = output_errors[j] * sigmoid_derivative(output[j]) * hidden[i];
            weights_ho[i][j] += LEARNING_RATE * delta;
        }
    }

    // Update weights between input and hidden layers
    for (int i = 0; i < INPUT_SIZE; i++) {
        for (int j = 0; j < HIDDEN_SIZE; j++) {
            double delta = hidden_errors[j] * sigmoid_derivative(hidden[j]) * input[i];
            weights_ih[i][j] += LEARNING_RATE * delta;
        }
    }
}

void test_nn() {
    // Initialize weights randomly
    double weights_ih[INPUT_SIZE][HIDDEN_SIZE];
    double weights_ho[HIDDEN_SIZE][OUTPUT_SIZE];
    for (int i = 0; i < INPUT_SIZE; i++) {
        for (int j = 0; j < HIDDEN_SIZE; j++) {
            weights_ih[i][j] = ((double)mt19937_extract_number(&mt_state) / UINT32_MAX) * 2 - 1; // Random weights between -1 and 1
        }
    }
    for (int i = 0; i < HIDDEN_SIZE; i++) {
        for (int j = 0; j < OUTPUT_SIZE; j++) {
            weights_ho[i][j] = ((double)mt19937_extract_number(&mt_state) / UINT32_MAX) * 2 - 1; // Random weights between -1 and 1
        }
    }

    // Training data (XOR problem)
    double inputs[4][INPUT_SIZE] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    double targets[4][OUTPUT_SIZE] = {{0}, {1}, {1}, {0}};

    // Train the neural network
    // Parallelize the training loop using OpenMP
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int iter = 0; iter < NUM_EPOCHS; iter++) {
        for (int i = 0; i < 4; i++) {
            double input[INPUT_SIZE];
            double hidden[HIDDEN_SIZE];
            double output[OUTPUT_SIZE];
            for (int j = 0; j < INPUT_SIZE; j++) {
                input[j] = inputs[i][j];
            }
            forward_pass(input, hidden, output, weights_ih, weights_ho);
            backpropagation(input, hidden, output, targets[i], weights_ih, weights_ho);
        }
    }

    // Test the neural network
    printf("Testing the neural network:\n");
    for (int i = 0; i < 4; i++) {
        double input[INPUT_SIZE];
        double hidden[HIDDEN_SIZE];
        double output[OUTPUT_SIZE];
        for (int j = 0; j < INPUT_SIZE; j++) {
            input[j] = inputs[i][j];
        }
        forward_pass(input, hidden, output, weights_ih, weights_ho);
        printf("Input: %d %d, Output: %f\n", (int)input[0], (int)input[1], output[0]);
    }

}

#endif
