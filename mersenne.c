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
#include <omp.h>

#ifdef _WIN32
#include <io.h>
#define F_OK 0
#define access _access
#else
#include <unistd.h>
#endif

#define MAX_SIZE          1024
#define MAX_PATH          1024
#define MAX_NOISE         1024
#define MAX_PROCESS       1024
#define MAX_BUFF_SIZE     8096
#define MAX_NUMBER       65535
#define MAX_CHAR           256
#define MAX_BIT              8
#define MAX_DIFF           256 // min 512 difference
#define MAX_PRIMES        1000
#define MAX_NOISES         100
#define MAX_DATE            26
#define MAX_INTEGER   10000000

#define MIN_CANDIDATE         37
#define MAX_CANDIDATE 1000000000
#define MIN_PRIME             37
#define MAX_PRIME     1000000000
#define MAX_CHUNK           1024

#define MERSENNE_NUMBER        0
#define MERSENNE_CANDIDATE     1
#define KNOWN_MERSENNE 51
//#define PI 3.14159265358979
#define PI 3.14159265358979323846

// Encryption
#define blockSize 16
#define KeySize 32
#define KeyExpSize 240
#define Nb 4
#define Nk 8
#define Nr 14

#define E_VALUE 3 // 65535

#define get_sb(num) (sb[(num)])
#define get_sb_inv(num) (rsb[(num)])

#define multiply(x, y)                                \
      (  ((y & 1) * x) ^                              \
      ((y>>1 & 1) * xtime(x)) ^                       \
      ((y>>2 & 1) * xtime(xtime(x))) ^                \
      ((y>>3 & 1) * xtime(xtime(xtime(x)))) ^         \
      ((y>>4 & 1) * xtime(xtime(xtime(xtime(x))))))   \

// Main Process
#define RUN_ALL       0
#define RUN_TFM       1
#define RUN_LLT       2
#define RUN_SERVER    3
#define RUN_CONSOLE   4
#define RUN_TEST      5
#define RUN_KILL      6

#define FINDER_NORMAL   0
#define FINDER_NOMOD    1
#define FINDER_SKIP     2
#define FINDER_EULER    3

#define PROCESS_ID        3
#define PROCESS_ID_CLASS  0
#define PROCESS_ID_PID    1
#define PROCESS_ID_STATUS 2

#define PRIME_ID       3
#define PRIME_ID_USE   0
#define PRIME_ID_VALUE 1
#define PRIME_ID_PRIME 2

#define ALGO_KJ       0
#define ALGO_MULMOD   1

#define USE_NOFFT     0 
#define USE_FFT       1 

#define USE_NOFILE    0
#define USE_FILE      1 

#define USE_NOOMP     0 
#define USE_OMP       1 

#define USE_DIRECT    0 
#define USE_SOCKET    1 

#define LOG_NONE     0
#define LOG_ERROR    1
#define LOG_INFO     2
#define LOG_DEBUG    3
#define LOG_DUMP     4
#define LOG_FILE     5
#define LOG_NETWORK  6
#define LOG_VERBOSE  7

#define TFM_PROCESS_NAME "TFM"
#define LLT_PROCESS_NAME "LLT"

#define PROCESS_FILE_NAME "process"
#define STATUS_FILE_NAME "status"
#define LOG_FILE_NAME "mersenne"
#define LOCK_FILE_NAME "lock.pid"
#define DB_FILE_NAME "db.dat"
#define CANDIDATE_FILE_NAME "candidate.dat"
#define LAST_FILE_NAME "candidate.last"

#define STATUS_INIT "init"
#define STATUS_READY "ready"
#define STATUS_RUN "run"
#define STATUS_STOP "stop"
#define STATUS_KILL "kill"
#define STATUS_EXIT "exit"

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

// Options
#define USE_COMPRESSION
#define USE_ENCRYPTION

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
void readCommand(int argc,char ** argv);
char g_achOptionPath[MAX_SIZE];

uint64_t g_nSystemMinCandidate = MIN_CANDIDATE;
uint64_t g_nSystemMaxCandidate = MAX_CANDIDATE;
uint64_t g_nMinCandidate = MIN_CANDIDATE;
uint64_t g_nMaxCandidate = MAX_CANDIDATE;
uint64_t g_nMinPrime = MIN_PRIME;
uint64_t g_nMaxPrime = MAX_PRIME;

int g_nRunMode = RUN_TFM;
int g_nAlgorithm = ALGO_KJ;
int g_nUseFFT = USE_FFT;
int g_nUseOMP = USE_OMP;
int g_nUseSocket = USE_SOCKET;
int g_nNoFile = USE_FILE;
int g_nFinder = FINDER_NORMAL;
int g_nTFMProcessCount = 1;
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

// Big Integer
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

// Big Integer Library
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
void square(uint64_t * src,uint64_t * dst);
void square_fft(uint64_t * src,uint64_t * dst);
void mulop(uint64_t op, uint64_t * dst);
void leftshift(uint64_t * src);
void leftshiftop(uint64_t nBits, uint64_t * src);
void rightshift(uint64_t * src);
void rightshiftop(uint64_t nBits, uint64_t * src);
void PrintArray(uint64_t * src);
void PrintBits(uint64_t * src);
uint64_t getDigits(uint64_t * arr);
uint64_t getBits(uint64_t * arr);
uint64_t getBitsop(uint64_t value);
int compare(uint64_t * src, uint64_t * dst);
int compareop(uint64_t op, uint64_t * dst);
bool iszero();
void propa_carrier(uint64_t nIndex, uint64_t * dst);

// Zeta Function ( need to convert big number version )
void zeta(long double s, long double si, long double * r, long double * ri);
double _pi(int accuracy);
double _log10(double x);
double _log(double x);
double _cos(double x, double y);
double _sin(double x, double y, int iy);

// Main Process
void runAll();
void runTFM();
void runLLT();
void runServer();
void runConsole();
void runTest();

// Mersenne Process
void launchLLT();
void launchTFM();

// LLT Process
void setMersenneNumber(uint64_t p);
bool PrimalityTesting(uint64_t p);
void LLTmethod(uint64_t * src, uint64_t *dst);
void LLTmulmod(uint64_t * src, uint64_t *dst);

// TFM Process
void findMersenne();
void findMersenneOMP();
void findChunkNormal();
void findChunkNomod();
void findChunkSkip();
void findChunkEuler();
void setLastTFM(long p);
long getLastTFM();

// Management Candidate
void initCandidate();
bool isPrime(uint64_t num);
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
void waitTFMReady();
bool isTFMReady();
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
// Math Function
static const double
ivln10hi  = 4.34294481878168880939e-01, /* 0x3fdbcb7b, 0x15200000 */
ivln10lo  = 2.50829467116452752298e-11, /* 0x3dbb9438, 0xca9aadd5 */
log10_2hi = 3.01029995663611771306e-01, /* 0x3FD34413, 0x509F6000 */
log10_2lo = 3.69423907715893078616e-13, /* 0x3D59FEF3, 0x11F12B36 */
Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

static const double
ln2_hi = 6.93147180369123816490e-01,  /* 3fe62e42 fee00000 */
ln2_lo = 1.90821492927058770002e-10;  /* 3dea39ef 35793c76 */

static const double
C1  =  4.16666666666666019037e-02, /* 0x3FA55555, 0x5555554C */
C2  = -1.38888888888741095749e-03, /* 0xBF56C16C, 0x16C15177 */
C3  =  2.48015872894767294178e-05, /* 0x3EFA01A0, 0x19CB1590 */
C4  = -2.75573143513906633035e-07, /* 0xBE927E4F, 0x809C52AD */
C5  =  2.08757232129817482790e-09, /* 0x3E21EE9E, 0xBDB4B1C4 */
C6  = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

static const double
S1  = -1.66666666666666324348e-01, /* 0xBFC55555, 0x55555549 */
S2  =  8.33333333332248946124e-03, /* 0x3F811111, 0x1110F8A6 */
S3  = -1.98412698298579493134e-04, /* 0xBF2A01A0, 0x19C161D5 */
S4  =  2.75573137070700676789e-06, /* 0x3EC71DE3, 0x57B1FE7D */
S5  = -2.50507602534068634195e-08, /* 0xBE5AE5E6, 0x8A2B9CEB */
S6  =  1.58969099521155010221e-10; /* 0x3DE5D93A, 0x5ACFD57C */

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

// Process Entry Point
int main(int argc, char ** argv) {
  readCommand(argc, argv);
  logging(LOG_INFO, "Platform %s\n", PLATFORM_NAME);

  if(g_nRunMode == RUN_ALL) {
    runAll();
  } else if(g_nRunMode == RUN_TFM) {
    runTFM();
  } else if(g_nRunMode == RUN_LLT) {
    runLLT();
  } else if(g_nRunMode == RUN_SERVER) {
    runServer();
  } else if(g_nRunMode == RUN_CONSOLE) {
    runConsole();
  } else if(g_nRunMode == RUN_TEST) {
    runTest();
  }else {
    logging(LOG_INFO, "Unknown Run Mode [%d].\n", g_nRunMode);
  }
}

void runAll() {
  logging(LOG_DEBUG, "This is launcher process. pid=[%ld]\n", (long)getpid()); 
  logging(LOG_INFO, "All process START!\n");
  launchTFM();
  launchLLT();
}

void runTFM(){
  logging(LOG_INFO, "TFM Process START!\n");

  char achPID[MAX_SIZE];
  memset(achPID, 0x00, MAX_SIZE);
  sprintf(achPID, "%ld", (long)getpid());
  updateProcess(TFM_PROCESS_NAME, achPID, STATUS_INIT);
  initCandidate();
  updateProcess(TFM_PROCESS_NAME, achPID, STATUS_READY);
  waitTFMReady();

  logging(LOG_INFO, "runTFM : TFM Ready!\n");
  g_lstPrime = getPrimeList();
  g_lstCandidate = getCandidateList();
  g_pNextPrime = g_lstPrime;
  findMersenne();
  freeList(g_lstPrime);
  freeList(g_lstCandidate);
}

void runLLT(){
  logging(LOG_INFO, "LLT Process START!\n");

  char achPID[MAX_SIZE];
  memset(achPID, 0x00, MAX_SIZE);
  sprintf(achPID, "%ld", (long)getpid());
  updateProcess(LLT_PROCESS_NAME, achPID, STATUS_INIT);
  waitTFMReady();
  updateProcess(LLT_PROCESS_NAME, achPID, STATUS_READY);
  logging(LOG_INFO, "runLLT : TFM Ready!\n");
  readBits();

  if(g_nRunMode == RUN_ALL || g_nRunMode == RUN_LLT) {
    while(1) {
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
        readBits();
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

void launchTFM() {
  pid_t pid;

  g_nSystemMinCandidate = g_nMinCandidate;
  g_nSystemMaxCandidate = g_nMaxCandidate;

  for(int i = 1; i <= g_nTFMProcessCount; i++) {
    pid = fork();

    int nMinPrime = 0;
    int nMaxPrime = 0;
    nMinPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTFMProcessCount)*(i-1);
    if(i != g_nTFMProcessCount) {
      nMaxPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTFMProcessCount)*i - 1;
    } else {
      nMaxPrime = g_nMaxPrime;
    }
    int nMinCandidate = 0;
    int nMaxCandidate = 0;
    nMinCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nTFMProcessCount)*(i-1);
    if(i != g_nTFMProcessCount) {
      nMaxCandidate = g_nMinCandidate + ((g_nMaxCandidate - g_nMinCandidate)/g_nTFMProcessCount)*i - 1;
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
        logging(LOG_DEBUG, "This is child TFM process. pid=[%ld]\n", (long)getpid());
        g_nChildProcessID = pid;
        g_nMinPrime = nMinPrime;
        g_nMaxPrime = nMaxPrime;
        g_nMinCandidate = nMinCandidate;
        g_nMaxCandidate = nMaxCandidate;
        runTFM();
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
  if(g_nUseOMP == USE_OMP) {
    mul_omp(src1,src2,dst);
    return;
  }
  if(g_nUseFFT == USE_FFT) {
    mul_fft(src1,src2,dst);
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

  #pragma omp parallel for num_threads(10)
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
  mul(src, src, dst);
}

void mul_fft(uint64_t * src1,uint64_t * src2, uint64_t * dst) {
  int size = getDigits(src1) + getDigits(src2) +1;

  memset(x, 0x00, sizeof(double)*MAX_INTEGER);
  memset(xi, 0x00, sizeof(double)*MAX_INTEGER);
  memset(rx, 0x00, sizeof(double)*MAX_INTEGER);
  memset(rxi, 0x00, sizeof(double)*MAX_INTEGER);
  memset(y, 0x00, sizeof(double)*MAX_INTEGER);
  memset(yi, 0x00, sizeof(double)*MAX_INTEGER);
  memset(ry, 0x00, sizeof(double)*MAX_INTEGER);
  memset(ryi, 0x00, sizeof(double)*MAX_INTEGER);

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

  memset(x, 0x00, sizeof(double)*MAX_INTEGER);
  memset(xi, 0x00, sizeof(double)*MAX_INTEGER);

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
    uint64_t value = dst[i] + (uint64_t)llround(rx[i]);
    dst[i] = value % bits_max;
    dst[i+1] = value / bits_max;
  }
}

void square_fft(uint64_t * src, uint64_t * dst) {
  int size = getDigits(src)*2 + 1;
  memset(x, 0x00, sizeof(double)*MAX_INTEGER);
  memset(xi, 0x00, sizeof(double)*MAX_INTEGER);
  memset(rx, 0x00, sizeof(double)*MAX_INTEGER);
  memset(rxi, 0x00, sizeof(double)*MAX_INTEGER);

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

  memset(x, 0x00, sizeof(double)*MAX_INTEGER);
  memset(xi, 0x00, sizeof(double)*MAX_INTEGER);

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

  memset(dst, 0x00, sizeof(uint64_t)*MAX_INTEGER);

  // compute carry
  for(int i = 0; i < size; i++) {
    uint64_t value = dst[i] + (uint64_t)llround(rx[i]);
    dst[i] = value % bits_max;
    dst[i+1] = value / bits_max;
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
      logging(LOG_INFO, "TFM Last Prime : [%ld]\n", getLastTFM());
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

void test() {
  initInteger();
  PrintArray(bigS);

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
  logging(LOG_DEBUG, "initCandidate Begin\n");
  if(g_nNoFile == USE_FILE && access(CANDIDATE_FILE_NAME, F_OK) == 0) {
    logging(LOG_DEBUG, "initCandidate From File Min[%llu]~Max[%llu]\n", g_nMinCandidate, g_nMaxCandidate);
    readBits();
  } else {
    logging(LOG_DEBUG, "initCandidate From Process Min[%llu]~Max[%llu]\n", g_nMinCandidate, g_nMaxCandidate);
    for(uint64_t i = g_nMinCandidate; i < g_nMaxCandidate; i+=2) {
      if(isPrime(i)) {
        setNumberCache(i, MERSENNE_CANDIDATE);
      }
    }
    updateBits();
  }
  logging(LOG_DEBUG, "initCandidate End\n");
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
    writeBits();
    oldexponent = exponent;
    pNextCandidate=pNextCandidate->next;
  }
}

void findChunkNomod() {
  uint64_t chunk[g_nChunk][PRIME_ID];
  memset(chunk, 0x00, PRIME_ID*g_nChunk*sizeof(uint64_t));
  int count = getNextChunk(chunk);

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
    writeBits();
    oldexponent = exponent;
    pNextCandidate=pNextCandidate->next;
  }
}

void findChunkSkip() {
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
            chunk[j][PRIME_ID_USE] = 0;
          }
          comparer <<= 1;
        }
      }
    }
    writeBits();
    oldexponent=pNextCandidate->data;
    pNextCandidate=pNextCandidate->next;
  }
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
    writeBits();
    oldexponent = pNextCandidate->data;
    pNextCandidate=pNextCandidate->next;
  }
}

struct Node* getPrimeList() {
  g_nMinPrime = getLastTFM();

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
  struct Node* candidate = (struct Node*)malloc(sizeof(struct Node *));
  candidate->data = g_nMinCandidate;
  candidate->next = NULL;

  struct Node* pNextCandidate = candidate;
  for(int i = g_nMinCandidate + 2; i < g_nMaxCandidate; i+=2) {
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
  if(num == 2) return true;
  if(num % 2 == 0) return false;
  for(uint64_t i = 3; i*i<=num; i+=2) {
    if(num % i == 0) return false;
  }
  return true;
}

int test_get_primes() {
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
  exit(0);

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
        printf("[i=%llu] maxDiff : [%llu] primeCount : [%llu]\n", i, maxDiff, primeCount);
        for(int j = 0; j < MAX_DIFF; j++) {
          printf("[%llu]", arr[j]);
          //reset
          //arr[j] = 0;
        }
        printf("\n");
      }
      ppos = i;
    }
  }
  printf("[i=1000000000] maxDiff : [%llu] primeCount : [%llu]\n", maxDiff, primeCount);
  for(int i = 1; i < MAX_DIFF; i++) {
    printf("[%llu]", arr[i]);
  }
}

int test_compress() {
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

  printf("Count of jobs on the queue : %d\n", pq->usedSize);
  while(!isEmpty(pq)) {
    dequeue(pq, &result);
    printNode(&result);
  }

  memset(&codeTable, 0, sizeof(HFCode)*MAX_CHAR);
  encode(&tree, (unsigned char*)source, &encoded, codeTable);

  printf("Original Size: %d-bit, encoded Size: %d-bit\n",
    (int)((strlen(source) + 1) * sizeof(char) * 8), encoded.size);

  decoded = (char*)malloc(sizeof(char) * (strlen(source) + 1));
  decode(tree, &encoded, (unsigned char*)decoded);

  printf("Original: %s\n", source);
  printf("encoded: ");
  printBinary(&encoded);
  printf("\ndecoded: %s\n", decoded);

  free(decoded);
  destroyTree(tree);

  return 0;
}

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
    g_nTFMProcessCount = 1;
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
    g_nTFMProcessCount = configInt("TFMProcess");
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
  memcpy(cache, candidate, g_nSystemMaxCandidate/64 + 1);
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
    for(int i = 0; i < (g_nMaxCandidate/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
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
  memcpy(cache, candidate, g_nSystemMaxCandidate/64 + 1);
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
    for(int i = 0; i < (g_nMaxCandidate/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
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
  memcpy(cache, candidate, g_nSystemMaxCandidate/64 + 1);
  readBits();
  int shift = g_nMinCandidate%64;
  candidate[g_nMinCandidate/64] = ((candidate[g_nMinCandidate/64] >> shift) << shift) + ((cache[g_nMinCandidate/64] << shift) >> shift);
  for(int i = g_nMinCandidate/64 + 1; i < g_nMaxCandidate/64; i++) {
    candidate[i] = cache[i];
  }
  candidate[g_nMaxCandidate/64 + 1] = ((cache[g_nMaxCandidate/64 + 1] >> shift) << shift) + ((candidate[g_nMaxCandidate/64 + 1] << shift) >> shift);
  writeBits();
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

void getPrime(uint64_t maxPrime, uint64_t * primes) {
  uncompressPrimes(maxPrime, primes);
}

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
        writePrime(maxPrime, primes);
        writeNoise(maxPrime, noises);
        memset(primes, 0x00, sizeof(uint8_t)*MAX_PRIMES);
        memset(noises, 0x00, sizeof(uint8_t)*MAX_NOISES);
        nPrimeIndex = 0;
      }
    }
  }
  compressPrimes(maxPrime);
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

#ifdef USE_COMPRESSION
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

int test_sign () {
    uint64_t e = E_VALUE;
    uint64_t phi = 0;

    uint64_t d, n, p, q, h, m, qInv, m1m2;
    uint64_t c, dP, dQ, m1, m2;

    printf("p and q must be prime numbers, e must be coprime to (p - 1)*(q - 1)\n");
    setprimes(e, &p, &q, &n, &phi);

    printf("p: %llu q: %llu n: %llu phi: %llu\n", p, q, n, phi);

    d = findD(e,phi);
    printf("find d: %llu\n", d);

    printf("Public Key:  (n,e) = (%llu, %llu)\n", n, e);
    printf("Private Key: (n,d) = (%llu, %llu)\n", n, d);
    printf("p = %llu, q = %llu, phi = %llu\n", p, q, phi);
    printf("Key generation success!\n\n");

    char message[MAX_SIZE];
    memset(message, 0x00, MAX_SIZE);
    strcpy(message, "This is a plain Text!");
    printf("original message [%s]\n", message);
 
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
    printf("decrypted [%s]\n", decrypted);
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
  for(int i = g_nMinCandidate; i < g_nMaxCandidate; i+=2) {
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
  return (g_nTFMProcessCount == 1);
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

void waitTFMReady() {
  do {
#ifdef _WIN32
    Sleep(1000);
#else
    sleep(1);
#endif
  } while(!isTFMReady());
}

bool isTFMReady() {
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
    if(strcmp(TFM_PROCESS_NAME, processStatus[i].class) == 0 && strcmp(STATUS_READY, processStatus[i].status) != 0) {
      lock_leave();
      return false;
    }
    if(strcmp(TFM_PROCESS_NAME, processStatus[i].class) == 0 && strcmp(STATUS_READY, processStatus[i].status) == 0) {
      count++;
    } 
  }
  if(g_nTFMProcessCount == count) {
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
}

void setLastTFM(long p) {
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

long getLastTFM() {
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

// Zeta Function
void zeta(long double s, long double si, long double * r, long double * ri) {
  long double a_arr[MAXNUM + 1];
  long double a_arri[MAXNUM + 1];
  long double half = 0.5;
  long double halfi = 0.0;
  long double one = 1.0;
  long double onei = 0.0;
  long double two = 2.0;
  long double twoi = 0.0;
  long double rev = -1.0;
  long double revi = 0.0;
  long double sum = 0.0;
  long double sumi = 0.0;
  long double prev = 1.0e+20;
  long double previ = 0.0;

  memset(a_arr, 0x00, sizeof(long double)*(MAXNUM + 1));
  memset(a_arri, 0x00, sizeof(long double)*(MAXNUM + 1)); 

  // initialize with a_0 = 0.5 / (1 - 2^(1-s))
  a_arr[0] = half / (one - pow(two, (one - s)));
  a_arri[0] = halfi / (onei - pow(twoi, (onei - s)));
  
  sum += a_arr[0]; 
  sumi += a_arri[0]; 

  for(int n = 1; n <= MAXNUM; n++) {
    long double nCplx = n;
    long double nCplxi = 0.0;

    for(int k = 0; k < n; k++) {
      // complex index
      long double kCplx = k;
      long double kCplxi = 0.0;
   
      a_arr[k] *= half * (nCplx / (nCplx - kCplx));
      a_arri[k] *= halfi * (nCplxi / (nCplxi - kCplxi));
      sum += a_arr[k];
      sumi += a_arri[k];
    }

    a_arr[n] = (rev * a_arr[n-1]* pow((nCplx / (nCplx + one)), s) / nCplx);
    a_arri[n] = (revi * a_arri[n-1]* pow((nCplxi / (nCplxi + onei)), si) / nCplxi);
    sum += a_arr[n];
    sumi += a_arr[n];

   
    if( abs(prev - sum) < LOWER_THRESHOLD && abs(previ - sumi) < LOWER_THRESHOLD) {
    // If the differences is less than or equal to the threshold value, it is considered to be convergent and the calculation is terminated.
       break;
     } 
     if( abs(sum) > UPPER_BOUND && abs(sumi) > UPPER_BOUND) {
       // doesn't work for large values, so it gets terminated when it exceeds UPPER_BOUND
       break;
     }
     prev = sum;
     previ = sumi;
   }
   *r = sum;
   *ri = sumi;
}

double _pi(int accuracy){
     double result = 1;
     int a = 2;
     int b = 1;

     for(int i = 0;i < accuracy; i ++){
          result = a/b * result;
          if(a < b){
               a = a + 2;
          }
          else if(b < a){
               b = b + 2;
          }
     }

     return result * 2;
}

double _log10(double x) {
	union {double f; uint64_t i;} u = {x};
	double_t hfsq,f,s,z,R,w,t1,t2,dk,y,hi,lo,val_hi,val_lo;
	uint32_t hx;
	int k;
	hx = u.i>>32;
	k = 0;
	if (hx < 0x00100000 || hx>>31) {
		if (u.i<<1 == 0)
			return -1/(x*x);  /* log(+-0)=-inf */
		if (hx>>31)
			return (x-x)/0.0; /* log(-#) = NaN */
		/* subnormal number, scale x up */
		k -= 54;
		x *= 0x1p54;
		u.f = x;
		hx = u.i>>32;
	} else if (hx >= 0x7ff00000) {
		return x;
	} else if (hx == 0x3ff00000 && u.i<<32 == 0)
		return 0;
	/* reduce x into [sqrt(2)/2, sqrt(2)] */
	hx += 0x3ff00000 - 0x3fe6a09e;
	k += (int)(hx>>20) - 0x3ff;
	hx = (hx&0x000fffff) + 0x3fe6a09e;
	u.i = (uint64_t)hx<<32 | (u.i&0xffffffff);
	x = u.f;
	f = x - 1.0;
	hfsq = 0.5*f*f;
	s = f/(2.0+f);
	z = s*s;
	w = z*z;
	t1 = w*(Lg2+w*(Lg4+w*Lg6));
	t2 = z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
	R = t2 + t1;
	/* See log2.c for details. */
	/* hi+lo = f - hfsq + s*(hfsq+R) ~ log(1+f) */
	hi = f - hfsq;
	u.f = hi;
	u.i &= (uint64_t)-1<<32;
	hi = u.f;
	lo = f - hi - hfsq + s*(hfsq+R);
	/* val_hi+val_lo ~ log10(1+f) + k*log10(2) */
	val_hi = hi*ivln10hi;
	dk = k;
	y = dk*log10_2hi;
	val_lo = dk*log10_2lo + (lo+hi)*ivln10lo + lo*ivln10hi;
	/*
	 * Extra precision in for adding y is not strictly needed
	 * since there is no very large cancellation near x = sqrt(2) or
	 * x = 1/sqrt(2), but we do it anyway since it costs little on CPUs
	 * with some parallelism and it reduces the error for many args.
	 */
	w = y + val_hi;
	val_lo += (y - w) + val_hi;
	val_hi = w;
	return val_lo + val_hi;
}

double _log(double x)
{
	union {double f; uint64_t i;} u = {x};
	double_t hfsq,f,s,z,R,w,t1,t2,dk;
	uint32_t hx;
	int k;
	hx = u.i>>32;
	k = 0;
	if (hx < 0x00100000 || hx>>31) {
		if (u.i<<1 == 0)
			return -1/(x*x);  /* log(+-0)=-inf */
		if (hx>>31)
			return (x-x)/0.0; /* log(-#) = NaN */
		/* subnormal number, scale x up */
		k -= 54;
		x *= 0x1p54;
		u.f = x;
		hx = u.i>>32;
	} else if (hx >= 0x7ff00000) {
		return x;
	} else if (hx == 0x3ff00000 && u.i<<32 == 0)
		return 0;
	/* reduce x into [sqrt(2)/2, sqrt(2)] */
	hx += 0x3ff00000 - 0x3fe6a09e;
	k += (int)(hx>>20) - 0x3ff;
	hx = (hx&0x000fffff) + 0x3fe6a09e;
	u.i = (uint64_t)hx<<32 | (u.i&0xffffffff);
	x = u.f;
	f = x - 1.0;
	hfsq = 0.5*f*f;
	s = f/(2.0+f);
	z = s*s;
	w = z*z;
	t1 = w*(Lg2+w*(Lg4+w*Lg6));
	t2 = z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
	R = t2 + t1;
	dk = k;
	return s*(hfsq+R) + dk*ln2_lo - hfsq + f + dk*ln2_hi;
}

double _cos(double x, double y)
{
	double_t hz,z,r,w;
	z  = x*x;
	w  = z*z;
	r  = z*(C1+z*(C2+z*C3)) + w*w*(C4+z*(C5+z*C6));
	hz = 0.5*z;
	w  = 1.0-hz;
	return w + (((1.0-w)-hz) + (z*r-x*y));
}

double _sin(double x, double y, int iy)
{
	double_t z,r,v,w;
	z = x*x;
	w = z*z;
	r = S2 + z*(S3 + z*S4) + z*w*(S5 + z*S6);
	v = z*x;
	if (iy == 0)
		return x + v*(S1 + z*r);
	else
		return x - ((z*(0.5*y - v*r) - y) - v*S1);
}
