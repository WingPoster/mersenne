// compile : gcc mersenne.c -o mersenne
// run : ./mersenne option.config

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <arpa/inet.h> // socket
#include <sys/types.h> // socket
#include <sys/socket.h> // socket

#ifdef WIN32
#include <io.h>
#define F_OK 0
#define access _access
#else
#include <unistd.h>
#endif

#define RUN_ALL       0
#define RUN_TFM       1
#define RUN_LLT       2
#define RUN_SERVER    3
#define RUN_CONSOLE   4
#define RUN_TEST      5

#define ALGO_KJ       0
#define ALGO_MULMOD   1

#define USE_NOFFT     0 
#define USE_FFT       1 

#define USE_DIRECT    0 
#define USE_SOCKET    1 

#define MAX_SIZE      1024
#define MAX_BUFF_SIZE 8092
#define MAX_INTEGER   10000000

#define MIN_MERSENNE  2
#define MAX_MERSENNE  320000000 // 3.2 billion
#define MIN_PRIME     3
#define MAX_PRIME     1000000000

#define MERSENNE_NUMBER 0
#define MERSENNE_PRIME  1

#define FULL_LIST  0
#define EULER_LIST 1

#define LOG_ERROR 0
#define LOG_INFO  1
#define LOG_DEBUG 2
#define LOG_DUMP  3

#define KNOWN_MERSENNE 51

#define STATUS_FILE_NAME "status"
#define LOCK_FILE_NAME "lock.pid"
#define DB_FILE_NAME "db.dat"
#define MERSENNE_FILE_NAME "mersenne.dat"

#define PI 3.14159265358979

#define STATUS_INIT "init"
#define STATUS_READY "ready"
#define STATUS_RUN "run"
#define STATUS_STOP "stop"
#define STATUS_EXIT "exit"

int g_nAlgorithmType = FULL_LIST;

uint64_t mersenne[MAX_MERSENNE/64 + 1];

uint64_t g_nMinMersenne = MIN_MERSENNE;
uint64_t g_nMaxMersenne = MAX_MERSENNE;
uint64_t g_nMinPrime = MIN_PRIME;
uint64_t g_nMaxPrime = MAX_PRIME;
uint64_t g_nStartPrime = 2;
uint64_t g_nEndPrime = MAX_PRIME;

// sequence A000043 in the OEIS
int knownMersenne[] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091, 756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951, 30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933};

int g_nRunMode = RUN_TFM;
int g_nAlgorithm = ALGO_KJ;
int g_nUseFFT = USE_FFT;
int g_nUseSocket = USE_SOCKET;
int g_nNoFile = 1;
int g_nTFMProcessCount = 1;
int g_nLLTProcessCount = 1;
int g_nPort = 8080;
int g_nLogLevel = LOG_DEBUG;

char g_achOptionPath[MAX_SIZE];
char g_achAddress[MAX_SIZE];
int g_nChildProcessID = 0;
int g_nIsMain = 0;

uint64_t g_nMersenne = 3;
uint64_t bits_max = 4294967296;
uint64_t p = 0;
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

void setMersenne(uint64_t p);
int  PrimalityTesting(uint64_t p);
void test();
void server();
void agent(int sock);
void console();
void logging(int level, const char *format, ...);
void printUsage();
void readCommand(int argc,char ** argv);
void readBits();
void writeBits();
void readDB();
void writeDB();
void setBit(uint64_t nIndex, int bit);
int setBitOnly(uint64_t nIndex, int bit);
int getBit(uint64_t nIndex);
int getCandidateCount();
int readConfigInteger(char *pKey);
void readConfig(char *pKey, char *pValue);

void add(uint64_t * src, uint64_t * dst);
void addop(uint64_t op, uint64_t * dst);
void set(uint64_t * dst);
void subtract(uint64_t * src, uint64_t * dst);
void subtractop(uint64_t op, uint64_t * dst);
void copy(uint64_t * src, uint64_t * dst);
void copyBits(uint64_t * src, uint64_t * dst, uint64_t nStart, uint64_t nEnd);
void mod(uint64_t * src, uint64_t * dst);
void mul(uint64_t * src1, uint64_t * src2, uint64_t * dst);
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
void PrintMersenneCount();

uint64_t getDigits(uint64_t * arr);
uint64_t getBits(uint64_t * arr);
int compare(uint64_t * src, uint64_t * dst);
int compareop(uint64_t op, uint64_t * dst);
int iszero();
void propa_carrier(uint64_t nIndex, uint64_t * dst);

void LLTmethod(uint64_t * src, uint64_t *dst);
void LLTmulmod(uint64_t * src, uint64_t *dst);

void initMersenne();
void findMersenne();
int isPrime(uint64_t num);
int isMersenne(uint64_t num);
int getKnownMersenne(int number);

void runAll();
void runTFM();
void runLLT();
void runServer();
void runConsole();
void runTest();

void launchLLT();
void launchTFM();

time_t tStart = 0;
void error_handling(char *message);

void lock_enter();
void lock_leave();

void getStatus(char * pchStatus);
void setStatus(char * pchStatus);

int main(int argc, char ** argv) {
  readCommand(argc, argv);
  tStart = time(NULL);

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
  if(g_nIsMain) {
    logging(LOG_INFO, "Processing Time : %ld\n", time(NULL) - tStart);
  }
}

void runAll() {
  logging(LOG_INFO, "All process START!\n");
  logging(LOG_INFO, "TFM process START!\n");
  launchTFM();
  launchLLT();
}

void runTFM(){
  logging(LOG_INFO, "TFM Mode START!\n");
  setStatus(STATUS_INIT);
  initMersenne();
  setStatus(STATUS_READY);
  findMersenne();
}

void runLLT(){
  logging(LOG_INFO, "LLT Mode START!\n");

  char achStatus[MAX_SIZE];
  do {
    memset(achStatus, 0x00, MAX_SIZE);
    getStatus(achStatus);
    sleep(1);
  } while(strcmp(achStatus, STATUS_READY) != 0);
  readBits();

  if(g_nRunMode == RUN_ALL) {
    for(int i = 0; i < g_nMaxMersenne; i++) {
      if(getBit(i) == MERSENNE_NUMBER || getKnownMersenne(i) == 1) {
        continue;
      }
      g_nMersenne = i;
      if(PrimalityTesting(i)) {
        logging(LOG_INFO, "2^%d - 1 is prime!!\n", i);
      } else {
        setBit(i, MERSENNE_NUMBER);
      }
    }
  } else if(g_nRunMode == RUN_LLT) {
    initInteger();
    for(int i = g_nStartPrime; i <= g_nEndPrime; i+=2) {
      if(getBit(i) == MERSENNE_NUMBER) {
        continue;
      }
      g_nMersenne = i;
      if(PrimalityTesting(i)) {
        logging(LOG_INFO, "2^%d - 1 is prime!!\n", i);
      } else {
        setBit(i, MERSENNE_NUMBER);
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
  for(int i = 1; i <= g_nLLTProcessCount; i++) {
    pid = fork();
    int nStartPrime = 0;
    int nEndPrime = 0;
    nStartPrime = g_nStartPrime + ((g_nEndPrime - g_nStartPrime)/g_nTFMProcessCount)*i;
    if(i != g_nTFMProcessCount -1) {
      nEndPrime = g_nStartPrime + ((g_nEndPrime - g_nStartPrime)/g_nTFMProcessCount)*(i+1) - 1;
    } else {
      nEndPrime = g_nEndPrime;
    }
    switch(pid) {
      case -1:
        logging(LOG_INFO, "fork failed, reason=failed create child process.\n");
        return;
      case 0:
        g_nIsMain = 1;
        logging(LOG_DEBUG, "This is parent process. pid=[%ld]\n", (long)getpid()); 
        break;
      default:
        logging(LOG_DEBUG, "This is child process. pid=[%ld]\n", (long)getpid());
        g_nChildProcessID = pid;
        g_nStartPrime = nStartPrime;
        g_nEndPrime = nEndPrime;
        runLLT();
        logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
        break;
    }
  }
}

void launchTFM() {
  pid_t pid;
  for(int i = 0; i < g_nTFMProcessCount; i++) {
    pid = fork();
    int nMinPrime = 0;
    int nMaxPrime = 0;
    nMinPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTFMProcessCount)*i;
    if(i != g_nTFMProcessCount -1) {
      nMaxPrime = g_nMinPrime + ((g_nMaxPrime - g_nMinPrime)/g_nTFMProcessCount)*(i+1) - 1;
    } else {
      nMaxPrime = g_nMaxPrime;
    }
    switch(pid) {
      case -1:
        logging(LOG_DEBUG, "fork failed, reason=failed create child process.\n");
        return;
      case 0:
        g_nIsMain = 1;
        logging(LOG_DEBUG, "This is parent process. pid=[%ld]\n", (long)getpid()); 
        break;
      default:
        logging(LOG_DEBUG, "This is child process. pid=[%ld]\n", (long)getpid());
        g_nChildProcessID = pid;
        g_nMinPrime = nMinPrime;
        g_nMaxPrime = nMaxPrime;
        runTFM();
        logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
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
int PrimalityTesting(uint64_t p) {
  time_t tStepStart = 0;
  time_t tStepEnd = 0;
  set(bigS);
  addop(4, bigS);

  setMersenne(p);
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
      return 1;
    } else if(i == p-2 && !iszero(bigS)) {
      return 0;
    }
  }
  return 0;
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
  uint64_t nPartition = g_nMersenne + (4 - g_nMersenne%4);
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
  copyBits(dst, bigA, nSubPart*2,  nSubPart*4 -1);
  set(bigL);
  square(bigA, bigL);
  if(nSubPart*4 > g_nMersenne) {
    leftshiftop(nSubPart*4 - g_nMersenne, bigL);
  }
  add(bigL, bigBuff);

  // A1*B1
  set(bigA);
  set(bigB);
  copyBits(dst, bigA, nSubPart*3,  nSubPart*4 -1);
  copyBits(dst, bigB, nSubPart, nSubPart*2 -1);
  set(bigL);
  mul(bigA, bigB, bigL);
  if(nSubPart*4 + 1 > g_nMersenne) {
    leftshiftop(nSubPart*4 - g_nMersenne + 1, bigL);
  } else {
    leftshiftop(nSubPart*4 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A1*B2
  set(bigA);
  set(bigB);
  copyBits(dst, bigA, nSubPart*3,  nSubPart*4 -1);
  copyBits(dst, bigB, 0,  nSubPart -1);
  set(bigL);
  mul(bigA, bigB, bigL);
  if(nSubPart*3 + 1 > g_nMersenne) {
    leftshiftop(nSubPart*3 - g_nMersenne + 1, bigL);
  } else {
    leftshiftop(nSubPart*3 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A2*B1
  set(bigA);
  set(bigB);
  copyBits(dst, bigA, nSubPart*2,  nSubPart*3 -1);
  copyBits(dst, bigB, nSubPart, nSubPart*2 -1);
  set(bigL);
  mul(bigA, bigB, bigL);
  if(nSubPart*3 + 1 > g_nMersenne) {
    leftshiftop(nSubPart*3 - g_nMersenne + 1, bigL);
  } else {
    leftshiftop(nSubPart*3 + 1, bigL);
  }
  add(bigL, bigBuff);

  // A2*B2
  set(bigA);
  set(bigB);
  copyBits(dst, bigA, nSubPart*2,  nSubPart*3 -1);
  copyBits(dst, bigB, 0,  nSubPart -1);
  set(bigL);
  mul(bigA, bigB, bigL);
  if(nSubPart*2 + 1 > g_nMersenne) {
    leftshiftop(nSubPart*2 - g_nMersenne + 1, bigL);
  } else {
    leftshiftop(nSubPart*2 + 1, bigL);
  }
  add(bigL, bigBuff);

  // B^2
  copyBits(dst, bigA, 0,  nSubPart*2 -1);
  set(bigL);
  square(bigA, bigL);
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

// The result must be unsigned integer.
void subtract(uint64_t * src, uint64_t * dst) {
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

void setMersenne(uint64_t p) {
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

int iszero(uint64_t * arr) {
  for(int i = 0; i < MAX_INTEGER; i++) {
    if(arr[i] != 0) {
      return 0;
    }
  }
  return 1;
}

void error_handling(char *message) {
  logging(LOG_ERROR,"socket error\n");
  exit(1);
}

void server() {
  int serv_sock;
  int clnt_sock;

  char message[MAX_SIZE];
  int str_len;

  struct sockaddr_in serv_addr;
  struct sockaddr_in clnt_addr;
  int clnt_addr_size;

  serv_sock = socket(PF_INET, SOCK_STREAM, 0);

  if(serv_sock == -1)
    error_handling("socket() error");
  memset(&serv_addr, 0, sizeof(serv_addr));

  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr=htonl(INADDR_ANY);
  serv_addr.sin_port = htons(g_nPort);

  if(bind(serv_sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr))==-1)
    error_handling("bind() error");
  if(listen(serv_sock, 5) == -1)
    error_handling("listen() error");
  clnt_addr_size = sizeof(clnt_addr);

  readDB();
  while(1) {
    clnt_sock = accept(serv_sock, (struct sockaddr*)&clnt_addr,(socklen_t *) &clnt_addr_size);
    if(clnt_sock == -1)
      error_handling("accept() error");
    agent(clnt_sock);
  }
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
      g_nIsMain = 1;
      logging(LOG_DEBUG, "This is parent process. pid=[%ld]\n", (long)getpid());
      break;
    default:
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
        write(sock, "HTTP1.1 200 OK\n\n", strlen("HTTP1.1 200 OK\n\n"));
      }
      close(sock);
      logging(LOG_DEBUG, "Child process id=[%d]\n", pid);
      break;
  }
}

void console() {
  if(g_nUseSocket == USE_DIRECT) {
    readBits();
    int nCount = 0;
    int nMinMersenne = 0;
    for(int i = 0; i < g_nMaxMersenne; i++) {
      if(getBit(i) == 1 && getKnownMersenne(i) == 0) {
        nCount++;
        if(nMinMersenne == 0) {
          nMinMersenne = i;
        }
      }
    }
    logging(LOG_INFO, "Candidate Count : [%d] (M%llu ~ M%llu)\n", nCount, g_nMinMersenne, g_nMaxMersenne);
    logging(LOG_INFO, "Min Mersenne Candidate : [%d]\n", nMinMersenne);

    for(int i = 0; i < 1000; i++) {
      if(getBit(i) == 1 && getKnownMersenne(i) == 0) {
        logging(LOG_DEBUG, "x"); 
      } else {
        logging(LOG_DEBUG, " ");
      }
      if(i%100 == 0) {
        logging(LOG_DEBUG, "\n");
      }
    }
  } else if(g_nUseSocket == USE_SOCKET) {
    int sock;
    struct sockaddr_in serv_addr;
    char message[MAX_SIZE];
    int str_len;

    sock = socket(PF_INET, SOCK_STREAM, 0);
    if(sock == -1)
      error_handling("socket() error");

    memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr(g_achAddress);
    serv_addr.sin_port = htons(g_nPort);

    if(connect(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) == -1)
      error_handling("connect() error");

    while(1) {
      fgets(message, MAX_SIZE, stdin);
      write(sock, message, strlen(message));

      str_len = read(sock, message, MAX_SIZE-1);
      message[str_len] = 0;
      logging(LOG_DUMP, "%s\n", message);
    }
    close(sock);
  } else {
    logging(LOG_ERROR, "Unknown UserSocket Option [%d]\n", g_nUseSocket);
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
  g_nMersenne = 107;
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

  if(g_nLogLevel > level) {
    printf("%s", buf);
  }
}

void initMersenne() {
  logging(LOG_DEBUG, "initMersenne Begin\n");
  if(g_nNoFile == 0 && access(MERSENNE_FILE_NAME, F_OK) == 0) {
    logging(LOG_DEBUG, "initMersenne From File\n");
    readBits();
  } else {
    logging(LOG_DEBUG, "initMersenne From Process[%llu]~[%llu]\n", g_nMinMersenne, g_nMaxMersenne);
    for(uint64_t i = g_nMinMersenne; i < g_nMaxMersenne; i++) {
      if(isPrime(i)) {
        setBitOnly(i, MERSENNE_PRIME);
      }
    }
    writeBits();
  }
  logging(LOG_DEBUG, "initMersenne End\n");
}

void findMersenne() {
  int nPrime = 0;

  logging(LOG_DEBUG, "MinPrime[%llu] ~ MaxPrime[%llu] Start\n", g_nMinPrime, g_nMaxPrime);
  // Except p=2 and check only odds
  // if Mersenne Prime is 2^p-1, p is 8k+1/-1.
  for(uint64_t p=g_nMinPrime; p < g_nMaxPrime; p+=2) {
    if(getBit(p) == MERSENNE_NUMBER || getKnownMersenne(p) == 1) {
      continue;
    }
    int remain = p%8;
    if(((remain == 1 || remain == 7) || (p == 2 || p == 3 || p == 5))) {
      if(g_nAlgorithmType == FULL_LIST) {
        // First of all, prime check performed before check mersenne.
        if(isPrime(p)) {
          nPrime++;
          isMersenne(p);
        }
      } else if(g_nAlgorithmType == EULER_LIST) {
        if(isPrime(p)) {
          nPrime++;
          for(int i=2; i<=2; i+=2) {
            // p > 3
            if((p-1)/i <= 3) break;
            if((p-1)/i > g_nMaxMersenne) continue;
            // k < 100, k not a multiple of three.
            if( i%3 ==0 ) continue;

            // if p ≡ 3 (mod 4) and p > 3, then the prime 2p+1 divides the Mersenne number Mp.
            if(i==2 && ((p-1)/i)%4 == 3) {
              setBit((p -1)/i, MERSENNE_NUMBER);
              continue;
            }
          }
        }
      } else {
        logging(LOG_ERROR, "Unknown Algorithm Type [%d]\n", g_nAlgorithmType);
        exit(0);
      }
    }

    PrintMersenneCount();
  }
  logging(LOG_DEBUG, "MinPrime[%llu] ~ MaxPrime[%llu] End\n", g_nMinPrime, g_nMaxPrime);
}

void PrintMersenneCount() {
  static int oldcount = MAX_MERSENNE;
  int candidateCount = getCandidateCount();
  if(oldcount > candidateCount) {
    logging(LOG_DEBUG, "Mersenne Candidate : %d\n", candidateCount);
    oldcount = candidateCount;
  }
}

int isPrime(uint64_t num) {
  if(num == 2) return 1;
  if(num % 2 == 0) return 0;
  for(uint64_t i = 3; i*i<=num; i+=2) {
    if(num % i == 0) return 0;
  }
  return 1;
}

// TFM Mersenne Prime
int isMersenne(uint64_t num) {
  uint64_t value = 1;
  for(uint64_t i=1; i < g_nMaxMersenne; i++) {
    value = (value+value) % num;
    //value = (value << 1) - num;
    if(value == 1 && i < 30 &&  pow(2, i) -1 == num) {
      break;
    } else {
      if(value == 1) {
        setBit(i, MERSENNE_NUMBER);
        return 0;
      }
    }
  }
  return 1;
}

void readCommand(int argc, char ** argv) {
  // if app was crashed, we clear a previous lock.
  lock_leave();
  memset(g_achOptionPath, 0x00, MAX_SIZE);

  if(argv[1] == NULL) {
    g_nMinMersenne = 2;
    g_nMaxMersenne = 100000000;
    g_nMinPrime    = 3;
    g_nMaxPrime    = 1000000000;
    g_nAlgorithmType = FULL_LIST;
    g_nStartPrime = 3;
    g_nEndPrime = 1000000000;
    g_nRunMode = RUN_ALL;
    g_nAlgorithm = ALGO_KJ;
    g_nUseFFT = USE_NOFFT;
    g_nUseSocket = USE_DIRECT;
    g_nNoFile = 0;
    g_nTFMProcessCount = 1;
    g_nLLTProcessCount = 1;
    g_nLogLevel = LOG_DEBUG;
  } else {
    strcpy(g_achOptionPath, argv[1]);

    g_nMinMersenne = readConfigInteger("MinMersenne");
    g_nMaxMersenne = readConfigInteger("MaxMersenne");
    g_nMinPrime = readConfigInteger("MinPrime");
    g_nMaxPrime = readConfigInteger("MaxPrime");
    g_nAlgorithmType = readConfigInteger("AlgorithmType");
    g_nStartPrime = readConfigInteger("StartPrime");
    g_nEndPrime = readConfigInteger("EndPrime");
    g_nRunMode = readConfigInteger("RunMode");
    g_nAlgorithm = readConfigInteger("Algorithm");
    g_nUseFFT = readConfigInteger("UseFFT");
    g_nUseSocket = readConfigInteger("UseSocket");
    g_nNoFile = readConfigInteger("NoFile");
    g_nTFMProcessCount = readConfigInteger("TFMProcess");
    g_nLLTProcessCount = readConfigInteger("LLTProcess");
    g_nLogLevel = readConfigInteger("LogLevel");

    logging(LOG_DEBUG, "COMMAND [%s]\n", argv[1]);
  }
  memset(mersenne, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nMaxMersenne/64 + 1));
}

int readConfigInteger(char *pKey) {
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
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;

  memset(mersenne, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nMaxMersenne/64 + 1));

  if((fp =fopen(MERSENNE_FILE_NAME, "r")) != NULL ){
    while(!feof(fp)){
      memset(buff, 0x00, MAX_SIZE);
      size = fread(&buff, 1, MAX_SIZE, fp);
      //logging(LOG_DEBUG, "read size=[%d]\n", size);
      if(size <= 0) break;
      memcpy(mersenne+pos,buff, size);
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
  lock_leave();
}

void writeBits() {
  lock_enter();
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;
  if((fp = fopen(MERSENNE_FILE_NAME, "w")) != NULL ) {
    for(int i = 0; i < (g_nMaxMersenne/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
      memset(buff, 0x00, MAX_SIZE);
      memcpy(buff, mersenne+pos, MAX_SIZE);
      size = fwrite(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
  lock_leave();
}

void readDB() {
  lock_enter();
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;

  memset(mersenne, MERSENNE_NUMBER, sizeof(uint64_t)*(g_nMaxMersenne/64 + 1));

  if((fp =fopen(DB_FILE_NAME, "r")) != NULL ){
    while(!feof(fp)){
      memset(buff, 0x00, MAX_SIZE);
      size = fread(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      memcpy(mersenne+pos,buff, size);
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
  lock_leave();
}

void writeDB() {
  lock_enter();
  FILE *fp = NULL;
  char buff[MAX_SIZE];
  int pos = 0;
  int size = 0;
  if((fp = fopen(DB_FILE_NAME, "w")) != NULL ) {
    for(int i = 0; i < (g_nMaxMersenne/64 + 1)/(MAX_SIZE/sizeof(uint64_t)) + 1;i++) {
      memset(buff, 0x00, MAX_SIZE);
      memcpy(buff, mersenne+pos, MAX_SIZE);
      size = fwrite(&buff, 1, MAX_SIZE, fp);
      if(size <= 0) break;
      pos += (size/sizeof(uint64_t));
    }
    fclose(fp);
  }
  lock_leave();
}

int setBitOnly(uint64_t nIndex, int newBit) {
  uint64_t value = mersenne[nIndex/64];
  uint64_t residue = nIndex % 64;
  int flag = 0;
  int bit = (value >> residue) & 0x01;
  if(newBit == 0) {
    if(bit != 0) {
      mersenne[nIndex/64]-=(uint64_t)pow(2, residue);
      flag = 1;
    }
  } else {
    if(bit == 0) {
      mersenne[nIndex/64]+=(uint64_t)pow(2, residue);
      flag = 1;
    }
  }
  return flag;
}

void setBit(uint64_t nIndex, int bit) {
  if(setBitOnly(nIndex, bit)) {
    writeBits();
  }
}

int getBit(uint64_t nIndex) {
  uint64_t value = mersenne[nIndex/64];
  uint64_t residue = nIndex % 64;
  return (value >> residue) & 0x01;
}

int getKnownMersenne(int number) {
  for(int i = 0; i < KNOWN_MERSENNE; i++) {
    if(knownMersenne[i] == number) {
      return 1;
    } 
  }
  return 0;
}

int getCandidateCount() {
  int nCount = 0;
  for(int i = 0; i < g_nMaxMersenne; i++) {
    if(getBit(i) == 1 && getKnownMersenne(i) == 0) {
      nCount++; 
    }
  }
  return nCount;
}

void lock_enter() {
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
