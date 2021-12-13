/******************************************/
/*
  duplex.cpp
  by Gary P. Scavone, 2006-2007.

  This program opens a duplex stream and passes
  input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double *_sintbl = 0;
int maxfftsize = 0;
int fft(double *x, double *y, const int m);
  int ifft(double *x, double *y, const int m);
int fftr(double *x, double *y, const int m);
int ifftr(double *x, double *y, const int l);
  static int checkm(const int m);
int get_nextpow2(int n);
char *getmem(int leng, unsigned size);
double *dgetmem(int leng);
double get_process_time();

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/

//typedef signed short MY_TYPE;
//#define FORMAT RTAUDIO_SINT16

/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32
*/

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64



//! The structure for specifying input or output stream parameters.
struct dataConvTemp {
  double* repImpul;     /*pointer to impulse response*/
  unsigned int repImpulLength; /*Length of the impulse response*/
  //unsigned int samplingFrequency;    /*!< Sampling frequency */
  double* bufferInter; /*Intermediary buffer*/
  unsigned int inputBufferLength; /*L input buffer length*/
  double* conv; /* array with the convolution between impulse and input*/
};

struct dataConvFreq{
  double* repImpulFreq_x; /*pointer to impulse response already in frequency domain*/
  double* repImpulFreq_y; /*pointer to impulse response imaginary part*/
  unsigned int repImpulFreqLength; /*Length of the impulse response*/
  //unsigned int samplingFrequency;    /*!< Sampling frequency */
  double* bufferInter; /*Intermediary buffer*/
  unsigned int inputBufferLength; /*L input buffer length*/
  double* resultProduct_x; /* array with the product between the input and the impulse response in frequency domain real part*/
  double* resultProduct_y; /* array with the product between the input and the impulse response in frequency domain imaginary part*/
  double* resultIfft;/* array ifft of resultProduct*/
  double* inputFreq; /* array with each input buffer will be stored in frequency domain*/
  double* x;/* real part fft (also input of fft)*/
  double* y;/* imaginary part fft */
};




int max(int a, int b){
  if (a<= b){
    return b;
  }
  return a;
  }



///////////////////////////////
// FFT functions
int get_nextpow2(int n)
{
  int k = 1;
  while (k < n){
    k *= 2;
  }

  return k;
}

int fftr(double *x, double *y, const int m)
{
   int i, j;
   double *xp, *yp, *xq;
   double *yq;
   int mv2, n, tblsize;
   double xt, yt, *sinp, *cosp;
   double arg;

   mv2 = m / 2;

   /* separate even and odd  */
   xq = xp = x;
   yp = y;
   for (i = mv2; --i >= 0;) {
      *xp++ = *xq++;
      *yp++ = *xq++;
   }

   if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
      return (-1);


   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = M_PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }
   //printf("Debug: m=%i, maxfftsize=%i\n",m,maxfftsize);

   n = maxfftsize / m;
   sinp = _sintbl;
   cosp = _sintbl + maxfftsize / 4;

   xp = x;
   yp = y;
   xq = xp + m;
   yq = yp + m;
   *(xp + mv2) = *xp - *yp;
   *xp = *xp + *yp;
   *(yp + mv2) = *yp = 0;

   for (i = mv2, j = mv2 - 2; --i; j -= 2) {
      ++xp;
      ++yp;
      sinp += n;
      cosp += n;
      yt = *yp + *(yp + j);
      xt = *xp - *(xp + j);
      *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
      *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
   }

   xp = x + 1;
   yp = y + 1;
   xq = x + m;
   yq = y + m;

   for (i = mv2; --i;) {
      *xp++ = *(--xq);
      *yp++ = -(*(--yq));
   }

   return (0);
}

int ifftr(double *x, double *y, const int l)
{
   int i;
   double *xp, *yp;

   fftr(x, y, l);

   xp = x;
   yp = y;
   i = l;
   while (i--) {
      *xp++ /= l;
      *yp++ /= -l;
   }

   return (0);
}


static int checkm(const int m)
{
   int k;

   for (k = 4; k <= m; k <<= 1) {
      if (k == m)
         return (0);
   }
   fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);

   return (-1);
}

int fft(double *x, double *y, const int m)
{
   int j, lmx, li;
   double *xp, *yp;
   double *sinp, *cosp;
   int lf, lix, tblsize;
   int mv2, mm1;
   double t1, t2;
   double arg;
   int checkm(const int);

   /**************
   * RADIX-2 FFT *
   **************/

   if (checkm(m))
      return (-1);

   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = M_PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   lf = maxfftsize / m;
   lmx = m;

   for (;;) {
      lix = lmx;
      lmx /= 2;
      if (lmx <= 1)
         break;
      sinp = _sintbl;
      cosp = _sintbl + maxfftsize / 4;
      for (j = 0; j < lmx; j++) {
         xp = &x[j];
         yp = &y[j];
         for (li = lix; li <= m; li += lix) {
            t1 = *(xp) - *(xp + lmx);
            t2 = *(yp) - *(yp + lmx);
            *(xp) += *(xp + lmx);
            *(yp) += *(yp + lmx);
            *(xp + lmx) = *cosp * t1 + *sinp * t2;
            *(yp + lmx) = *cosp * t2 - *sinp * t1;
            xp += lix;
            yp += lix;
         }
         sinp += lf;
         cosp += lf;
      }
      lf += lf;
   }

   xp = x;
   yp = y;
   for (li = m / 2; li--; xp += 2, yp += 2) {
      t1 = *(xp) - *(xp + 1);
      t2 = *(yp) - *(yp + 1);
      *(xp) += *(xp + 1);
      *(yp) += *(yp + 1);
      *(xp + 1) = t1;
      *(yp + 1) = t2;
   }

   /***************
   * bit reversal *
   ***************/
   j = 0;
   xp = x;
   yp = y;
   mv2 = m / 2;
   mm1 = m - 1;
   for (lmx = 0; lmx < mm1; lmx++) {
      if ((li = lmx - j) < 0) {
         t1 = *(xp);
         t2 = *(yp);
         *(xp) = *(xp + li);
         *(yp) = *(yp + li);
         *(xp + li) = t1;
         *(yp + li) = t2;
      }
      li = mv2;
      while (li <= j) {
         j -= li;
         li /= 2;
      }
      j += li;
      xp = x + j;
      yp = y + j;
   }

   return (0);
}

int ifft(double *x, double *y, const int m)
{
   int i;

   if (fft(y, x, m) == -1)
      return (-1);

   for (i = m; --i >= 0; ++x, ++y) {
      *x /= m;
      *y /= m;
   }

   return (0);
}

double *dgetmem(int leng)
{
    return ( (double *)getmem(leng, sizeof(double)) );
}

char *getmem(int leng, unsigned size)
{
    char *p = NULL;

    if ((p = (char *)calloc(leng, size)) == NULL){
        fprintf(stderr, "Memory allocation error !\n");
        exit(3);
    }
    return (p);
}

double get_process_time() {
    struct rusage usage;
    if( 0 == getrusage(RUSAGE_SELF, &usage) ) {
        return (double)(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
               (double)(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1.0e6;
    }
    return 0;
}


void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset> type\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n";
  std::cout << "    type : freq or temp.\n\n";
  exit( 0 );
}

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
  dataConvTemp* d = (dataConvTemp*)data;
  int L = d->inputBufferLength;
  memcpy( outputBuffer, inputBuffer, L * sizeof(double));
 // printf("%d\n",*bytes);
  return 0;
}
size_t getFileSize(size_t blocSize,const char* filename){
  FILE* stream = fopen(filename,"r");
  if ( stream == NULL ) {
        fprintf( stderr, "Cannot open file for reading\n" );
        exit( -1 );
    }
  // Search size of file
  size_t blocCount = 0;
  double test_read;


  while(fread(&test_read,blocSize,1,stream) == 1){
    blocCount ++;
  }
  fclose(stream);
  return blocCount;
}


int readFile(dataConvTemp* data,size_t blocSize, size_t blocCount,const char* filename){
  FILE* stream = fopen(filename,"r");
  if ( stream == NULL ) {
        fprintf( stderr, "Cannot open file for reading\n" );
        exit( -1 );
  }

  printf("BlocCOunt : %d \n",blocCount);

  data->repImpulLength = blocCount;
  printf("M2 : %d \n", data->repImpulLength);
  data->repImpul = (double*)malloc(blocSize * blocCount);
  size_t x = fread(data->repImpul, blocSize, blocCount, stream) ;
  if( x != blocCount){
    
    printf("ERROR FREAD\n");
    exit(-1);
  }
  return 0;

}


/*
Callback function doing temporal convolution 
*/
int callback_conv_freq( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data ){
  
  double tic = get_process_time();


  dataConvFreq* d = (dataConvFreq*)data;
  int L = d->inputBufferLength;
  int M = d->repImpulFreqLength;
  int maxLM = get_nextpow2(max(L,M));

  memset(d->x, 0, maxLM * sizeof(double));
  memset(d->y, 0, maxLM * sizeof(double));
  memcpy(d->x,inputBuffer,L * sizeof(double));
  //printf("AAAAAAAAAAAAAAAAAAAAA : %d \n",maxLM);
  fftr(d->x,d->y,maxLM);

  for (int i = 0; i < maxLM; i ++){
    d->resultProduct_x[i] = d->x[i]  * d->repImpulFreq_x[i] - d->y[i]  * d->repImpulFreq_y[i];
    d->resultProduct_y[i] = d->x[i]  * d->repImpulFreq_y[i] + d->y[i]  * d->repImpulFreq_x[i];
  }

  ifft(d->resultProduct_x,d->resultProduct_y,maxLM);

  for(int i = 0;i< maxLM - L; i++){
    d->resultProduct_x[i] += d->bufferInter[i];
  }

  /*for(int j= 0; j< L;j++){
    ((double*)outputBuffer)[j] += d->resultProduct_x[j]; 
  }*/

  memcpy(outputBuffer,d->resultProduct_x,L * sizeof(double));
  memcpy(d->bufferInter,d->resultProduct_x + L, (maxLM - L) * sizeof(double));

  /*for (int i = 0;i<512; i ++ ){
      printf("index : %d \t valeur outputbuffer : %f \n",i,((double*)outputBuffer)[i]);
  }*/

double toc = get_process_time();
  

printf("Elapsed time: %f\n",(toc-tic));
return 0;
}

/*
Callback function doing temporal convolution 
*/
int callback_conv_temp( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data ){
  dataConvTemp* d = (dataConvTemp*)data;
  int L = d->inputBufferLength;
  int M = d->repImpulLength;
double tic = get_process_time();
  // Create an array initialy with null values
  
  for(int i = 0; i < L + M - 1; i++){
    double sum = 0;
    for(int j = 0; j < L; j++){
      if(i >= j && (i - j) < M){
        //printf("B j : %d \n",j);
        sum += ((double*)inputBuffer)[j] * d->repImpul[i - j];
      }
      //printf("B j : %d \n",j);
    }
    //printf("B\n");
    d->conv[i] = sum;
    // Add overlapp-add
    if(i <= M - 1){
      //printf("C %d  \n", i);
      d->conv[i] += d->bufferInter[i];
      //printf("C2   \n");
    }
  }
  
  // Update outputBuffer
  
  memcpy(outputBuffer,d->conv,L * sizeof(double));
/*
  for (int i = 0; i < d->inputBufferLength; i++)
  {
   
    printf("i: %d\t",i);
    printf("input: %f\t",((double*)inputBuffer)[i]);
    printf("conv: %f\t",(d->conv)[i]);
    if(i < M - 1) printf("overlapp add : %f\t", d->bufferInter[i]);
    printf("output: %f\n",((double*)outputBuffer)[i]);
  }*/

  // Update interBuffer
  
  memcpy(d->bufferInter, d->conv + L, (M - 1) * sizeof(double));
  

  /*
  for(int i = 0; i < M - 1; i++){
    printf("%f\n", conv[L + i]);
  }*/
double toc = get_process_time();
printf("Elapsed time: %f\n",(toc-tic));
  return 0;
}


int main( int argc, char *argv[] )
{
  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;
  char type[6];
  dataConvTemp bufferBytesTemp;
  dataConvFreq bufferBytesFreq;
  // Minimal command-line checking
  if (argc < 4 || argc > 8 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  strcpy(type, argv[3]);
  if ( argc > 4 )
    iDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    oDevice = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    iOffset = (unsigned int) atoi(argv[6]);
  if ( argc > 7 )
    oOffset = (unsigned int) atoi(argv[7]);
  

  printf("%s\n",type);
  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 1100;
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = iDevice;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.deviceId = oDevice;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;

  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();

  RtAudio::StreamOptions options;
  //options.flags |= RTAUDIO_NONINTERLEAVED;
  int blocSize = sizeof(double);
  char* fileName = "../../ressources_tstr_v1_1/c/impres";
  size_t blocCount = 64000 ;//getFileSize(blocSize,fileName);
  readFile(&bufferBytesTemp,blocSize,blocCount,fileName);

  if(strcmp(type,"temp") == 0){
    // TEMPOREL
    bufferBytesTemp.bufferInter = (double*) malloc(sizeof(double) * (bufferBytesTemp.repImpulLength - 1 ));
    bufferBytesTemp.inputBufferLength = bufferFrames * channels;
    bufferBytesTemp.conv = (double*) malloc(sizeof(double) * (bufferBytesTemp.inputBufferLength + bufferBytesTemp.repImpulLength - 1));
    printf("REpImpLength = %d \n", bufferBytesTemp.repImpulLength);
    printf("InputLength = %d \n", bufferBytesTemp.inputBufferLength);


    try {
      adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &callback_conv_temp, (void *)&bufferBytesTemp, &options );

    }
    catch ( RtAudioError& e ) {
      std::cout << '\n' << e.getMessage() << '\n' << std::endl;
      exit( 1 );
    }
  }
  if(strcmp(type,"freq") == 0){
        // FREQUENTIEL
    bufferBytesFreq.inputBufferLength = bufferFrames * channels;
    bufferBytesFreq.repImpulFreqLength = bufferBytesTemp.repImpulLength;

    int maxLM = get_nextpow2( max(bufferBytesFreq.inputBufferLength ,bufferBytesFreq.repImpulFreqLength));

    bufferBytesFreq.x = (double*)calloc(maxLM, maxLM * sizeof(double));
    bufferBytesFreq.y = (double*)calloc(maxLM, maxLM * sizeof(double));
    
    bufferBytesFreq.repImpulFreq_x = (double*)calloc(maxLM, maxLM  * sizeof(double));
    bufferBytesFreq.repImpulFreq_y = (double*)calloc(maxLM, maxLM  * sizeof(double));
    memcpy(bufferBytesFreq.repImpulFreq_x,bufferBytesTemp.repImpul, bufferBytesFreq.repImpulFreqLength);
    fft(bufferBytesFreq.repImpulFreq_x, bufferBytesFreq.repImpulFreq_y, maxLM);

    bufferBytesFreq.resultProduct_x = (double*)calloc(maxLM,maxLM * sizeof(double));
    bufferBytesFreq.resultProduct_y = (double*)calloc(maxLM,maxLM * sizeof(double));
    bufferBytesFreq.resultIfft = (double*)calloc(maxLM,maxLM * sizeof(double));
    
    bufferBytesFreq.bufferInter = (double*) malloc(sizeof(double) * (maxLM - bufferBytesFreq.inputBufferLength ));
    
    try {
      adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &callback_conv_freq, (void *)&bufferBytesFreq, &options );

    }
    catch ( RtAudioError& e ) {
      std::cout << '\n' << e.getMessage() << '\n' << std::endl;
      exit( 1 );
    }

  }

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  try {
    adac.startStream();

    char input;
    std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
    std::cin.get(input);

    // Stop the stream.
    adac.stopStream();
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    goto cleanup;
  }

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  return 0;
}
