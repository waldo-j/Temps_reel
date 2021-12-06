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
struct dataConv {
  double* repImpul;     /*pointer to impulse response*/
  unsigned int repImpulLength; /*Length of the impulse response*/
  //unsigned int samplingFrequency;    /*!< Sampling frequency */
  double* bufferInter; /*Intermediary buffer*/
  unsigned int inputBufferLength; /*L input buffer length*/
  double* conv; /* array with the convolution between impulse and input*/
};


void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  unsigned int *bytes = (unsigned int *) data;
  memcpy( outputBuffer, inputBuffer, *bytes );
  printf("%d\n",*bytes);
  return 0;
}

int readFile(dataConv data, size_t blocSize,const char* filename){
  FILE* stream = fopen(filename,"r");
  if ( stream == NULL ) {
        fprintf( stderr, "Cannot open file for reading\n" );
        exit( -1 );
    }
  // Search size of file
  int blocCount = 0;
  double test_read;


  while(fread(&test_read,blocSize,1,stream) == 1){
    blocCount ++;
  }
  fclose(stream);

  stream = fopen(filename,"r");
  if ( stream == NULL ) {
        fprintf( stderr, "Cannot open file for reading\n" );
        exit( -1 );
  }

  printf("BlocCOunt : %d \n",blocCount);

  data.repImpulLength = blocCount;
  data.repImpul = (double*)malloc(blocSize * blocCount);
  size_t x = fread(data.repImpul, blocSize, blocCount, stream) ;
  if( x != blocCount){
    printf("X = %d; BlocCount = %d ; BlocSize = %d \n",x, blocCount,blocSize);
    printf("ERROR FREAD\n");
    exit(-1);
  }





/*
  if(blocCount != fread(outputBuffer,blocSize, blocCount,stream)){
    printf("ERROR READING FILE");
    return -1;
  } */
  return 0;

}

/*
Callback function doing convolution 
*/
int callback_conv( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data ){
  dataConv* d = (dataConv*)data;
  int L = d->inputBufferLength;
  int M = d->repImpulLength;

  // Create an array initialy with null values


  for(int i = 0; i < L + M - 1; i++){
    double sum = 0;
    for(int j = 0; j < L; j++){
      if(i >= j && i - j < M){
        sum += ((double*)inputBuffer)[j] * d->repImpul[i - j];
      }
    }
    d->conv[i] = sum;
    // Add overlapp-add
    if(i <= M - 1){
      d->conv[i] += d->bufferInter[i];
    }
  }
  
  // Update outputBuffer
  memcpy(outputBuffer,d->conv,L * sizeof(double));
/*
  for (int i = 0; i < d->inputBufferLength; i++)
  {
   
    printf("i: %d\t",i);
    printf("input: %f\t",((double*)inputBuffer)[i]);
    printf("conv: %f\t",((double*)conv)[i]);
    if(i < M - 1) printf("overlapp add : %f\t", d->bufferInter[i]);
    printf("output: %f\n",((double*)outputBuffer)[i]);
  }*/

  // Update interBuffer
  memcpy(d->bufferInter, d->conv + L, (M - 1) * sizeof(double));
  /*
  for(int i = 0; i < M - 1; i++){
    printf("%f\n", conv[L + i]);
  }*/

  return 0;
}

int main( int argc, char *argv[] )
{
  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;
  dataConv bufferBytes;
  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 512;
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
  int repImpulLength = 5;
  double* tab_tempo = (double*)malloc(sizeof(double)  * repImpulLength);
  for(int i = 0; i < repImpulLength; i++){
    tab_tempo[i] =  0;
  }
  tab_tempo[0] = 1;

  int blocSize = sizeof(double);
  char* fileName = "../../ressources_tstr_v1_1/c/impres";
  readFile(bufferBytes,blocSize,fileName);

  bufferBytes.bufferInter = (double*) malloc(sizeof(double) * (repImpulLength - 1 ));
  bufferBytes.inputBufferLength = bufferFrames * channels;
  bufferBytes.conv = (double*) malloc(sizeof(double) * (bufferBytes.inputBufferLength + bufferBytes.repImpulLength - 1));

  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &callback_conv, (void *)&bufferBytes, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  //bufferBytes = bufferFrames * channels * sizeof( MY_TYPE );

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
