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



/*
Callback function doing convolution 
*/
int callback_conv( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data ){

  dataConv* d = (dataConv*)data;
  unsigned int L = d->inputBufferLength;
  unsigned int M = d->repImpulLength;

  double conv[L + M - 1];
  for (unsigned int i = 0; i < L; i++){
    for (unsigned int j = 0; j < M; j++){
      conv[i + j] += d->repImpul[j] * ((double*)inputBuffer)[L - i - 1];
      //printf("j : %d \t valeur : %f\n",j,((double*)d->repImpul)[j]);
    } 
  }
  


  memcpy(outputBuffer,conv,(size_t)L * sizeof(double));
  for (size_t i = 0; i < d->inputBufferLength; i++)
  {
    /* code */
    printf("i: %d\t",i);
    printf("input: %f\t",((double*)inputBuffer)[i]);
    printf("conv: %f\t",((double*)conv)[i]);
    printf("output: %f\n",((double*)outputBuffer)[i]);
  }
  

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
  double* tab_tempo = (double*)malloc(sizeof(double)  * 5);
  for(int i = 0; i < 5; i++){
    tab_tempo[i] =  0;
  }
  tab_tempo[0] = 1;
  bufferBytes.bufferInter = NULL;
  bufferBytes.inputBufferLength = bufferFrames * channels;
  bufferBytes.repImpul = tab_tempo;
  bufferBytes.repImpulLength = 5; 


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
